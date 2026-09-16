import re
import os
import sys
import json
import types
import hashlib
import platform
import subprocess
import traceback
import importlib.metadata
from datetime import datetime
from typing import Any, Dict, List, Optional, Tuple
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib
matplotlib.use("Agg")  # headless: figures are saved, never shown
import matplotlib.pyplot as plt

from langchain_openai import ChatOpenAI
from langchain_core.messages import SystemMessage, HumanMessage

from src.tools.data_openers import dataset_notes, list_datasets
from src.tools.registry import TOOLS, auto_import_tools
from src.models.prompt import PromptTemplates

auto_import_tools()

# Specs always handed to the code-generation step, whether or not the planning
# step listed them: the prompt tells the model to load through these loaders and
# to weight averages with compute_basic_stat, so it must be able to see their
# signatures even when the plan forgot to name them.
ALWAYS_AVAILABLE_TOOLS = (
    "load_regular_data",
    "load_climatology_data",
    "load_data_with_geolocator",
    "compute_basic_stat",
    "transform_unit",
)


REPO_ROOT = Path(__file__).resolve().parents[2]

# Harness written into every bundle so the recorded code can be re-executed
# without the agent. It rebuilds the namespace the agent executes code in: every
# registered tool by name, the scientific stack, and the turn's dataset ids.
RERUN_HARNESS = '''
"""Re-run this turn's computation exactly as the agent executed it.

    cd __REPO__ && python <bundle>/code/run.py

The agent runs generated code with every registered tool already in scope and
the selected dataset ids bound to `selected_datasets`, so this script rebuilds
that namespace instead of importing the code as a module. Figures are written to
figures/rerun/ so they can be compared against the originals in figures/.
"""
import sys
from pathlib import Path

BUNDLE = Path(__file__).resolve().parent.parent
REPO = Path(r"__REPO__")
sys.path.insert(0, str(REPO))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

# The same rcParams the agent's sandbox sets: different font metrics shift the
# figures by a pixel or two and an exact reproduction then looks like a mismatch.
plt.rcParams["font.sans-serif"] = ["Noto Sans CJK SC"]
plt.rcParams["axes.unicode_minus"] = False

from src.tools.registry import TOOLS, auto_import_tools
auto_import_tools()

namespace = {name: meta["func"] for name, meta in TOOLS.items()}
namespace.update({"np": np, "pd": pd, "xr": xr, "plt": plt})
namespace["selected_datasets"] = __DATASETS__
__ALIAS__

with xr.set_options(keep_attrs=True):
    exec(compile((BUNDLE / "code" / "compute.py").read_text(), "compute.py", "exec"), namespace)
    visualize = BUNDLE / "code" / "visualize.py"
    if visualize.exists():
        exec(compile(visualize.read_text(), "visualize.py", "exec"), namespace)

out = BUNDLE / "figures" / "rerun"
out.mkdir(parents=True, exist_ok=True)
for index, number in enumerate(plt.get_fignums(), 1):
    plt.figure(number).savefig(out / ("figure_%d.png" % index), dpi=150, bbox_inches="tight")
print("re-ran %s; %d figure(s) written to %s" % (BUNDLE.name, len(plt.get_fignums()), out))
'''


class RunBundle:
    """
    A self-contained, auditable record of one turn, written as the turn happens.

    Everything a reviewer needs to judge or reuse the result lands on disk here:
    the plan, the exact prompts the model saw, every code attempt with its
    traceback, the final runnable code, the observation, the answer, the figures,
    and a manifest pinning the code version, the dataset registry and the library
    versions. Written step by step rather than at the end, so a crash mid-turn
    still leaves behind the evidence gathered so far.
    """

    def __init__(self, root: Path):
        self.root = root
        for sub_dir in ("prompts", "code/attempts", "figures"):
            (self.root / sub_dir).mkdir(parents=True, exist_ok=True)
        self._prompt_count = 0

    def write(self, relative: str, text: str) -> Path:
        path = self.root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")
        return path

    def write_json(self, relative: str, payload: Any) -> Path:
        return self.write(relative, json.dumps(payload, indent=2, ensure_ascii=False, default=str))

    def record_exchange(self, name: str, system: str, human: str, response: str,
                        seconds: float) -> Path:
        """Store one LLM exchange verbatim: what was sent, and what came back."""
        self._prompt_count += 1
        body = (
            f"# exchange {self._prompt_count}: {name}\n"
            f"# elapsed: {seconds:.1f} s\n"
            f"\n===== SYSTEM =====\n{system}\n"
            f"\n===== USER =====\n{human}\n"
            f"\n===== RESPONSE =====\n{response}\n"
        )
        return self.write(f"prompts/{self._prompt_count:02d}_{name}.txt", body)

    def record_attempt(self, code_type: str, iteration: int, code: str,
                       result: Dict[str, Any]) -> None:
        """Store one execution attempt, with its traceback when it failed."""
        stem = f"code/attempts/{code_type}_{iteration:02d}"
        self.write(f"{stem}.py", code if code.endswith("\n") else code + "\n")
        if not result.get("success"):
            self.write(
                f"{stem}.error.txt",
                f"{result.get('error_type', 'Error')}: {result.get('error', '')}\n"
                f"{result.get('hint', '')}\n\n{result.get('traceback', '')}\n"
            )

    def write_rerun_harness(self, datasets: List[str]) -> None:
        alias = f"namespace[\"selected_data\"] = {datasets[0]!r}" if len(datasets) == 1 else ""
        script = (RERUN_HARNESS
                  .replace("__REPO__", str(REPO_ROOT))
                  .replace("__DATASETS__", repr(datasets))
                  .replace("__ALIAS__", alias))
        self.write("code/run.py", script)


class ConversationSession:
    """Manage a multi-turn conversation session."""
    def __init__(self, session_id: str, max_turns: int = 10):
        self.session_id = session_id
        self.max_turns = max_turns
        self.turn_count = 0
        self.history = []  # Record of all turns
        
    def add_turn(self, turn_record: Dict[str, Any]):
        """Add a turn record to the session history."""
        self.turn_count += 1
        turn_record["turn_number"] = self.turn_count
        turn_record["timestamp"] = datetime.now().isoformat()
        self.history.append(turn_record)
    
    def is_expired(self) -> bool:
        return self.turn_count >= self.max_turns
    
    def get_context_summary(self) -> str:
        """Generate a brief summary of the conversation history for context."""
        if not self.history:
            return ""
        
        summary = f"Previous conversation (Turn {self.turn_count}):\n"
        for i, turn in enumerate(self.history[-3:], 1):  # only summarize last 3 turns
            summary += f"Q{i}: {turn['user_query']}\n"
            if 'solution' in turn:
                solution_match = re.search(r"<solution>(.*?)</solution>", turn['solution'], re.S)
                solution_text = solution_match.group(1).strip() if solution_match else turn['solution']
                summary += f"A{i}: {solution_text}\n"
        return summary


class AgentWithTools:
    """An agent that can select and call provided tools to answer scientific questions with multi-turn conversation support."""
    
    def __init__(self, tools: Dict[str, Dict], data_list: List[Dict[str, str]],
                 api_key: str, model_name: str = "deepseek-reasoner", temperature: float = 0.0, 
                 result_dir: str = "./results", max_turns_per_session: int = 10, max_code_iterations: int = 3):
        self.tools_meta = tools
        self.data_list = data_list
        self.data_index = {data["name"]: data for data in data_list}
        self.dataset_notes = dataset_notes()
        self.model_name = model_name
        self.temperature = temperature
        self.base_url = "https://api.deepseek.com"
        self.llm = ChatOpenAI(model=model_name, temperature=temperature, api_key=api_key, base_url=self.base_url)
        self.result_dir = Path(result_dir)
        self.result_dir.mkdir(parents=True, exist_ok=True)
        
        self.sessions: Dict[str, ConversationSession] = {}
        self.max_turns_per_session = max_turns_per_session
        self.max_code_iterations = max_code_iterations
    
    def _environment_snapshot(self) -> Dict[str, Any]:
        """
        What a later reader needs to tell whether a rerun stands on the same ground:
        the code version, the dataset registry's fingerprint, the library versions
        and the model settings. A result whose numbers cannot be tied to these is
        not reproducible, only repeatable by luck.
        """
        def git(*args) -> Optional[str]:
            try:
                done = subprocess.run(["git", "-C", str(REPO_ROOT), *args],
                                      capture_output=True, text=True, timeout=10)
                return done.stdout.strip() or None
            except Exception:
                return None

        versions = {}
        for package in ("xarray", "numpy", "pandas", "dask", "scipy", "gsw", "netCDF4",
                        "scikit-learn", "matplotlib", "pyproj", "langchain-openai", "openai"):
            try:
                versions[package] = importlib.metadata.version(package)
            except Exception:
                continue

        config_digests = {}
        for name in ("datasets.yaml", "default_paths.yaml"):
            path = REPO_ROOT / "configs" / name
            if path.exists():
                config_digests[name] = hashlib.sha256(path.read_bytes()).hexdigest()

        return {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "packages": versions,
            "repo_root": str(REPO_ROOT),
            "git_commit": git("rev-parse", "HEAD"),
            "git_describe": git("describe", "--always", "--dirty"),
            "git_uncommitted_changes": bool(git("status", "--porcelain")),
            "config_sha256": config_digests,
            "llm": {"model": self.model_name, "temperature": self.temperature,
                    "base_url": self.base_url},
        }

    def _ask(self, bundle: Optional["RunBundle"], name: str, system: str, human: str):
        """One LLM exchange, timed and recorded in the bundle."""
        started = datetime.now()
        response = self.llm.invoke([SystemMessage(content=system), HumanMessage(content=human)])
        seconds = (datetime.now() - started).total_seconds()
        if bundle is not None:
            bundle.record_exchange(name, system, human, response.content, seconds)
        return response, seconds

    def _specs_for(self, names) -> Dict[str, Any]:
        """Full specification of each registered tool in names, in order."""
        specs = {}
        for name in names:
            if name in self.tools_meta and name not in specs:
                meta = self.tools_meta[name]
                specs[name] = {
                    "description": meta["description"],
                    "args": meta.get("args", {}),
                    "returns": meta.get("returns", "")
                }
        return specs

    def _specs_used_by(self, code: str) -> Dict[str, Any]:
        """
        Specs of every registered tool the code actually calls.

        The sandbox exposes all registered tools, but the code-generation step is
        only shown the ones the plan named. When generated code reaches for another
        one, the fix step needs its signature to repair the call.
        """
        called = [name for name in self.tools_meta
                  if re.search(rf"\b{re.escape(name)}\s*\(", code)]
        return self._specs_for(called)

    def _detect_language(self, text: str) -> str:
        """Detect if the input text is Chinese or English.
        Returns 'zh' for Chinese, 'en' for English."""
        chinese_char_count = len(re.findall(r'[\u4e00-\u9fff]', text))
        english_char_count = len(re.findall(r'[a-zA-Z]', text))

        if chinese_char_count > english_char_count:
            return 'zh'
        else:
            return 'en'

    def _safe_execute_tool(self, code: str, context: dict = None, code_type: str = "compute") -> Dict[str, Any]:
        """Safely execute generated code in a sandboxed environment with registered tools."""
        sandbox_globals = {"__builtins__": __builtins__}
        
        registered_count = 0
        for name, meta in self.tools_meta.items():
            if not isinstance(meta, dict):
                print(f"[Warning] Tool '{name}' meta is not a dict: {type(meta)}")
                continue
                
            func = meta.get("func")
            if func is None:
                print(f"[Warning] Tool '{name}' has no 'func' key in meta")
                continue
            if not callable(func):
                print(f"[Warning] Tool '{name}' func is not callable: {type(func)}")
                continue
            
            sandbox_globals[name] = func
            registered_count += 1
        
        # The scientific stack the generated code is told it may use. Every tool
        # takes and returns xr.DataArray, so xarray and pandas belong here.
        sandbox_globals.update({"np": np, "pd": pd, "xr": xr, "plt": plt})
        plt.rcParams['font.sans-serif'] = ['Noto Sans CJK SC']
        plt.rcParams['axes.unicode_minus'] = False
        
        print(f"[Info] Successfully registered {registered_count} tools to sandbox")

        if code_type == "compute":
            print(f"[Info] Executing computation code with {registered_count} registered tools")
        else:
            print(f"[Info] Executing visualization code")

        if context:
            sandbox_globals.update(context)
            
        # One namespace for globals and locals. With two dicts, a function defined
        # in the generated code cannot see names bound at top level (those land in
        # locals, which a function body never consults) and dies with NameError.
        predefined = set(sandbox_globals)
        try:
            # xarray drops attrs on reductions and arithmetic by default, which
            # strips the source_dataset/grid_type tags the weighted statistics and
            # derived-quantity tools need to find the right grid. Generated code
            # routinely does `da.mean(dim="time")` before a weighted mean, so the
            # tags are kept for the duration of the generated code.
            with xr.set_options(keep_attrs=True):
                exec(code, sandbox_globals)
            output = {k: v for k, v in sandbox_globals.items() if k not in predefined}
            return {"success": True, "output": output}
        except Exception as e:
            tb = traceback.format_exc()
            error_info = {
                "success": False,
                "error": str(e),
                "error_type": type(e).__name__,
                "traceback": tb,
                "registered_tools": list(self.tools_meta.keys())
            }
            
            if isinstance(e, NameError):
                error_info["hint"] = "Function or variable not found. Check spelling and availability."
            elif isinstance(e, TypeError):
                error_info["hint"] = "Type error. Check function arguments and data types."
            elif isinstance(e, KeyError):
                error_info["hint"] = "Key not found in data structure. Check data keys."
            elif isinstance(e, IndexError):
                error_info["hint"] = "Index out of range. Check array dimensions."
            
            return error_info
    
    @staticmethod
    def _describe_array(values, max_full: int = 50) -> Dict[str, Any]:
        """Numbers from an array: the values themselves if few, statistics if many."""
        array = np.asarray(values)
        described: Dict[str, Any] = {"size": int(array.size)}
        if array.size == 0:
            return described

        if not np.issubdtype(array.dtype, np.number):
            described["values"] = [str(v) for v in array.ravel()[:max_full]]
            return described

        finite = array[np.isfinite(array)]
        described["n_valid"] = int(finite.size)
        if finite.size:
            described.update({
                "min": float(np.min(finite)),
                "max": float(np.max(finite)),
                "mean": float(np.mean(finite)),
                "std": float(np.std(finite)),
            })
        if array.size <= max_full:
            described["values"] = (np.round(array, 6).tolist()
                                   if np.issubdtype(array.dtype, np.floating) else array.tolist())
        return described

    def _summarize_value(self, value: Any, max_chars: int = 800) -> Any:
        """
        Render one result variable for the solution step.

        The loaders hand back lazy dask-backed DataArrays, whose repr shows the task
        graph and no values at all, so a plain str() would leave the final answer
        with no numbers in it. Arrays are therefore materialised and reduced to a
        compact numeric summary, with small ones listed in full.
        """
        if isinstance(value, (xr.DataArray, xr.Dataset)):
            try:
                value = value.compute()
            except Exception as e:
                return f"<{type(value).__name__} could not be computed: {e}>"

        if isinstance(value, xr.Dataset):
            return {name: self._summarize_value(value[name], max_chars) for name in value.data_vars}

        if isinstance(value, xr.DataArray):
            summary = {
                "type": "DataArray",
                "dims": {str(k): int(v) for k, v in value.sizes.items()},
                "units": value.attrs.get("units"),
                "long_name": value.attrs.get("long_name") or value.attrs.get("standard_name"),
                "source_dataset": value.attrs.get("source_dataset"),
            }
            summary = {k: v for k, v in summary.items() if v is not None}
            summary.update(self._describe_array(value.values))
            return summary

        if isinstance(value, np.ndarray):
            return {"type": "ndarray", "shape": list(value.shape), **self._describe_array(value)}
        if isinstance(value, np.generic):
            return value.item()
        if isinstance(value, (pd.DataFrame, pd.Series)):
            return str(value)[:max_chars]
        if value is None or isinstance(value, (bool, int, float, str)):
            return value
        if isinstance(value, dict):
            return {str(k): self._summarize_value(v, max_chars) for k, v in value.items()}
        if isinstance(value, (list, tuple)):
            if len(value) <= 20:
                return [self._summarize_value(v, max_chars) for v in value]
            return {"type": type(value).__name__, "length": len(value),
                    "head": [self._summarize_value(v, max_chars) for v in value[:5]]}
        return str(value)[:max_chars]

    def _summarize_outputs(self, output: Dict[str, Any]) -> Dict[str, Any]:
        """Summarize every result variable the executed code left behind."""
        return {
            name: self._summarize_value(value)
            for name, value in output.items()
            if not name.startswith("_")
            and not callable(value)
            and not isinstance(value, types.ModuleType)
        }

    def _manifest(self, record: Dict[str, Any], started_at: datetime,
                  compute_iterations: List[Dict[str, Any]],
                  vis_iterations: List[Dict[str, Any]]) -> Dict[str, Any]:
        """The bundle's index: what was asked, what answered it, and on what footing."""
        names = record.get("selected_data", {}).get("name", [])
        if isinstance(names, str):
            names = [names]

        def steps(iterations):
            return [{"iteration": item["iteration"], "success": item["success"],
                     "error_type": item.get("error_type")} for item in iterations]

        return {
            "session_id": record.get("session_id"),
            "turn_number": record.get("turn_number"),
            "user_query": record.get("user_query"),
            "language": record.get("language"),
            "started_at": started_at.isoformat(),
            "finished_at": datetime.now().isoformat(),
            "elapsed_seconds": round((datetime.now() - started_at).total_seconds(), 1),
            "status": "error" if record.get("error") else "ok",
            "error": record.get("error"),
            "selected_datasets": [
                {"name": name,
                 "root": self.data_index[name]["path"],
                 "description": self.data_index[name].get("description", "")}
                for name in names if name in self.data_index
            ],
            "selected_tools": record.get("selected_tools", []),
            "code_iterations": {"compute": steps(compute_iterations),
                                "visualization": steps(vis_iterations)},
            "figures": [os.path.basename(path) for path in record.get("visualization_paths", [])],
            "files": {
                "query": "query.txt",
                "plan": "plan.md",
                "prompts_sent_to_the_model": "prompts/",
                "final_compute_code": "code/compute.py",
                "final_visualization_code": "code/visualize.py",
                "every_attempt_with_traceback": "code/attempts/",
                "rerun_harness": "code/run.py",
                "observation": "observation.json",
                "answer": "solution.md",
                "figures": "figures/",
            },
            "environment": self._environment_snapshot(),
        }

    def _save_visualizations(self, turn_id: str, out_dir: Optional[Path] = None) -> List[str]:
        """
        Save every matplotlib figure the visualization code left open.

        Questions that want one panel per basin or per model produce several
        figures, and savefig writes only the current one, so they are walked here.
        """
        paths = []
        try:
            destination = Path(out_dir) if out_dir is not None else self.result_dir
            destination.mkdir(parents=True, exist_ok=True)
            fignums = plt.get_fignums()
            for i, num in enumerate(fignums, 1):
                suffix = "" if len(fignums) == 1 else f"_{i}"
                vis_path = destination / f"visualization_{turn_id}{suffix}.png"
                plt.figure(num).savefig(vis_path, dpi=150, bbox_inches='tight')
                paths.append(str(vis_path))
                print(f"[System] Visualization saved to {vis_path}")
        except Exception as e:
            print(f"[Warning] Failed to save visualization: {e}")
        finally:
            plt.close('all')
        return paths
    
    def _save_session(self, session: ConversationSession):
        """Save session history to a JSON file."""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        fname = self.result_dir / f"session_{session.session_id}_{timestamp}.json"
        
        session_data = {
            "session_id": session.session_id,
            "total_turns": session.turn_count,
            "max_turns": session.max_turns,
            "history": session.history
        }
        
        with open(fname, "w", encoding="utf-8") as f:
            json.dump(session_data, f, indent=2, ensure_ascii=False)
        print(f"[System] Session saved to {fname}")
    
    def _execute_with_retry(self, initial_code: str, context: dict, code_type: str,
                           full_specs: Dict[str, Any], user_query: str, language: str = 'en',
                           bundle: Optional["RunBundle"] = None) -> Tuple[Dict[str, Any], List[Dict[str, Any]]]:
        """
        Execute code with retries upon failure, regenerating code as needed.
        
        Returns:
            - Dict: Final execution result dict
            - List: List of all execution attempts with code and results
        """
        iteration_history = []
        current_code = initial_code

        for iteration in range(self.max_code_iterations):
            print(f"\n--- Iteration {iteration + 1}/{self.max_code_iterations} ({code_type}) ---")
            
            exec_result = self._safe_execute_tool(current_code, context, code_type)
            if bundle is not None:
                bundle.record_attempt(code_type, iteration + 1, current_code, exec_result)

            iteration_record = {
                "iteration": iteration + 1,
                "code": current_code,
                "success": exec_result["success"]
            }
            
            if exec_result["success"]:
                # Summarised once here and reused by the caller, so lazy arrays are
                # materialised a single time rather than on every read.
                exec_result["summary"] = self._summarize_outputs(exec_result["output"])
                iteration_record["output"] = exec_result["summary"]
                iteration_history.append(iteration_record)
                print(f"[Success] {code_type.capitalize()} code executed successfully!")
                return exec_result, iteration_history
            
            # If failed, record error info. The traceback is kept too: without it
            # an auditor cannot tell which call inside the generated code failed.
            iteration_record["error"] = exec_result["error"]
            iteration_record["error_type"] = exec_result["error_type"]
            iteration_record["traceback"] = exec_result.get("traceback")
            iteration_history.append(iteration_record)
            
            print(f"[Error] Execution failed: {exec_result['error']}")
            
            if iteration < self.max_code_iterations - 1:
                print("[System] Requesting code fix from LLM...")
                
                # Include the spec of any tool the failed code reached for beyond
                # the ones the plan selected, otherwise the fix is blind to it.
                specs_for_fix = dict(full_specs)
                specs_for_fix.update(self._specs_used_by(current_code))

                fix_prompt = PromptTemplates.get_fix_prompt(
                    language, current_code, exec_result, user_query, specs_for_fix, code_type
                )
                system_content = PromptTemplates.get_fix_system_message(language)
                
                fix_response, _ = self._ask(
                    bundle, f"fix_{code_type}_iter{iteration + 2}", system_content, fix_prompt
                )
                
                code_match = re.search(r"<execute>(.*?)</execute>", fix_response.content, re.S)
                if code_match:
                    current_code = code_match.group(1).strip()
                    print(f"[System] Code fixed, retrying...")
                else:
                    print("[Warning] LLM did not return fixed code in proper format")
                    break
            else:
                print(f"[System] Max iterations reached for {code_type}")
        
        return exec_result, iteration_history
    
    def create_session(self, session_id: Optional[str] = None) -> str:
        """Create a new conversation session."""
        if session_id is None:
            session_id = f"session_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        
        self.sessions[session_id] = ConversationSession(session_id, self.max_turns_per_session)
        print(f"[System] Created new session: {session_id}")
        return session_id
    
    def end_session(self, session_id: str):
        """End a conversation session and save its history."""
        if session_id in self.sessions:
            self._save_session(self.sessions[session_id])
            del self.sessions[session_id]
            print(f"[System] Session {session_id} ended and saved")
    
    def run(self, user_query: str, session_id: Optional[str] = None) -> Dict[str, Any]:
        """
        Run the agent to answer a user query within a conversation session.
        
        Args:
            - user_query, str: User's question
            - session_id: Optional session ID for multi-turn conversation
            
        Returns:
            Dict, the result record including plan, selected data/tools, generated code, observation, and solution.
        """
        language = self._detect_language(user_query)
        print(f"[System] Detected language: {'Chinese' if language == 'zh' else 'English'}")

        if session_id is None:
            session_id = self.create_session()
        elif session_id not in self.sessions:
            print(f"[Warning] Session {session_id} not found, creating new session")
            session_id = self.create_session(session_id)
        
        session = self.sessions[session_id]
        
        if session.is_expired():
            print(f"[Warning] Session {session_id} has reached max turns, creating new session")
            self.end_session(session_id)
            session_id = self.create_session()
            session = self.sessions[session_id]
        
        # The bundle is opened before the first model call, so that even a turn
        # that crashes halfway leaves its prompts and attempts on disk.
        turn_number = session.turn_count + 1
        turn_id = f"{session_id}_turn{turn_number}"
        bundle = RunBundle(self.result_dir / session_id / f"turn{turn_number:02d}")
        started_at = datetime.now()

        result_record = {
            "session_id": session_id,
            "user_query": user_query,
            "language": language,
            "turn_number": turn_number,
            "bundle_path": str(bundle.root)
        }
        bundle.write("query.txt", user_query + "\n")

        context_summary = session.get_context_summary()
        
        # The planner only has to choose; it is shown one line per tool. The full
        # entry — arguments, return contract, caveats — goes to the code step for
        # the tools the plan selected, which is where it can actually be used.
        tool_overview = {
            k: (v.get("summary") or v["description"].strip().split("\n")[0])
            for k, v in self.tools_meta.items()
        }
        
        # Datasets are offered by registry id only: an id is what the loaders take,
        # and showing the root path alongside it invites a path in the answer.
        data_overview = {
            data["name"]: {"description": data.get("description", "")}
            for data in self.data_list
        }

        # Step 1: Planning, including selecting data and tools
        system_prompt_1 = PromptTemplates.get_planning_prompt(
            language, tool_overview, data_overview, context_summary
        )
        planning_input = PromptTemplates.get_planning_input(language, user_query)

        planning_response, planning_seconds = self._ask(
            bundle, "planning", system_prompt_1, planning_input
        )
        print("\n=== Plan & Selection ===\n", planning_response.content)
        result_record["plan"] = planning_response.content

        plan_match = re.search(r"<plan>(.*?)</plan>", planning_response.content, re.S)
        plan_text = plan_match.group(1).strip() if plan_match else ""
        bundle.write("plan.md", f"# Plan\n\n{plan_text or planning_response.content}\n")

        # Extract the selected datasets. A dataset is a registry id (woa23,
        # cmip6_historical, ...) holding many variables and, for CMIP6, several
        # models; the id is what the loaders take, so the id is what is carried
        # forward rather than a file path.
        selected_data_match = re.search(r"<selected_data>(.*?)</selected_data>", planning_response.content, re.S)
        requested = [d.strip() for d in selected_data_match.group(1).split(",")] if selected_data_match else []
        requested = [d for d in requested if d]

        selected_datasets = [name for name in requested if name in self.data_index]
        unknown = [name for name in requested if name not in self.data_index]
        if unknown:
            print(f"[Warning] Planning step named unregistered datasets: {unknown}")

        if not selected_datasets:
            # The turn stops rather than substituting another dataset: against a
            # registry of unlike products, a substitution quietly answers a
            # different question than the one that was asked.
            message = (f"No registered dataset was selected (requested: {requested or 'none'}; "
                       f"available: {sorted(self.data_index)})")
            print(f"[Error] {message}")
            result_record["error"] = message
            bundle.write_json("manifest.json", self._manifest(result_record, started_at, [], []))
            session.add_turn(result_record)
            return result_record

        # Only the description: the code-generation step addresses datasets by id,
        # so handing it the root path would contradict the prompt's own rule.
        dataset_specs = {
            name: {"description": self.data_index[name].get("description", "")}
            for name in selected_datasets
        }

        print(f"[System] Selected datasets: {selected_datasets}")
        result_record["selected_data"] = {
            "name": selected_datasets,
            "path": [self.data_index[name]["path"] for name in selected_datasets]
        }

        selected_tool_match = re.search(r"<selected_tools>(.*?)</selected_tools>", planning_response.content, re.S)
        selected_tool_names = [t.strip() for t in selected_tool_match.group(1).split(",")] if selected_tool_match else []
        result_record["selected_tools"] = selected_tool_names

        # Step 2: Generate executable code
        full_specs = self._specs_for(list(selected_tool_names) + list(ALWAYS_AVAILABLE_TOOLS))

        system_prompt_2 = PromptTemplates.get_code_generation_prompt(
            language, full_specs, user_query, dataset_specs, self.dataset_notes, plan=plan_text
        )
        exec_input = PromptTemplates.get_code_generation_input(language)
        
        exec_response, _ = self._ask(bundle, "code_generation", system_prompt_2, exec_input)
        print("\n=== Generated Code ===\n", exec_response.content)
        result_record["generated_code"] = exec_response.content

        # Step 3: Extract and execute the computation code
        compute_match = re.search(r"<compute>(.*?)</compute>", exec_response.content, re.S)
        if not compute_match:
            print("[Error] No <compute> code block found.")
            result_record["error"] = "No computation code found"
            bundle.write("code/raw_response.txt", exec_response.content)
            bundle.write_json("manifest.json", self._manifest(result_record, started_at, [], []))
            session.add_turn(result_record)
            return result_record
        
        compute_code = compute_match.group(1).strip()

        # The generated code addresses datasets by id. `selected_data` is a
        # convenience alias for the single-dataset case.
        code_context = {"selected_datasets": selected_datasets}
        if len(selected_datasets) == 1:
            code_context["selected_data"] = selected_datasets[0]

        bundle.write_rerun_harness(selected_datasets)

        compute_result, compute_iterations = self._execute_with_retry(
            compute_code,
            context=code_context,
            code_type="compute",
            full_specs=full_specs,
            user_query=user_query,
            language=language,
            bundle=bundle
        )
        # The code that actually ran, not the first draft: this is what a reader
        # must audit and what code/run.py re-executes.
        final_compute = compute_iterations[-1]["code"] if compute_iterations else compute_code
        bundle.write("code/compute.py", final_compute + "\n")

        result_record["compute_iterations"] = compute_iterations

        observation_parts = []
        
        if compute_result["success"]:
            print("\n=== Computation Observation ===")
            compute_obs = compute_result.get("summary", {})
            print(json.dumps(compute_obs, indent=2))
            observation_parts.append({
                "type": "computation",
                "status": "success",
                "results": compute_obs
            })

            # Step 4: Extract and execute the visualization code (if any)
            visualize_match = re.search(r"<visualize>(.*?)</visualize>", exec_response.content, re.S)
            if visualize_match:
                visualize_code = visualize_match.group(1).strip()
                
                if visualize_code:
                    vis_context = dict(code_context)
                    vis_context.update(compute_result["output"])
                    
                    vis_result, vis_iterations = self._execute_with_retry(
                        visualize_code,
                        context=vis_context,
                        code_type="visualization",
                        full_specs=full_specs,
                        user_query=user_query,
                        language=language,
                        bundle=bundle
                    )

                    result_record["visualization_iterations"] = vis_iterations
                    bundle.write("code/visualize.py",
                                 (vis_iterations[-1]["code"] if vis_iterations else visualize_code) + "\n")
                    
                    if vis_result["success"]:
                        print("\n=== Visualization Observation ===")
                        print("Visualization code executed successfully")
                        
                        vis_paths = self._save_visualizations(turn_id, bundle.root / "figures")
                        if vis_paths:
                            result_record["visualization_paths"] = vis_paths
                            observation_parts.append({
                                "type": "visualization",
                                "status": "success",
                                "paths": vis_paths
                            })
                        else:
                            observation_parts.append({
                                "type": "visualization",
                                "status": "no_figure",
                                "message": "Code executed but no figure was generated"
                            })
                    else:
                        print("\n=== Visualization Failed ===")
                        print(f"Error: {vis_result['error']}")
                        observation_parts.append({
                            "type": "visualization",
                            "status": "failed",
                            "error": vis_result['error'],
                            "message": "Visualization failed but computation succeeded"
                        })
                else:
                    print("[Info] No visualization code provided")
                    observation_parts.append({
                        "type": "visualization",
                        "status": "not_applicable",
                        "message": "No visualization was generated for this query"
                    })
            else:
                observation_parts.append({
                    "type": "visualization",
                    "status": "not_provided",
                    "message": "No visualization block found in generated code"
                })
        else:
            print("\n=== Computation Failed ===")
            print(compute_result["traceback"])
            observation_parts.append({
                "type": "computation",
                "status": "failed",
                "error": compute_result['error'],
                "iterations": len(compute_iterations)
            })

        observation_text = json.dumps(observation_parts, indent=2, ensure_ascii=False)
        result_record["observation"] = observation_text
        bundle.write_json("observation.json", observation_parts)

        # Step 4: Generate final answer. The question, the plan and the code that
        # actually ran all go in, so that the methodology is reported from the code
        # rather than guessed from the names of the result variables.
        has_computation_success = observation_parts[0]["status"] == "success"
        has_visualization = any(p.get("type") == "visualization" and p.get("status") == "success" 
                               for p in observation_parts)

        executed_code = final_compute
        if result_record.get("visualization_iterations"):
            executed_code += "\n\n# --- visualization ---\n" + \
                result_record["visualization_iterations"][-1]["code"]

        solution_prompt = PromptTemplates.get_solution_prompt(
            language, observation_text, has_computation_success, has_visualization,
            user_query=user_query, plan=plan_text, executed_code=executed_code
        )
        final_input = PromptTemplates.get_solution_input(language)

        final_response, _ = self._ask(bundle, "solution", solution_prompt, final_input)
        print("\n=== Final Solution ===\n", final_response.content)
        result_record["solution"] = final_response.content

        solution_match = re.search(r"<solution>(.*?)</solution>", final_response.content, re.S)
        bundle.write("solution.md",
                     (solution_match.group(1).strip() if solution_match else final_response.content) + "\n")

        session.add_turn(result_record)
        bundle.write_json("manifest.json", self._manifest(
            result_record, started_at, compute_iterations,
            result_record.get("visualization_iterations", [])
        ))
        print(f"[System] Audit bundle written to {bundle.root}")

        return result_record
    

if __name__ == "__main__":
    # Built from configs/datasets.yaml so this catalogue cannot drift from what
    # the loader tools actually read.
    data_list = list_datasets()


    def run_tests_from_folder(folder_path, agent):
        """
        Read all JSON files inside folder_path.
        Each JSON file contains an array of test items.
        For each item: create a new agent session and run the query.
        """

        json_files = [
            os.path.join(folder_path, f)
            for f in os.listdir(folder_path)
            if f.endswith(".json")
        ]

        print(f"Found {len(json_files)} test files.")

        all_results = []   # collect everything

        for json_file in json_files:
            print(f"\n=== Loading test file: {json_file} ===")
            with open(json_file, "r") as f:
                test_items = json.load(f)

            if not isinstance(test_items, list):
                raise ValueError(f"JSON file {json_file} must contain a list/array.")

            for idx, item in enumerate(test_items):
                query = item.get("query")

                if query is None:
                    print(f"Skipping item {idx} because no 'query' field.")
                    continue

                print(f"\n--- Running test case {idx+1}/{len(test_items)} ---")
                print(f"Query: {query}")

                session_id = agent.create_session()

                try:
                    result = agent.run(query, session_id=session_id)
                    all_results.append({
                        "file": json_file,
                        "index": idx,
                        "query": query,
                        "result": result
                    })
                except Exception as e:
                    print(f"Error running test case {idx}: {e}")
                finally:
                    agent.end_session(session_id)

        return all_results

    # Read from the environment rather than living in the source, so the key is
    # never committed: export DEEPSEEK_API_KEY before running.
    api_key = os.environ.get("DEEPSEEK_API_KEY", "")
    if not api_key:
        raise SystemExit("DEEPSEEK_API_KEY is not set in the environment.")

    agent = AgentWithTools(TOOLS, data_list, api_key, model_name="deepseek-reasoner",
                           max_turns_per_session=10, max_code_iterations=3)
    
    folder = os.path.join(os.path.dirname(__file__), "..", "..", "data", "web_applications")
    for _ in range(10):
        results = run_tests_from_folder(folder, agent)

    print("Finished all tests.")
    
