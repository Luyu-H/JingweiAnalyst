import re
import os
import ast
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

import yaml
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib
matplotlib.use("Agg")  # headless: figures are saved, never shown
import matplotlib.pyplot as plt

from langchain_core.messages import SystemMessage, HumanMessage

from src.tools.data_openers import dataset_notes, list_datasets
from src.tools.registry import TOOLS, auto_import_tools
from src.models.prompt import PromptTemplates
from src.models.llm_providers import LLM_KEYS, build_llm

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

# Every setting a run depends on - the model and its endpoint, the loop limits,
# where results land - lives in this file, so a run is changed by editing YAML
# rather than the source.
CONFIG_PATH = REPO_ROOT / "configs" / "agent_with_tools.yaml"


def load_config(path: Optional[Path] = None) -> Dict[str, Any]:
    """Read the runtime configuration, configs/agent_with_tools.yaml by default."""
    path = Path(path) if path is not None else CONFIG_PATH
    if not path.exists():
        raise FileNotFoundError(f"Agent configuration not found: {path}")
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f) or {}


def resolve_path(value: str) -> Path:
    """A configured path is absolute, or relative to the repo root - never to the cwd."""
    path = Path(value).expanduser()
    return path if path.is_absolute() else REPO_ROOT / path


def slugify(text: str) -> str:
    """
    A model id as a directory name.

    Vendors write model ids with characters a path cannot take: an
    openai_compatible endpoint may serve "Qwen/Qwen3-235B", and a version may
    carry a colon. Everything outside letters, digits, dot, dash and underscore
    becomes a single dash.
    """
    return re.sub(r"[^A-Za-z0-9._-]+", "-", text).strip("-") or "model"


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
    plt.figure(number).savefig(out / ("figure_%d.png" % index), dpi=__DPI__, bbox_inches="tight")
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

    def write_rerun_harness(self, datasets: List[str], figure_dpi: int = 150) -> None:
        alias = f"namespace[\"selected_data\"] = {datasets[0]!r}" if len(datasets) == 1 else ""
        script = (RERUN_HARNESS
                  .replace("__REPO__", str(REPO_ROOT))
                  .replace("__DATASETS__", repr(datasets))
                  .replace("__ALIAS__", alias)
                  .replace("__DPI__", str(figure_dpi)))
        self.write("code/run.py", script)


class ConversationSession:
    """Manage a multi-turn conversation session."""
    def __init__(self, session_id: str, max_turns: int = 10, context_turns: int = 3):
        self.session_id = session_id
        self.max_turns = max_turns
        self.context_turns = context_turns
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
        for i, turn in enumerate(self.history[-self.context_turns:], 1):  # only the most recent turns
            summary += f"Q{i}: {turn['user_query']}\n"
            if 'solution' in turn:
                solution_match = re.search(r"<solution>(.*?)</solution>", turn['solution'], re.S)
                solution_text = solution_match.group(1).strip() if solution_match else turn['solution']
                summary += f"A{i}: {solution_text}\n"
        return summary


class AgentWithTools:
    """An agent that can select and call provided tools to answer scientific questions with multi-turn conversation support."""
    
    def __init__(self, tools: Dict[str, Dict], data_list: List[Dict[str, str]],
                 api_key: Optional[str] = None, config: Optional[Dict[str, Any]] = None,
                 config_path: Optional[Path] = None, **overrides):
        """
        Args:
            - tools, data_list: the registry and dataset catalogue the agent works with
            - api_key: the provider key; read from the environment variable named in
              the configuration when left out
            - config / config_path: an already-loaded configuration, or another YAML
              file to read instead of configs/agent_with_tools.yaml
            - overrides: single settings from the 'agent' or 'llm' section, by their
              configured name, for a caller that needs to deviate for one run
        """
        self.config = config if config is not None else load_config(config_path)
        self.config_path = Path(config_path) if config_path is not None else CONFIG_PATH
        self.config_from_caller = config is not None
        self.config_overrides = dict(overrides)
        llm_cfg = dict(self.config.get("llm") or {})
        agent_cfg = dict(self.config.get("agent") or {})
        observation_cfg = dict(agent_cfg.pop("observation", None) or {})
        for key, value in overrides.items():
            if key in LLM_KEYS:
                llm_cfg[key] = value
            elif key in observation_cfg or key in ("max_array_elements", "max_value_chars"):
                observation_cfg[key] = value
            else:
                agent_cfg[key] = value

        self.tools_meta = tools
        self.data_list = data_list
        self.data_index = {data["name"]: data for data in data_list}
        self.dataset_notes = dataset_notes()

        # Which vendor, which endpoint and which environment variable holds the key
        # is decided in src/models/llm_providers.py from the 'llm' section; the agent
        # itself only ever sees a chat model it can invoke.
        self.llm, self.llm_settings = build_llm(llm_cfg, api_key, hint=str(CONFIG_PATH))
        self.provider = self.llm_settings["provider"]
        self.model_name = self.llm_settings["model"]
        self.temperature = self.llm_settings["temperature"]
        self.base_url = self.llm_settings["base_url"]

        self.result_dir = resolve_path(agent_cfg.get("result_dir", "./results"))
        self.result_dir.mkdir(parents=True, exist_ok=True)

        self.sessions: Dict[str, ConversationSession] = {}
        self.max_turns_per_session = agent_cfg.get("max_turns_per_session", 10)
        self.max_code_iterations = agent_cfg.get("max_code_iterations", 3)
        self.language = agent_cfg.get("language", "auto")
        self.context_history_turns = agent_cfg.get("context_history_turns", 3)
        self.figure_dpi = agent_cfg.get("figure_dpi", 150)
        self.max_array_elements = observation_cfg.get("max_array_elements", 50)
        self.max_value_chars = observation_cfg.get("max_value_chars", 800)

        # The settings as they ended up, overrides included. Written to each
        # turn's bundle so a reader never has to reconstruct them from the YAML.
        self.agent_settings = {
            "result_dir": str(self.result_dir),
            "max_turns_per_session": self.max_turns_per_session,
            "max_code_iterations": self.max_code_iterations,
            "language": self.language,
            "context_history_turns": self.context_history_turns,
            "figure_dpi": self.figure_dpi,
            "observation": {"max_array_elements": self.max_array_elements,
                            "max_value_chars": self.max_value_chars},
        }

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
                        "scikit-learn", "matplotlib", "pyproj", "langchain-core",
                        "langchain-openai", "openai", "langchain-anthropic", "anthropic",
                        "langchain-google-genai"):
            try:
                versions[package] = importlib.metadata.version(package)
            except Exception:
                continue

        config_digests = {}
        for name in ("datasets.yaml", "default_paths.yaml", "agent_with_tools.yaml"):
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
            # Everything that reached the client, so a later reader knows which
            # vendor, endpoint and sampling settings produced these numbers.
            "llm": {key: value for key, value in self.llm_settings.items() if key != "kind"},
        }

    def _configuration_snapshot(self) -> Dict[str, Any]:
        """
        Everything the agent was configured with, as a record of its own.

        The manifest says what the turn did; this says under which settings. It is
        kept in a separate file because it is the same for every turn of a session
        and because it is what a reader copies to repeat the run: the resolved
        model settings, every agent-side limit, the overrides a caller passed, and
        the configuration file itself - content as well as digest, so the record
        can still be read on a machine that does not have that file.
        """
        config_file: Dict[str, Any] = {"path": str(self.config_path)}
        if self.config_path.exists():
            config_file["sha256"] = hashlib.sha256(self.config_path.read_bytes()).hexdigest()
        if self.config_from_caller:
            config_file["note"] = ("the configuration was supplied by the caller; the file above is "
                                   "only where it would have been read from")
        config_file["content"] = self.config

        return {
            "recorded_at": datetime.now().isoformat(),
            # What actually reached the client, after the provider's defaults and
            # any override were applied - not what the file says.
            "llm": self.llm_settings,
            "agent": self.agent_settings,
            "overrides_from_caller": self.config_overrides,
            "config_file": config_file,
            "always_available_tools": list(ALWAYS_AVAILABLE_TOOLS),
            "registered_tools": sorted(self.tools_meta),
            "registered_datasets": sorted(self.data_index),
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
        
        # Every attempt starts from a clean canvas. A failed attempt can leave a
        # half-drawn figure open, and _save_visualizations walks whatever is open
        # when it runs, so without this the discarded attempt's figure would be
        # saved beside the good one and reported as a result.
        plt.close("all")

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

    def _summarize_value(self, value: Any, max_chars: Optional[int] = None) -> Any:
        """
        Render one result variable for the solution step.

        The loaders hand back lazy dask-backed DataArrays, whose repr shows the task
        graph and no values at all, so a plain str() would leave the final answer
        with no numbers in it. Arrays are therefore materialised and reduced to a
        compact numeric summary, with small ones listed in full.
        """
        if max_chars is None:
            max_chars = self.max_value_chars

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
            summary.update(self._describe_array(value.values, self.max_array_elements))
            return summary

        if isinstance(value, np.ndarray):
            return {"type": "ndarray", "shape": list(value.shape),
                    **self._describe_array(value, self.max_array_elements)}
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

    @staticmethod
    def _loop_local_names(code: str) -> set:
        """
        Names that only a loop assigned, which therefore hold its last iteration.

        `for m in models: dm = compute_difference_metrics(...)` leaves `dm` bound to
        the last model alone. Reported as an ordinary result it reads as a property
        of the analysis, and the final answer then states one model's number as if
        it held for all of them. A name the code also assigns outside the loop is a
        different matter - `metrics = {}` filled in the loop keeps every model - so
        those are kept.
        """
        try:
            tree = ast.parse(code)
        except SyntaxError:
            return set()

        loop_names, outer_names = set(), set()

        def collect(node, inside_loop: bool) -> None:
            for child in ast.iter_child_nodes(node):
                target = loop_names if inside_loop else outer_names
                if isinstance(child, ast.Name) and isinstance(child.ctx, ast.Store):
                    target.add(child.id)
                elif isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                    # A function's own locals never reach the module namespace, so
                    # only the name it binds matters.
                    target.add(child.name)
                elif isinstance(child, (ast.Lambda, ast.ListComp, ast.SetComp,
                                        ast.DictComp, ast.GeneratorExp)):
                    # A lambda's arguments and a comprehension's loop variable have
                    # their own scope and never bind in the namespace the code runs
                    # in, so they must not count as an assignment outside a loop.
                    continue
                else:
                    collect(child, inside_loop or isinstance(child, (ast.For, ast.AsyncFor, ast.While)))

        collect(tree, False)
        return loop_names - outer_names

    # Names a generated script uses when it parks a caught exception instead of
    # raising it. Only string-like values are inspected, so a numeric 'error'
    # (a standard error, a fit residual) is never mistaken for a failure.
    _ERROR_KEYS = ("error", "errors", "exception", "exceptions", "err",
                   "traceback", "failure", "failures", "failed")
    _ERROR_TEXT = re.compile(
        r"Traceback \(most recent call last\)"
        r"|^\s*[A-Za-z_][A-Za-z0-9_.]*(?:Error|Exception|Warning)\s*\("
        r"|^\s*[A-Za-z_][A-Za-z0-9_.]*(?:Error|Exception)\s*:\s*\S"
    )

    @classmethod
    def _looks_like_error_text(cls, value: Any) -> bool:
        return isinstance(value, str) and bool(cls._ERROR_TEXT.search(value))

    @classmethod
    def _detect_swallowed_errors(cls, output: Dict[str, Any], max_depth: int = 4,
                                 limit: int = 12) -> List[str]:
        """
        Find failures the generated code caught and stored instead of raising.

        A `try/except` that records the exception and carries on leaves the
        interpreter exit status clean, so the run reads as a success and the
        repair budget is never spent - while part of the answer is quietly
        missing. This walks the result namespace for the two shapes that
        actually occur: an exception object kept in a variable, and an error
        string parked under a key such as 'error'.
        """
        found: List[str] = []

        def walk(node: Any, path: str, depth: int) -> None:
            if len(found) >= limit or depth > max_depth:
                return
            if isinstance(node, BaseException):
                found.append(f"{path} holds {type(node).__name__}: {node}")
                return
            if isinstance(node, dict):
                for key, value in node.items():
                    if len(found) >= limit:
                        return
                    where = f"{path}[{key!r}]"
                    if (isinstance(key, str) and key.strip().lower() in cls._ERROR_KEYS
                            and cls._looks_like_error_text(value)):
                        found.append(f"{where} = {value.strip().splitlines()[-1][:200]}")
                    else:
                        walk(value, where, depth + 1)
                return
            if isinstance(node, (list, tuple, set)):
                for i, value in enumerate(list(node)[:50]):
                    walk(value, f"{path}[{i}]", depth + 1)
                return

        for name, value in output.items():
            if name.startswith("_") or callable(value) or isinstance(value, types.ModuleType):
                continue
            if cls._looks_like_error_text(value) and name.strip().lower() in cls._ERROR_KEYS:
                found.append(f"{name} = {value.strip().splitlines()[-1][:200]}")
            else:
                walk(value, name, 1)
        return found[:limit]

    def _summarize_outputs(self, output: Dict[str, Any], code: str = "") -> Dict[str, Any]:
        """Summarize every result variable the executed code left behind."""
        loop_locals = self._loop_local_names(code) if code else set()
        return {
            name: self._summarize_value(value)
            for name, value in output.items()
            if not name.startswith("_")
            and name not in loop_locals
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
                "agent_configuration": "agent_config.json",
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
                plt.figure(num).savefig(vis_path, dpi=self.figure_dpi, bbox_inches='tight')
                paths.append(str(vis_path))
                print(f"[System] Visualization saved to {vis_path}")
        except Exception as e:
            print(f"[Warning] Failed to save visualization: {e}")
        finally:
            plt.close('all')
        return paths
    
    def _save_session(self, session: ConversationSession):
        """
        Save session history beside the turns it summarises.

        The session id already carries the case, the model and the time, so the
        history goes inside that session's own directory rather than back at the
        top of result_dir under a second timestamp.
        """
        fname = self.result_dir / session.session_id / "session_history.json"
        fname.parent.mkdir(parents=True, exist_ok=True)
        
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
                           bundle: Optional["RunBundle"] = None, upstream_code: str = "",
                           available_variables: Optional[Dict[str, Any]] = None
                           ) -> Tuple[Dict[str, Any], List[Dict[str, Any]]]:
        """
        Execute code with retries upon failure, regenerating code as needed.

        Args:
            - upstream_code: code that already ran and whose variables this block
              builds on; shown to the repair step so a broken visualization is not
              repaired by recomputing the answer a second time
            - available_variables: what those variables actually hold, as the same
              summary the final answer is written from

        Returns:
            - Dict: Final execution result dict
            - List: List of all execution attempts with code and results
        """
        iteration_history = []
        current_code = initial_code
        # A run that only swallowed errors is still worth more than nothing. It is
        # kept here so that if every later repair raises outright, the partial
        # answer is returned rather than lost to the retry it triggered.
        best_partial: Optional[Dict[str, Any]] = None

        for iteration in range(self.max_code_iterations):
            print(f"\n--- Iteration {iteration + 1}/{self.max_code_iterations} ({code_type}) ---")
            
            exec_result = self._safe_execute_tool(current_code, context, code_type)

            # A caught-and-stored exception leaves the exit status clean, so without
            # this the run counts as a success and the repair budget goes unspent
            # while part of the answer is missing. Detected before the attempt is
            # recorded, so the audit bundle does not file it as a clean success.
            swallowed = (self._detect_swallowed_errors(exec_result["output"])
                         if exec_result["success"] else [])
            retries_left = iteration < self.max_code_iterations - 1
            if swallowed:
                print(f"[Warning] {len(swallowed)} failure(s) were caught and stored by the "
                      f"generated code rather than raised:")
                for item in swallowed:
                    print(f"          - {item}")

            if exec_result["success"] and swallowed and retries_left:
                # Rewritten as a failure so the existing repair path handles it.
                exec_result["swallowed_errors"] = swallowed
                exec_result["summary"] = self._summarize_outputs(exec_result["output"],
                                                                 current_code)
                best_partial = exec_result
                detail = "\n".join(f"- {item}" for item in swallowed)
                exec_result = {
                    "success": False,
                    "error": ("The code caught these failures and stored them instead of "
                              "letting them raise, so part of the analysis is missing:\n"
                              + detail),
                    "error_type": "SwallowedError",
                    "traceback": detail,
                    "hint": ("Remove the try/except around the failing work, or fix the cause, "
                             "so every requested item is actually computed. Do not record an "
                             "error string and continue."),
                    "registered_tools": list(self.tools_meta.keys()),
                }
            if bundle is not None:
                bundle.record_attempt(code_type, iteration + 1, current_code, exec_result)

            iteration_record = {
                "iteration": iteration + 1,
                "code": current_code,
                "success": exec_result["success"]
            }
            if swallowed:
                iteration_record["swallowed_errors"] = swallowed

            if exec_result["success"]:
                if swallowed:
                    # Last iteration: keep the partial result rather than discarding
                    # it, but carry the failures into the observation so neither the
                    # answer nor a reader can miss them.
                    exec_result["swallowed_errors"] = swallowed
                # Summarised once here and reused by the caller, so lazy arrays are
                # materialised a single time rather than on every read.
                exec_result["summary"] = self._summarize_outputs(exec_result["output"], current_code)
                iteration_record["output"] = exec_result["summary"]
                omitted = sorted(self._loop_local_names(current_code)
                                 & set(exec_result["output"]))
                if omitted:
                    # Recorded, not hidden: an auditor can see which names were held
                    # back and why, without the model reading one model's number as
                    # if it were everyone's.
                    iteration_record["omitted_loop_variables"] = omitted
                    print(f"[Info] Held back {len(omitted)} loop-local variable(s) from the "
                          f"observation: {', '.join(omitted)}")
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
                    language, current_code, exec_result, user_query, specs_for_fix, code_type,
                    upstream_code=upstream_code, available_variables=available_variables
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

        if not exec_result.get("success") and best_partial is not None:
            print("[System] Every repair failed outright; returning the earlier partial "
                  "result, whose caught failures are carried into the observation.")
            return best_partial, iteration_history

        return exec_result, iteration_history
    
    def create_session(self, session_id: Optional[str] = None,
                       label: Optional[str] = None) -> str:
        """
        Create a new conversation session.

        The session id is the directory every turn of the session is written to, so
        it carries the model that produced those turns: a run of the same case under
        two models is then two directories side by side, distinguishable without
        opening either. Pass `label` for what the session is (a case id, say) and
        the id becomes "<label>_<model>_<timestamp>"; pass `session_id` to name the
        directory yourself.
        """
        if session_id is None:
            session_id = "_".join([
                slugify(label) if label else "session",
                slugify(self.model_name),
                datetime.now().strftime("%Y%m%d_%H%M%S"),
            ])
        
        self.sessions[session_id] = ConversationSession(session_id, self.max_turns_per_session,
                                                       self.context_history_turns)
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
        # The configuration can pin the answer's language; 'auto' reads it off the question.
        language = self.language if self.language in ("zh", "en") else self._detect_language(user_query)
        print(f"[System] Language: {'Chinese' if language == 'zh' else 'English'}"
              f"{' (detected)' if self.language not in ('zh', 'en') else ' (from configuration)'}")

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
        bundle.write_json("agent_config.json", self._configuration_snapshot())

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

        bundle.write_rerun_harness(selected_datasets, self.figure_dpi)

        compute_result, compute_iterations = self._execute_with_retry(
            compute_code,
            context=code_context,
            code_type="compute",
            full_specs=full_specs,
            user_query=user_query,
            language=language,
            bundle=bundle,
            # Nothing ran before this block; the repair step is told only what the
            # sandbox pre-binds, so it does not invent a dataset id of its own.
            available_variables={name: self._summarize_value(value)
                                 for name, value in code_context.items()}
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
            computation_part = {
                "type": "computation",
                "status": "success",
                "results": compute_obs
            }
            # Reported, never dropped: the code ran to completion but caught these
            # on the way, so some requested quantity is missing. The answer step is
            # told, and told to say so.
            if compute_result.get("swallowed_errors"):
                computation_part["status"] = "success_with_caught_failures"
                computation_part["caught_failures"] = compute_result["swallowed_errors"]
                computation_part["message"] = (
                    "The code caught and stored these failures instead of raising them, so "
                    "the corresponding results are MISSING. Say so explicitly in the answer, "
                    "name what could not be computed, and do not present the remaining items "
                    "as if the whole question had been answered."
                )
            observation_parts.append(computation_part)

            # Step 4: Extract and execute the visualization code (if any)
            visualize_match = re.search(r"<visualize>(.*?)</visualize>", exec_response.content, re.S)
            if visualize_match:
                visualize_code = visualize_match.group(1).strip()
                
                if visualize_code:
                    vis_context = dict(code_context)
                    vis_context.update(compute_result["output"])
                    
                    # A visualization failure is nearly always a name or a shape:
                    # the block is handed variables it did not create. So the repair
                    # step gets the code that created them and a summary of each.
                    vis_variables = {name: self._summarize_value(value)
                                     for name, value in code_context.items()}
                    vis_variables.update(compute_result.get("summary", {}))

                    vis_result, vis_iterations = self._execute_with_retry(
                        visualize_code,
                        context=vis_context,
                        code_type="visualization",
                        full_specs=full_specs,
                        user_query=user_query,
                        language=language,
                        bundle=bundle,
                        upstream_code=final_compute,
                        available_variables=vis_variables
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

        # Step 5: Generate final answer. The question, the plan and the code that
        # actually ran all go in, so that the methodology is reported from the code
        # rather than guessed from the names of the result variables.
        # A run that caught failures still produced numbers, so the answer is written
        # from them - the observation carries the instruction to declare what is missing.
        has_computation_success = observation_parts[0]["status"] in (
            "success", "success_with_caught_failures")
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

                # The case's own id when it has one, so the directory says which
                # case, which model and when - T1_deepseek-reasoner_20260922_193101.
                label = item.get("case_id") or Path(json_file).stem
                session_id = agent.create_session(label=label)

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

    # Model, limits, output directory and the batch below all come from
    # configs/agent_with_tools.yaml; the API key comes from the environment
    # variable that file names, so no key is ever committed.
    config = load_config()
    agent = AgentWithTools(TOOLS, data_list, config=config)

    batch = config.get("batch_run") or {}
    folder = resolve_path(batch.get("test_folder", "cases"))
    repeats = int(batch.get("repeats", 1))
    print(f"[System] Running {folder} x{repeats} with {agent.model_name} "
          f"({agent.provider}; results -> {agent.result_dir})")

    for _ in range(repeats):
        results = run_tests_from_folder(str(folder), agent)

    print("Finished all tests.")
    
