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
import dask
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib
matplotlib.use("Agg")  # headless: figures are saved, never shown
import matplotlib.pyplot as plt

from langchain_core.messages import SystemMessage, HumanMessage

from src.tools.data_openers import dataset_notes, list_datasets
from src.tools.registry import TOOLS, auto_import_tools
from src.docs.knowledge_retrieval import KnowledgeBase, KnowledgeRetriever
from src.models.prompt import PromptTemplates
from src.models.llm_providers import (LLM_KEYS, build_chat_model, build_llm, resolve_provider,
                                      resolve_settings)

auto_import_tools()

# Specs always handed to the code-generation step, whether or not the planning step listed them.
ALWAYS_AVAILABLE_TOOLS = (
    "load_regular_data",
    "load_climatology_data",
    "load_data_with_geolocator",
    "compute_basic_stat",
    "transform_unit",
)


REPO_ROOT = Path(__file__).resolve().parents[2]

CONFIG_PATH = REPO_ROOT / "configs" / "agent_full.yaml"


def config_chain(path: Optional[Path] = None) -> List[Path]:
    """
    The configuration file and every file it `extends`, the base first.

    Each arm of the ablation has its own file that extends another and holds only
    what differs, so a setting shared by every arm is changed in one place.
    """
    path = (Path(path) if path is not None else CONFIG_PATH).resolve()
    chain: List[Path] = []
    while path is not None:
        if not path.exists():
            raise FileNotFoundError(f"Agent configuration not found: {path}")
        if path in chain:
            raise ValueError(f"Configuration files extend each other in a cycle: {path}")
        chain.insert(0, path)
        with open(path, "r", encoding="utf-8") as f:
            parent = (yaml.safe_load(f) or {}).get("extends")
        path = (path.parent / parent).resolve() if parent else None
    return chain


def _merge(base: Dict[str, Any], override: Dict[str, Any]) -> Dict[str, Any]:
    """Mappings merge key by key; any other value replaces the one it overrides."""
    merged = dict(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = _merge(merged[key], value)
        else:
            merged[key] = value
    return merged


def load_config(path: Optional[Path] = None) -> Dict[str, Any]:
    """
    Read the runtime configuration, configs/agent_full.yaml by default.
    """
    config: Dict[str, Any] = {}
    for file in config_chain(path):
        with open(file, "r", encoding="utf-8") as f:
            layer = yaml.safe_load(f) or {}
        layer.pop("extends", None)
        config = _merge(config, layer)
    return config


def config_path_from_argv(default: Path) -> Path:
    """The --config of a batch run's command line, the arm's default file without one."""
    import argparse
    parser = argparse.ArgumentParser(description="Run the case batch under one configuration.")
    parser.add_argument("--config", type=Path, default=default,
                        help=f"configuration file (default: {default.relative_to(REPO_ROOT)})")
    path = parser.parse_args().config
    return path if path.is_absolute() else Path.cwd() / path


def check_arm(config: Dict[str, Any], allowed: Tuple[str, ...], path: Path) -> str:
    """
    The configuration's `arm`, refused when this agent class does not implement it.

    Keeps a no-tools configuration from being run by the tool agent (or the other
    way round), which would file the results under an arm they do not belong to.
    """
    arm = config.get("arm")
    if arm not in allowed:
        raise ValueError(f"{path} declares arm {arm!r}; this entry point runs {', '.join(allowed)}. "
                         f"See configs/README.md for which module runs which arm.")
    return arm


def resolve_path(value: str) -> Path:
    """A configured path is absolute, or relative to the repo root - never to the cwd."""
    path = Path(value).expanduser()
    return path if path.is_absolute() else REPO_ROOT / path


def slugify(text: str) -> str:
    """
    A model id as a directory name.
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

    def write_rerun_harness(self, datasets: List[str], figure_dpi: int = 150,
                            template: str = RERUN_HARNESS) -> None:
        alias = f"namespace[\"selected_data\"] = {datasets[0]!r}" if len(datasets) == 1 else ""
        script = (template
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

    # The prompts of planning, code generation, repair and the answer. A class
    # attribute so that a variant of the agent (src/models/ablation_tools) can swap
    # them while every other step stays exactly this agent's.
    prompts = PromptTemplates

    def __init__(self, tools: Dict[str, Dict], data_list: List[Dict[str, str]],
                 api_key: Optional[str] = None, config: Optional[Dict[str, Any]] = None,
                 config_path: Optional[Path] = None, **overrides):
        """
        Args:
            - tools, data_list: the registry and dataset catalogue the agent works with
            - api_key: the provider key; read from the environment variable named in
              the configuration when left out
            - config / config_path: an already-loaded configuration, or another YAML
              file to read instead of configs/agent_full.yaml
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
            elif key in observation_cfg or key in ("max_array_elements", "max_value_chars",
                                                   "materialize_max_mb", "materialize_total_mb"):
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
        self.materialize_max_mb = observation_cfg.get("materialize_max_mb", 2048)
        self.materialize_total_mb = observation_cfg.get("materialize_total_mb", 16384)

        self._setup_knowledge(dict(self.config.get("knowledge") or {}), llm_cfg, api_key)

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
                            "max_value_chars": self.max_value_chars,
                            "materialize_max_mb": self.materialize_max_mb,
                            "materialize_total_mb": self.materialize_total_mb},
        }

    def _setup_knowledge(self, knowledge_cfg: Dict[str, Any], llm_cfg: Dict[str, Any],
                         api_key: Optional[str]) -> None:
        """
        Load the knowledge base and the model that selects from it (the 'knowledge'
        section of the configuration).
        """
        self.knowledge_enabled = bool(knowledge_cfg.get("enabled", False))
        self.knowledge_base: Optional[KnowledgeBase] = None
        self.knowledge_retriever: Optional[KnowledgeRetriever] = None
        self.selector_llm, self.selector_llm_settings = self.llm, self.llm_settings
        prefilter = dict(knowledge_cfg.get("prefilter") or {})
        selector = dict(knowledge_cfg.get("selector") or {})
        condense = dict(knowledge_cfg.get("condense") or {})
        injection = dict(knowledge_cfg.get("injection") or {})
        plan_check = dict(knowledge_cfg.get("plan_check") or {})

        self.knowledge_settings = {
            "enabled": self.knowledge_enabled,
            "root": str(resolve_path(knowledge_cfg.get("root", "docs/knowledge"))),
            "validate_on_start": bool(knowledge_cfg.get("validate_on_start", True)),
            "max_candidates": prefilter.get("max_candidates"),
            "severity_weights": prefilter.get("severity_weights"),
            "backtick_only_variables": list(prefilter.get("backtick_only_variables") or []),
            "dataset_aliases": dict(prefilter.get("dataset_aliases") or {}),
            "selector_mode": selector.get("mode", "llm"),
            "max_selected": int(selector.get("max_selected", 8)),
            "fallback_top_k": int(selector.get("fallback_top_k", 4)),
            "retries": int(selector.get("retries", 1)),
            "condense_enabled": bool(condense.get("enabled", True)),
            "source_sections": condense.get("source_sections") or {
                "must": ["key_points", "implications", "common_mistakes"],
                "should": ["key_points", "implications", "common_mistakes"],
                "info": ["key_points", "implications"]},
            "max_notes_per_fragment": int(condense.get("max_notes_per_fragment", 2)),
            "condense_max_chars": int(condense.get("max_chars", 3600)),
            "always_allowed_numbers": list(condense.get("always_allowed_numbers") or []),
            "max_chars": int(injection.get("max_chars", 16000)),
            "sections": injection.get("sections") or {},
            "on_uncited_must": plan_check.get("on_uncited_must", "warn"),
        }
        if not self.knowledge_enabled:
            return
        if self.knowledge_settings["selector_mode"] not in ("llm", "lexical"):
            raise ValueError(f"knowledge.selector.mode must be 'llm' or 'lexical', not "
                             f"{self.knowledge_settings['selector_mode']!r} ({CONFIG_PATH})")

        self.knowledge_base = KnowledgeBase(
            Path(self.knowledge_settings["root"]),
            severity_weights=self.knowledge_settings["severity_weights"],
            backtick_only_variables=self.knowledge_settings["backtick_only_variables"],
            dataset_aliases=self.knowledge_settings["dataset_aliases"],
            always_allowed_numbers=self.knowledge_settings["always_allowed_numbers"],
        )
        if self.knowledge_settings["validate_on_start"]:
            problems = self.knowledge_base.check()
            if problems:
                raise ValueError(
                    "The knowledge base failed its checks (python -m src.docs.knowledge_index "
                    "--check):\n" + "\n".join(f"  - {p}" for p in problems))
        self.knowledge_retriever = KnowledgeRetriever(self.knowledge_base, self.knowledge_settings)
        print(f"[System] Knowledge base: {len(self.knowledge_base.fragments)} fragments from "
              f"{self.knowledge_settings['root']}")

        override = selector.get("llm")
        if self.knowledge_settings["selector_mode"] == "llm" and override:
            self.selector_llm, self.selector_llm_settings = self._build_selector_llm(
                llm_cfg, dict(override), api_key)
            print(f"[System] Knowledge selector: {self.selector_llm_settings['model']} "
                  f"({self.selector_llm_settings['provider']})")

    def _build_selector_llm(self, llm_cfg: Dict[str, Any], override: Dict[str, Any],
                            api_key: Optional[str]):
        """
        The selector's model: the 'llm' section with the override's settings on top.
        """
        merged = dict(llm_cfg)
        if "provider" in override and resolve_provider(override["provider"])[0] != self.provider:
            for key in ("model_name", "base_url", "api_key_env", "extra_params"):
                if key not in override:
                    merged.pop(key, None)
            api_key = None
        merged.update(override)
        settings = resolve_settings(merged)
        if settings == self.llm_settings:
            return self.llm, self.llm_settings
        return build_chat_model(settings, api_key, hint=str(CONFIG_PATH)), settings

    def _retrieve_knowledge(self, user_query: str, context_summary: str, language: str,
                            bundle: Optional["RunBundle"]) -> Tuple[str, Dict[str, Any]]:
        """
        The knowledge for this question, condensed for the planner
        (src/docs/knowledge_retrieval.py does the work).

        Returns the text for the planning prompt ("" when nothing applies, which
        leaves that prompt as it is without retrieval) and the record written to the
        bundle's knowledge.json. The selection and condensing calls go through _ask,
        so they are recorded under prompts/ like every other exchange.
        """
        if not self.knowledge_enabled:
            return "", {"enabled": False}

        def ask(name: str, system: str, human: str):
            response, seconds = self._ask(bundle, name, system, human, llm=self.selector_llm)
            return response.content, seconds

        text, record = self.knowledge_retriever.retrieve(user_query, context_summary, language, ask)
        if self.knowledge_settings["selector_mode"] == "llm":
            record["selector_model"] = self.selector_llm_settings["model"]

        for stage in ("selection", "condensing_run"):
            reason = (record.get(stage) or {}).get("fallback_reason")
            if reason:
                print(f"[Warning] Knowledge {stage.split('_')[0]} fell back: {reason}")
        by_status: Dict[str, List[str]] = {}
        for entry in record.get("selected", []):
            by_status.setdefault(entry["status"], []).append(entry["id"])
        print(f"[System] Knowledge for planning ({record['injected_chars']} chars): "
              f"{', '.join(record['injected_ids']) or 'none'}")
        for status, label in (("not_applicable", "judged not applicable"),
                              ("summary_fallback", "given as summary only"),
                              ("dropped", "dropped for lack of room")):
            if by_status.get(status):
                print(f"[Info] Knowledge {label}: {', '.join(by_status[status])}")
        return text, record

    def _environment_snapshot(self) -> Dict[str, Any]:
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

        # The data registry, and every file of this arm's configuration chain.
        config_files = [REPO_ROOT / "configs" / name for name in ("datasets.yaml", "default_paths.yaml")]
        if self.config_path.exists():
            config_files += config_chain(self.config_path)
        config_digests = {path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                          for path in config_files if path.exists()}

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
            # The fragments the planner could have been shown, so a result can be
            # tied to the version of the knowledge base it was planned against.
            "knowledge": self.knowledge_base.digest() if self.knowledge_base else None,
        }

    def _configuration_snapshot(self) -> Dict[str, Any]:
        """
        Everything the agent was configured with, as a record of its own.
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
            # Which arm of the ablation this is (configs/README.md), as its file names it.
            "arm": self.config.get("arm"),
            # What actually reached the client, after the provider's defaults and
            # any override were applied - not what the file says.
            "llm": self.llm_settings,
            "agent": self.agent_settings,
            "knowledge": {**self.knowledge_settings,
                          "selector_llm": (self.selector_llm_settings
                                           if self.knowledge_enabled
                                           and self.knowledge_settings["selector_mode"] == "llm"
                                           else None)},
            "overrides_from_caller": self.config_overrides,
            "config_file": config_file,
            "always_available_tools": list(ALWAYS_AVAILABLE_TOOLS),
            "registered_tools": sorted(self.tools_meta),
            "registered_datasets": sorted(self.data_index),
        }

    def _ask(self, bundle: Optional["RunBundle"], name: str, system: str, human: str, llm=None):
        """One LLM exchange, timed and recorded in the bundle; the agent's own model unless llm is given."""
        started = datetime.now()
        response = (llm or self.llm).invoke([SystemMessage(content=system), HumanMessage(content=human)])
        if not isinstance(response.content, str):
            # A model that thinks before answering (claude-sonnet-5-5 by default) returns a
            # list of content blocks, thinking first; every caller parses the answer text.
            response.content = str(response.text)
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

    def _code_reference(self, selected_tools: List[str],
                        selected_datasets: List[str]) -> Dict[str, Any]:
        """
        What the code-generation and repair steps are given to program against: the
        full specs of the tools the plan selected, plus the ones always handed over.
        """
        return self._specs_for(list(selected_tools) + list(ALWAYS_AVAILABLE_TOOLS))

    def _write_rerun_harness(self, bundle: "RunBundle", datasets: List[str]) -> None:
        """The bundle's code/run.py, which re-executes the turn's code without the agent."""
        bundle.write_rerun_harness(datasets, self.figure_dpi)

    def _sandbox_builtins(self):
        """The builtins the generated code runs with."""
        return __builtins__

    def _specs_used_by(self, code: str) -> Dict[str, Any]:
        """
        Specs of every registered tool the code actually calls.
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
        sandbox_globals = {"__builtins__": self._sandbox_builtins()}
        
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
        
        sandbox_globals.update({"np": np, "pd": pd, "xr": xr, "plt": plt})
        plt.rcParams['font.sans-serif'] = ['Noto Sans CJK SC']
        plt.rcParams['axes.unicode_minus'] = False
        
        plt.close("all")

        print(f"[Info] Successfully registered {registered_count} tools to sandbox")

        if code_type == "compute":
            print(f"[Info] Executing computation code with {registered_count} registered tools")
        else:
            print(f"[Info] Executing visualization code")

        if context:
            sandbox_globals.update(context)
            
        predefined = set(sandbox_globals)
        try:
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

    def _materialize_outputs(self, output: Dict[str, Any]) -> Dict[str, Any]:
        """
        Compute the lazy arrays a successful compute block left behind, once, and
        put the computed arrays back in their place.
        """
        per_array = self.materialize_max_mb * 1024 ** 2
        budget = self.materialize_total_mb * 1024 ** 2
        found: Dict[int, Tuple[str, Any]] = {}

        def collect(obj: Any, path: str, depth: int) -> None:
            if depth > 4:
                return
            if isinstance(obj, (xr.DataArray, xr.Dataset)):
                if dask.is_dask_collection(obj) and id(obj) not in found:
                    found[id(obj)] = (path, obj)
            elif isinstance(obj, dict):
                for key, value in obj.items():
                    collect(value, f"{path}[{key!r}]", depth + 1)
            elif isinstance(obj, (list, tuple)):
                for i, value in enumerate(obj):
                    collect(value, f"{path}[{i}]", depth + 1)

        for name, value in output.items():
            if not name.startswith("_") and not callable(value) and not isinstance(value, types.ModuleType):
                collect(value, name, 0)
        if not found:
            return {"materialized": 0}

        chosen, kept_lazy, total = [], [], 0
        for path, obj in sorted(found.values(), key=lambda item: item[1].nbytes):
            if obj.nbytes <= per_array and total + obj.nbytes <= budget:
                chosen.append(obj)
                total += obj.nbytes
            else:
                kept_lazy.append(f"{path} ({obj.nbytes / 1024 ** 2:.0f} MB)")
        if not chosen:
            return {"materialized": 0, "kept_lazy": kept_lazy}

        started = datetime.now()
        try:
            computed = dask.compute(*chosen)
        except Exception as e:
            return {"materialized": 0, "kept_lazy": kept_lazy,
                    "failed": f"{type(e).__name__}: {e}"[:300]}
        replacement = {id(obj): new for obj, new in zip(chosen, computed)}

        def replace(obj: Any, depth: int) -> Any:
            if depth > 4:
                return obj
            if isinstance(obj, (xr.DataArray, xr.Dataset)):
                return replacement.get(id(obj), obj)
            if isinstance(obj, dict):
                for key in list(obj):
                    obj[key] = replace(obj[key], depth + 1)
                return obj
            if isinstance(obj, list):
                obj[:] = [replace(value, depth + 1) for value in obj]
                return obj
            if isinstance(obj, tuple):
                items = [replace(value, depth + 1) for value in obj]
                return type(obj)(*items) if hasattr(obj, "_fields") else tuple(items)
            return obj

        for name in list(output):
            output[name] = replace(output[name], 0)
        return {"materialized": len(chosen), "megabytes": round(total / 1024 ** 2, 1),
                "seconds": round((datetime.now() - started).total_seconds(), 1),
                "kept_lazy": kept_lazy}

    def _summarize_value(self, value: Any, max_chars: Optional[int] = None) -> Any:
        """
        Render one result variable for the solution step.
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
            "knowledge": record.get("knowledge"),
            "code_iterations": {"compute": steps(compute_iterations),
                                "visualization": steps(vis_iterations)},
            "figures": [os.path.basename(path) for path in record.get("visualization_paths", [])],
            "files": {
                "query": "query.txt",
                "knowledge_retrieval": "knowledge.json",
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
        best_partial: Optional[Dict[str, Any]] = None

        for iteration in range(self.max_code_iterations):
            print(f"\n--- Iteration {iteration + 1}/{self.max_code_iterations} ({code_type}) ---")
            
            exec_result = self._safe_execute_tool(current_code, context, code_type)

            # The results of a compute block are computed here, once, and those
            # computed arrays serve both the summary and the visualization block.
            materialized = None
            if exec_result["success"] and code_type == "compute":
                materialized = self._materialize_outputs(exec_result["output"])
                if materialized.get("materialized"):
                    print(f"[Info] Computed {materialized['materialized']} lazy result array(s) "
                          f"({materialized['megabytes']} MB) in {materialized['seconds']} s")
                if materialized.get("kept_lazy"):
                    print(f"[Info] Left lazy (over the size limit): {', '.join(materialized['kept_lazy'])}")
                if materialized.get("failed"):
                    print(f"[Warning] Computing the results together failed ({materialized['failed']}); "
                          f"they stay lazy and are computed one by one for the summary")

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
            if materialized:
                iteration_record["materialized_results"] = materialized

            if exec_result["success"]:
                if swallowed:
                    # Last iteration: keep the partial result rather than discarding
                    # it, but carry the failures into the observation so neither the
                    # answer nor a reader can miss them.
                    exec_result["swallowed_errors"] = swallowed

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

                fix_prompt = self.prompts.get_fix_prompt(
                    language, current_code, exec_result, user_query, specs_for_fix, code_type,
                    upstream_code=upstream_code, available_variables=available_variables
                )
                system_content = self.prompts.get_fix_system_message(language)
                
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
        session, bundle, result_record, started_at, context_summary = self._open_turn(
            user_query, session_id)
        language = result_record["language"]
        
        tool_overview = {
            k: (v.get("summary") or v["description"].strip().split("\n")[0])
            for k, v in self.tools_meta.items()
        }
        
        data_overview = {
            data["name"]: {"description": data.get("description", "")}
            for data in self.data_list
        }

        # Step 0: Knowledge retrieval. 
        knowledge_text, knowledge_record = self._retrieve_knowledge(
            user_query, context_summary, language, bundle)
        bundle.write_json("knowledge.json", knowledge_record)

        # Step 1: Planning, including selecting data and tools
        system_prompt_1 = self.prompts.get_planning_prompt(
            language, tool_overview, data_overview, context_summary, knowledge=knowledge_text
        )
        planning_input = self.prompts.get_planning_input(language, user_query)

        planning_response, planning_seconds = self._ask(
            bundle, "planning", system_prompt_1, planning_input
        )
        print("\n=== Plan & Selection ===\n", planning_response.content)
        result_record["plan"] = planning_response.content

        plan_match = (re.search(r"<plan>(.*?)</plan>", planning_response.content, re.S)
                      or re.search(r"<plan>(.*?)(?=<selected_data>|<selected_tools>|\Z)",
                                   planning_response.content, re.S))
        plan_text = plan_match.group(1).strip() if plan_match else ""
        bundle.write("plan.md", f"# Plan\n\n{plan_text or planning_response.content}\n")

        # Step 1b: Which of the fragments the planner was shown does the plan cite?
        # Checked against plan_text, since that alone travels on to the later steps.
        if knowledge_text:
            knowledge_record["plan_check"] = self.knowledge_base.check_plan(
                plan_text, knowledge_record["injected_ids"])
            uncited_must = knowledge_record["plan_check"]["uncited_must"]
            if uncited_must and self.knowledge_settings["on_uncited_must"] == "warn":
                print(f"[Warning] The plan does not cite these must-fragments it was shown: "
                      f"{', '.join(uncited_must)}")
            repeated = knowledge_record["plan_check"]["repeated_notes"]
            if repeated:
                print(f"[Info] The plan repeats notes for: "
                      f"{', '.join(f'{i} x{n}' for i, n in repeated.items())}")
            bundle.write_json("knowledge.json", knowledge_record)
        result_record["knowledge"] = {
            "selected": {s["id"]: s["status"] for s in knowledge_record.get("selected", [])},
            "selection_fallback": (knowledge_record.get("selection") or {}).get("fallback_reason"),
            "condensing_fallback": (knowledge_record.get("condensing_run") or {}).get("fallback_reason"),
            "injected_chars": knowledge_record.get("injected_chars"),
            "plan_check": knowledge_record.get("plan_check"),
        } if knowledge_record.get("enabled") else None

        # Extract the selected datasets.
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
            return self._close_turn_early(session, bundle, result_record, started_at)

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
        full_specs = self._code_reference(selected_tool_names, selected_datasets)

        system_prompt_2 = self.prompts.get_code_generation_prompt(
            language, full_specs, user_query, dataset_specs, self.dataset_notes, plan=plan_text
        )
        exec_input = self.prompts.get_code_generation_input(language)
        
        exec_response, _ = self._ask(bundle, "code_generation", system_prompt_2, exec_input)
        print("\n=== Generated Code ===\n", exec_response.content)
        result_record["generated_code"] = exec_response.content

        # The generated code addresses datasets by id. `selected_data` is a
        # convenience alias for the single-dataset case.
        code_context = {"selected_datasets": selected_datasets}
        if len(selected_datasets) == 1:
            code_context["selected_data"] = selected_datasets[0]

        executed = self._run_generated_code(exec_response.content, code_context, full_specs,
                                            user_query, language, bundle, result_record)
        if executed is None:
            return self._close_turn_early(session, bundle, result_record, started_at)
        observation_parts, compute_iterations, final_compute = executed
        return self._answer_and_finish(session, bundle, result_record, started_at,
                                       observation_parts, compute_iterations, final_compute,
                                       user_query, language, plan_text)

    def _open_turn(self, user_query: str, session_id: Optional[str]):
        """
        Start a turn: pick the language, find or open the session, and open the
        turn's bundle before the first model call.

        Returns the session, the bundle, the turn's result record (which carries
        the session id, the language and the turn number), the start time and the
        summary of the conversation so far.
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
        return session, bundle, result_record, started_at, context_summary

    def _close_turn_early(self, session: ConversationSession, bundle: RunBundle,
                          result_record: Dict[str, Any], started_at: datetime) -> Dict[str, Any]:
        """End a turn that stopped before any code ran; its record carries the reason."""
        bundle.write_json("manifest.json", self._manifest(result_record, started_at, [], []))
        session.add_turn(result_record)
        return result_record

    def _run_generated_code(self, generated: str, code_context: Dict[str, Any],
                            full_specs: Dict[str, Any], user_query: str, language: str,
                            bundle: RunBundle, result_record: Dict[str, Any]):
        """
        Execute the <compute> and <visualize> blocks of a code-generation answer,
        each with its repair loop, and gather what they produced.

        Returns the observation parts, the compute iterations and the compute
        code that finally ran; None when the answer has no <compute> block, with
        the reason in result_record.
        """
        turn_id = f"{result_record['session_id']}_turn{result_record['turn_number']}"
        # Step 3: Extract and execute the computation code
        compute_match = re.search(r"<compute>(.*?)</compute>", generated, re.S)
        if not compute_match:
            print("[Error] No <compute> code block found.")
            result_record["error"] = "No computation code found"
            bundle.write("code/raw_response.txt", generated)
            return None
        
        compute_code = compute_match.group(1).strip()

        self._write_rerun_harness(bundle, code_context.get("selected_datasets", []))

        compute_result, compute_iterations = self._execute_with_retry(
            compute_code,
            context=code_context,
            code_type="compute",
            full_specs=full_specs,
            user_query=user_query,
            language=language,
            bundle=bundle,
            available_variables={name: self._summarize_value(value)
                                 for name, value in code_context.items()}
        )
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
            visualize_match = re.search(r"<visualize>(.*?)</visualize>", generated, re.S)
            if visualize_match:
                visualize_code = visualize_match.group(1).strip()
                
                if visualize_code:
                    vis_context = dict(code_context)
                    vis_context.update(compute_result["output"])
                    
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
        return observation_parts, compute_iterations, final_compute

    def _answer_and_finish(self, session: ConversationSession, bundle: RunBundle,
                           result_record: Dict[str, Any], started_at: datetime,
                           observation_parts: List[Dict[str, Any]],
                           compute_iterations: List[Dict[str, Any]], final_compute: str,
                           user_query: str, language: str, plan_text: str = "") -> Dict[str, Any]:
        """Write the observation, ask for the final answer, and close the turn."""
        observation_text = json.dumps(observation_parts, indent=2, ensure_ascii=False)
        result_record["observation"] = observation_text
        bundle.write_json("observation.json", observation_parts)

        # Step 5: Generate final answer.
        has_computation_success = observation_parts[0]["status"] in (
            "success", "success_with_caught_failures")
        has_visualization = any(p.get("type") == "visualization" and p.get("status") == "success" 
                               for p in observation_parts)

        executed_code = final_compute
        if result_record.get("visualization_iterations"):
            executed_code += "\n\n# --- visualization ---\n" + \
                result_record["visualization_iterations"][-1]["code"]

        solution_prompt = self.prompts.get_solution_prompt(
            language, observation_text, has_computation_success, has_visualization,
            user_query=user_query, plan=plan_text, executed_code=executed_code
        )
        final_input = self.prompts.get_solution_input(language)

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
            # case, which model and when - T3_deepseek-reasoner_20260922_193101.
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


if __name__ == "__main__":
    # Built from configs/datasets.yaml.
    data_list = list_datasets()

    config_path = config_path_from_argv(CONFIG_PATH)
    config = load_config(config_path)
    arm = check_arm(config, ("full", "no_knowledge"), config_path)
    agent = AgentWithTools(TOOLS, data_list, config_path=config_path)

    batch = config.get("batch_run") or {}
    folder = resolve_path(batch.get("test_folder", "cases"))
    repeats = int(batch.get("repeats", 1))
    print(f"[System] Arm {arm}: running {folder} x{repeats} with {agent.model_name} "
          f"({agent.provider}; knowledge {'on' if agent.knowledge_enabled else 'off'}; "
          f"results -> {agent.result_dir})")

    for _ in range(repeats):
        results = run_tests_from_folder(str(folder), agent)

    print("Finished all tests.")
    
