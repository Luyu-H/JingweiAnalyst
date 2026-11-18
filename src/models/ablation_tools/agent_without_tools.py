"""
The tool ablation (-T): the full agent with none of the project's tools.

Run the case batch with

    python -m src.models.ablation_tools.agent_without_tools

configured by configs/agent_with_knowledge_no_tools.yaml, which extends
configs/agent_full.yaml and changes only where results are written.
"""
import builtins
from pathlib import Path
from typing import Any, Dict, List, Optional


from src.models.agent_with_tools import (AgentWithTools, REPO_ROOT, RERUN_HARNESS, RunBundle,
                                         check_arm, config_path_from_argv, load_config,
                                         resolve_path, run_tests_from_folder)
from src.models.ablation_tools.prompt import PromptTemplatesWithoutTools
from src.tools.data_openers import dataset_notes, list_datasets, load_dataset_config


CONFIG_PATH = REPO_ROOT / "configs" / "agent_with_knowledge_no_tools.yaml"


# The registry's cross_dataset_notes, restated for the raw files. The originals
# describe what the openers deliver ("the opener always outputs -180-180"); with
# no opener, what the code meets is the native state, so each note gives that
# instead. Same topics, same facts, taken from configs/datasets.yaml.
RAW_CROSS_DATASET_NOTES = {
    "longitude": "Native conventions differ: CMIP6 and ERSSTv5 are 0-360, WOA23 and SEANOE -180-180, GLODAP "
                 "20.5-379.5. Nothing converts them; bring products to one convention before comparing or "
                 "sub-setting across them.",
    "depth": "CESM2's lev is in centimetres (its lev_bnds declare metres), the rest in metres; all positive "
             "down. Nothing converts them.",
    "grid": "CESM2 and MPI-ESM1-2-LR are curvilinear and cannot be sel'd by lat/lon directly; use a 2-D mask. "
            "Area weighting of model fields uses the model's areacello.",
    "time": "WOA23 needs decode_times=False; SEANOE's time carries no CF units; CMIP6 uses three different "
            "calendars (CESM2 noleap / 365_day, GFDL-ESM4 noleap, MPI-ESM1-2-LR proleptic_gregorian), so "
            "their time axes decode to different types and must be aligned before they are combined.",
    "units": "CMIP6 uses mol m-3 (per volume), WOA23 and GLODAP use umol kg-1 (per mass); converting between "
             "them needs in situ density, mol m-3 = umol kg-1 * rho / 1000. Nothing converts them.",
    "climatology_period": "WOA23 temperature/salinity 1971-2000 or 1991-2020, oxygen/AOU/percent-saturation "
                          "1971-2000 or 1965-2022, nutrients (nitrate/phosphate/silicate) 1965-2022 ONLY; "
                          "GLODAP 1972-2013; ERSSTv5 internal climatology 1971-2000; CPC uses 1991-2020 from "
                          "2021-01 onward.",
}


RAW_RERUN_HARNESS = RERUN_HARNESS
for _old, _new in (
    ("The agent runs generated code with every registered tool already in scope and\n"
     "the selected dataset ids bound to `selected_datasets`,",
     "The tool-ablation agent (src/models/ablation_tools) runs generated code with no\n"
     "project tools in scope, only the scientific stack, and the selected dataset ids\n"
     "bound to `selected_datasets`,"),
    ("from src.tools.registry import TOOLS, auto_import_tools\nauto_import_tools()\n\n"
     "namespace = {name: meta[\"func\"] for name, meta in TOOLS.items()}\n",
     "namespace = {}\n"),
):
    if _old not in RAW_RERUN_HARNESS:
        raise ValueError("RERUN_HARNESS in src/models/agent_with_tools.py changed; update "
                         f"src/models/ablation_tools/agent_without_tools.py. Missing text:\n{_old}")
    RAW_RERUN_HARNESS = RAW_RERUN_HARNESS.replace(_old, _new, 1)


def _is_project_module(name: str) -> bool:
    return name == "src" or name.startswith("src.")


def _import_without_project(name, globals=None, locals=None, fromlist=(), level=0):
    """__import__ for the generated code: everything installed, nothing of this project."""
    if level == 0 and _is_project_module(name):
        raise ImportError(f"No module named {name!r} in this setting: the project's own modules are "
                          f"not available. Open the data files directly.")
    return builtins.__import__(name, globals, locals, fromlist, level)


class AgentWithoutTools(AgentWithTools):
    """AgentWithTools with no tools: the generated code reads and processes the raw files itself."""

    prompts = PromptTemplatesWithoutTools

    def __init__(self, data_list: List[Dict[str, str]], api_key: Optional[str] = None,
                 config: Optional[Dict[str, Any]] = None, config_path: Optional[Path] = None,
                 **overrides):
        """
        Args are AgentWithTools's without `tools`: the registry is always empty. The
        configuration defaults to configs/agent_with_knowledge_no_tools.yaml.
        """
        config_from_caller = config is not None
        config_path = Path(config_path) if config_path is not None else CONFIG_PATH
        if config is None:
            config = load_config(config_path)
        super().__init__({}, data_list, api_key=api_key, config=config,
                         config_path=config_path, **overrides)
        self.config_from_caller = config_from_caller

        # The raw facts replace what the openers would have normalised. Checked
        # against the registry's own topics, so a convention added there is not
        # silently missing here.
        registry_topics = set(dataset_notes())
        if registry_topics != set(RAW_CROSS_DATASET_NOTES):
            raise ValueError(
                "cross_dataset_notes in configs/datasets.yaml and RAW_CROSS_DATASET_NOTES in "
                f"{Path(__file__).relative_to(REPO_ROOT)} cover different topics: "
                f"{sorted(registry_topics ^ set(RAW_CROSS_DATASET_NOTES))}")
        self.dataset_notes = dict(RAW_CROSS_DATASET_NOTES)
        self.registry_entries = load_dataset_config()["datasets"]

        self.sandbox_builtins = dict(vars(builtins))
        self.sandbox_builtins["__import__"] = _import_without_project

    def _code_reference(self, selected_tools: List[str],
                        selected_datasets: List[str]) -> Dict[str, Any]:
        """The raw registry entry of each selected dataset, in place of tool specs."""
        return {name: self.registry_entries[name] for name in selected_datasets
                if name in self.registry_entries}

    def _write_rerun_harness(self, bundle: RunBundle, datasets: List[str]) -> None:
        bundle.write_rerun_harness(datasets, self.figure_dpi, template=RAW_RERUN_HARNESS)

    def _sandbox_builtins(self):
        return self.sandbox_builtins

    def _configuration_snapshot(self) -> Dict[str, Any]:
        snapshot = super()._configuration_snapshot()
        snapshot["always_available_tools"] = []
        snapshot["ablation"] = {
            "name": "without_tools",
            "agent": f"{type(self).__module__}.{type(self).__name__}",
            "prompts": f"{self.prompts.__module__}.{self.prompts.__name__}",
            "code_reference": "raw registry entries of the selected datasets (configs/datasets.yaml)",
            "cross_dataset_notes": self.dataset_notes,
            "project_imports_blocked": True,
            "knowledge": self.knowledge_enabled,
        }
        return snapshot


if __name__ == "__main__":
    data_list = list_datasets()

    config_path = config_path_from_argv(CONFIG_PATH)
    arm = check_arm(load_config(config_path), ("no_tools",), config_path)
    agent = AgentWithoutTools(data_list, config_path=config_path)

    batch = agent.config.get("batch_run") or {}
    folder = resolve_path(batch.get("test_folder", "cases"))
    repeats = int(batch.get("repeats", 1))
    print(f"[System] Arm {arm}: running {folder} x{repeats} with {agent.model_name} "
          f"({agent.provider}; no tools; knowledge {'on' if agent.knowledge_enabled else 'off'}; "
          f"results -> {agent.result_dir})")

    for _ in range(repeats):
        results = run_tests_from_folder(str(folder), agent)

    print("Finished all tests.")
