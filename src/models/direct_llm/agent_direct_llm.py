"""
The direct-LLM baseline: no tools, no knowledge.

Run the case batch with

    python -m src.models.direct_llm.agent_direct_llm

configured by configs/direct_llm.yaml.
"""
from pathlib import Path
from typing import Any, Dict, List, Optional

from src.models.agent_with_tools import (REPO_ROOT, check_arm, config_path_from_argv, load_config,
                                         resolve_path, run_tests_from_folder)
from src.models.ablation_tools.agent_without_tools import AgentWithoutTools
from src.models.direct_llm.prompt import DirectLLMPrompts
from src.tools.data_openers import list_datasets


CONFIG_PATH = REPO_ROOT / "configs" / "direct_llm.yaml"


class DirectLLMAgent(AgentWithoutTools):
    """Question and data description in, code out, run, repaired on failure, answered."""

    prompts = DirectLLMPrompts

    def __init__(self, data_list: List[Dict[str, str]], api_key: Optional[str] = None,
                 config: Optional[Dict[str, Any]] = None, config_path: Optional[Path] = None,
                 **overrides):
        """
        Args are AgentWithoutTools's. The configuration defaults to
        configs/direct_llm.yaml; its knowledge section is ignored, since this
        baseline is defined without knowledge.
        """
        config_from_caller = config is not None
        config_path = Path(config_path) if config_path is not None else CONFIG_PATH
        config = dict(config) if config is not None else load_config(config_path)
        if (config.get("knowledge") or {}).get("enabled"):
            print("[System] Direct-LLM baseline: the knowledge section is ignored (no retrieval)")
        config["knowledge"] = {"enabled": False}
        super().__init__(data_list, api_key=api_key, config=config, config_path=config_path,
                         **overrides)
        self.config_from_caller = config_from_caller

    def _configuration_snapshot(self) -> Dict[str, Any]:
        snapshot = super()._configuration_snapshot()
        snapshot["ablation"].update({
            "name": "direct_llm",
            "agent": f"{type(self).__module__}.{type(self).__name__}",
            "prompts": f"{self.prompts.__module__}.{self.prompts.__name__}",
            "code_reference": "raw registry entries of every dataset (configs/datasets.yaml)",
            "planning": False,
            "knowledge": False,
        })
        return snapshot

    def _datasets_in(self, code: str) -> List[str]:
        """The registered datasets whose root directory the code names."""
        return [name for name, data in self.data_index.items() if data["path"] in code]

    def run(self, user_query: str, session_id: Optional[str] = None) -> Dict[str, Any]:
        """
        Answer a user query within a conversation session.

        Returns the result record, shaped like the full agent's, without plan or
        tool selection.
        """
        session, bundle, result_record, started_at, context_summary = self._open_turn(
            user_query, session_id)
        language = result_record["language"]
        # Written like the other arms' bundles, so a reader finds the same files.
        bundle.write_json("knowledge.json", {"enabled": False})
        result_record["selected_tools"] = []
        result_record["knowledge"] = None

        # Step 1: code generation, straight from the question and the data.
        registry = dict(self.registry_entries)
        system_prompt = self.prompts.get_code_generation_prompt(
            language, registry, self.dataset_notes, context_summary)
        response, _ = self._ask(bundle, "code_generation", system_prompt,
                                self.prompts.get_code_generation_input(language, user_query))
        print("\n=== Generated Code ===\n", response.content)
        result_record["generated_code"] = response.content

        # Steps 2-3: execution with repair, then the answer - the full agent's own.
        # Nothing is pre-bound in the namespace: no dataset was selected beforehand.
        executed = self._run_generated_code(response.content, {}, registry, user_query,
                                            language, bundle, result_record)
        if executed is None:
            return self._close_turn_early(session, bundle, result_record, started_at)
        observation_parts, compute_iterations, final_compute = executed

        # Which datasets the code actually read, found from the roots it names, so
        # the manifest can say what the answer rests on as it does for the other arms.
        code_ran = final_compute + "".join(
            item["code"] for item in result_record.get("visualization_iterations", [])[-1:])
        used = self._datasets_in(code_ran)
        result_record["selected_data"] = {"name": used,
                                          "path": [self.data_index[name]["path"] for name in used],
                                          "how": "inferred from the root directories the executed code names"}

        return self._answer_and_finish(session, bundle, result_record, started_at,
                                       observation_parts, compute_iterations, final_compute,
                                       user_query, language)


if __name__ == "__main__":
    config_path = config_path_from_argv(CONFIG_PATH)
    check_arm(load_config(config_path), ("direct_llm",), config_path)
    agent = DirectLLMAgent(list_datasets(), config_path=config_path)

    batch = agent.config.get("batch_run") or {}
    folder = resolve_path(batch.get("test_folder", "cases"))
    repeats = int(batch.get("repeats", 1))
    print(f"[System] Direct-LLM baseline: running {folder} x{repeats} with {agent.model_name} "
          f"({agent.provider}; results -> {agent.result_dir})")

    for _ in range(repeats):
        results = run_tests_from_folder(str(folder), agent)

    print("Finished all tests.")
