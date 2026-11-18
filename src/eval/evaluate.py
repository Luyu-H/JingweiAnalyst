"""
Score a folder of agent runs against the reference answers.

    python -m src.eval.evaluate
    python -m src.eval.evaluate --config configs/evaluate.yaml --result-dir results/<run>

What is graded in each task is listed in cases/reference_ans/grading_spec.yaml: every
entry points at values in that task's answer_<task>.json and states the precision the
task asks for. The agent's answer (solution.md) is free text, so a language model reads
it and copies out the value it reports for each graded item. That model is never shown
the reference values. Everything after that is arithmetic done here:

  - An item is correct when the reported value, rounded to the precision the task asks
    for, equals the reference value, or when its relative error
    |reported - reference| / |reference| is below the tolerance. Rankings, yes/no
    answers, dates and month indices must match exactly. A value the answer does not
    report counts as wrong.
  - A run's accuracy is the fraction of its items that are correct, its answer rate the
    fraction it reports a value for at all, and its relative error the mean over its
    numeric items.
  - Each run's efficiency is read from its bundle: wall time and code attempts from
    manifest.json, the number, duration and size of the model calls from prompts/.
  - A task's numbers are the mean and standard deviation over its runs. The k-th run of
    every task forms repeat k, and the summary gives the mean over tasks together with
    its standard deviation across repeats. Each difficulty level configured in
    configs/evaluate.yaml gets the same two numbers for answer rate and accuracy, over
    its own tasks.

Every setting lives in configs/evaluate.yaml. The output file, written into the results
folder, runs from the summary down to every graded item: summary, by_difficulty, per_task,
repeats, runs; its `definitions` section says what each metric means.
"""

import argparse
import fnmatch
import hashlib
import json
import math
import re
import statistics
import sys
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from decimal import ROUND_HALF_UP, Decimal, InvalidOperation
from pathlib import Path
from typing import Any, Dict, List, Optional

import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    # `python src/eval/evaluate.py` puts src/eval/ on the path, not the root that
    # src.models needs.
    sys.path.insert(0, str(REPO_ROOT))

CONFIG_PATH = REPO_ROOT / "configs" / "evaluate.yaml"

# The batch runner names each session folder <case_id>_<model>_<timestamp>.
TASK_FOLDER = re.compile(r"^(T\d+)_")
RUN_STAMP = re.compile(r"_(\d{8}_\d{6})$")
TURN_FOLDER = re.compile(r"^turn(\d+)$")
WILDCARD = re.compile(r"^\{(\w+)\}$")
# RunBundle.record_exchange heads every saved model call with its duration.
EXCHANGE_SECONDS = re.compile(r"^# elapsed: ([\d.]+) s$", re.M)

SCORE_KEYS = ("accuracy", "answer_rate", "mean_relative_error", "median_relative_error")
# The scores reported for each difficulty level.
LEVEL_KEYS = ("answer_rate", "accuracy")
EFFICIENCY_KEYS = ("elapsed_seconds", "llm_seconds", "non_llm_seconds", "llm_calls",
                   "code_attempts", "llm_input_chars", "llm_output_chars")

# Written into the output file, so that it can be read without this source.
DEFINITIONS = {
    "accuracy": "fraction of the graded items that are correct: the reported value, rounded "
                "to the precision the task asks for, equals the reference, or its relative "
                "error is below the tolerance; non-numeric items must match exactly; an item "
                "the answer does not report is wrong",
    "answer_rate": "fraction of the graded items for which the answer reports a usable "
                   "value, right or wrong",
    "mean_relative_error": "mean of |reported - reference| / |reference| over the run's "
                           "numeric items that have one; a zero reference has none, and an "
                           "unreported value counts as settings.scoring.missing_relative_error",
    "median_relative_error": "median of the same relative errors",
    "elapsed_seconds": "wall-clock time of the agent's turn, from its manifest.json",
    "llm_seconds": "summed duration of the turn's language-model calls, from prompts/",
    "non_llm_seconds": "elapsed_seconds minus llm_seconds: code execution, data loading and "
                       "everything else",
    "llm_calls": "number of language-model calls in the turn",
    "code_attempts": "code blocks executed, compute and visualization stages together; one "
                     "per stage means no repair was needed",
    "llm_input_chars": "characters sent to the model over all calls; the agent does not "
                       "record token counts, so characters stand in for them",
    "llm_output_chars": "characters the model returned over all calls",
    "per_task": "every metric as mean and sample standard deviation (n - 1) over the task's "
                "runs, n being the runs counted; std is null for a single run. Scores count "
                "every run whose answer could be read (a run without an answer scores 0); "
                "efficiency counts the runs whose turn completed, i.e. wrote a manifest",
    "repeats": "the k-th run of every task, in chronological order, forms repeat k, as the "
               "batch runner replays the whole case folder; only repeats that every task has "
               "are formed. Each value is the mean over the repeat's tasks",
    "summary": "every metric as the mean over tasks of the per-task means, and its sample "
               "standard deviation across repeats (null with fewer than two); totals sum "
               "the efficiency of every completed run",
    "by_difficulty": "for each difficulty level of settings.difficulty_levels, answer rate "
                     "and accuracy as the mean over the level's scored tasks of the per-task "
                     "means, and the sample standard deviation across repeats of the level's "
                     "per-repeat means; tasks of the level without a run are listed as missing",
}

# Numbers are graded by precision and relative error, everything else by exact match.
EXACT_TYPES = ("integer", "boolean", "list", "text")
ITEM_TYPES = ("number",) + EXACT_TYPES

JSON_TYPE_HINTS = {
    "number": "a JSON number",
    "integer": "a JSON integer",
    "boolean": "true or false",
    "list": "a JSON array, in the order the answer gives",
    "text": "a JSON string, in the format the unit names",
}


class SpecError(ValueError):
    """The grading spec does not match the reference answer it points into."""


def resolve_path(value: Any) -> Path:
    """A configured path, relative to the repository root unless absolute."""
    path = Path(str(value)).expanduser()
    return path if path.is_absolute() else REPO_ROOT / path


def sha256_of(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


# ---------------------------------------------------------------------------
# Graded items: the grading spec expanded over a reference answer
# ---------------------------------------------------------------------------

def is_number(value: Any) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def expand_entry(task: str, answer: Dict[str, Any], entry: Dict[str, Any]) -> List[Dict[str, Any]]:
    """
    One spec entry as graded items, one per value its path reaches.
    """
    segments = entry["path"].split(".")
    items: List[Dict[str, Any]] = []

    def walk(node: Any, depth: int, labels: List[str], bindings: Dict[str, str]) -> None:
        if depth == len(segments):
            items.append(make_item(task, entry, node, ".".join(labels), bindings))
            return
        segment = segments[depth]
        wildcard = WILDCARD.match(segment)
        if wildcard:
            if isinstance(node, dict):
                children = [(str(key), child) for key, child in node.items()]
            elif isinstance(node, list):
                children = [(str(i + 1), child) for i, child in enumerate(node)]
            else:
                raise SpecError(f"{task}: {entry['path']}: nothing to expand at "
                                f"{'.'.join(labels) or 'the top level'}")
            for label, child in children:
                walk(child, depth + 1, labels + [label],
                     {**bindings, wildcard.group(1): label})
        elif isinstance(node, dict) and segment in node:
            walk(node[segment], depth + 1, labels + [segment], bindings)
        elif isinstance(node, list) and segment.isdigit() and int(segment) < len(node):
            walk(node[int(segment)], depth + 1, labels + [segment], bindings)
        else:
            raise SpecError(f"{task}: {entry['path']}: answer_{task}.json has no "
                            f"'{segment}' under {'.'.join(labels) or 'the top level'}")

    walk(answer, 0, [], {})
    return items


def make_item(task: str, entry: Dict[str, Any], reference: Any, item_id: str,
              bindings: Dict[str, str]) -> Dict[str, Any]:
    kind = entry.get("type", "number")
    if kind not in ITEM_TYPES:
        raise SpecError(f"{task}: {entry['path']}: unknown type {kind!r}; one of {ITEM_TYPES}")
    if kind == "number":
        if not is_number(reference):
            raise SpecError(f"{task}: {item_id}: the reference value {reference!r} is not a number")
        if ("decimals" in entry) == ("sig_figs" in entry):
            raise SpecError(f"{task}: {entry['path']}: a number needs exactly one of "
                            f"decimals or sig_figs")
    try:
        # Keys such as "GFDL-ESM4|historical" read better in prose without the bar.
        description = entry["description"].format(
            **{name: label.replace("|", " ") for name, label in bindings.items()})
    except KeyError as exc:
        raise SpecError(f"{task}: {entry['path']}: the description names {exc}, "
                        f"which the path does not bind") from None
    return {
        "id": item_id,
        "description": description,
        "unit": str(entry.get("unit", "")),
        "type": kind,
        "decimals": entry.get("decimals"),
        "sig_figs": entry.get("sig_figs"),
        "transform": entry.get("transform"),
        "reference": reference,
    }


def load_items(spec: Dict[str, List[Dict[str, Any]]], reference_dir: Path,
               tasks: List[str]) -> Dict[str, List[Dict[str, Any]]]:
    """Every graded item of every task, checked against its reference answer."""
    items = {}
    for task in tasks:
        if task not in spec:
            raise SpecError(f"{task} is not in the grading spec")
        answer_path = reference_dir / task / f"answer_{task}.json"
        if not answer_path.exists():
            raise SpecError(f"{task}: no reference answer at {answer_path}")
        answer = json.loads(answer_path.read_text())
        task_items = [item for entry in spec[task] for item in expand_entry(task, answer, entry)]
        ids = [item["id"] for item in task_items]
        duplicates = sorted({i for i in ids if ids.count(i) > 1})
        if duplicates:
            raise SpecError(f"{task}: items listed twice: {', '.join(duplicates)}")
        items[task] = task_items
    return items


# ---------------------------------------------------------------------------
# Runs in the results folder
# ---------------------------------------------------------------------------

def difficulty_levels(config: Dict[str, Any], known_tasks: List[str]) -> Dict[str, List[str]]:
    """The configured levels, each a list of task ids; every task in at most one level."""
    levels: Dict[str, List[str]] = {}
    level_of: Dict[str, str] = {}
    for level, members in (config.get("difficulty_levels") or {}).items():
        levels[str(level)] = []
        for task in members or []:
            task = str(task)
            if task not in known_tasks:
                raise SystemExit(f"difficulty_levels: {level} names {task}, which is not a "
                                 f"task of the grading spec")
            if task in level_of:
                raise SystemExit(f"difficulty_levels: {task} is in both {level_of[task]} "
                                 f"and {level}")
            level_of[task] = str(level)
            levels[str(level)].append(task)
    return levels


def run_order(session: Path):
    """Chronological, by the timestamp ending the folder name."""
    stamp = RUN_STAMP.search(session.name)
    return (stamp.group(1) if stamp else "", session.name)


def find_runs(result_dir: Path) -> Dict[str, List[Path]]:
    """Session folders by task, oldest first."""
    runs: Dict[str, List[Path]] = {}
    for folder in sorted((p for p in result_dir.iterdir() if p.is_dir()), key=run_order):
        match = TASK_FOLDER.match(folder.name)
        if match:
            runs.setdefault(match.group(1), []).append(folder)
    return runs


def last_turn(session: Path) -> Optional[Path]:
    """
    The session's last turn, whose answer is the session's final one.

    A turn that crashed leaves no answer file; it is not replaced by an earlier turn's
    answer, which answered a different question.
    """
    turns = [p for p in session.iterdir() if p.is_dir() and TURN_FOLDER.match(p.name)]
    return max(turns, key=lambda p: int(TURN_FOLDER.match(p.name).group(1)), default=None)


def run_efficiency(turn: Optional[Path]) -> Dict[str, Any]:
    """
    What the turn cost, read from its bundle.
    """
    manifest_path = turn / "manifest.json" if turn else None
    manifest = (json.loads(manifest_path.read_text())
                if manifest_path and manifest_path.exists() else None)

    llm_calls, llm_seconds, chars_in, chars_out = 0, 0.0, 0, 0
    prompts = turn / "prompts" if turn else None
    for path in sorted(prompts.glob("*.txt")) if prompts and prompts.is_dir() else []:
        text = path.read_text(errors="replace")
        llm_calls += 1
        seconds = EXCHANGE_SECONDS.search(text)
        llm_seconds += float(seconds.group(1)) if seconds else 0.0
        sent, _, response = text.partition("\n===== RESPONSE =====\n")
        sent = sent.partition("===== SYSTEM =====\n")[2].replace("\n===== USER =====\n", "")
        chars_in += len(sent)
        chars_out += len(response)

    stages = ("compute", "visualization")
    if manifest:
        iterations = manifest.get("code_iterations") or {}
        attempts = {stage: len(iterations.get(stage) or []) for stage in stages}
        succeeded = {stage: (iterations[stage][-1].get("success") if iterations.get(stage)
                             else None) for stage in stages}
    else:
        tried = turn / "code" / "attempts" if turn else None
        attempts = {stage: (len(list(tried.glob(f"{stage}_*.py")))
                            if tried and tried.is_dir() else 0) for stage in stages}
        succeeded = {stage: None for stage in stages}

    elapsed = manifest.get("elapsed_seconds") if manifest else None
    return {
        "turn_completed": manifest is not None,
        "agent_status": manifest.get("status") if manifest else None,
        "elapsed_seconds": elapsed,
        "llm_seconds": round(llm_seconds, 1),
        "non_llm_seconds": round(elapsed - llm_seconds, 1) if elapsed is not None else None,
        "llm_calls": llm_calls,
        "code_attempts": sum(attempts.values()),
        "llm_input_chars": chars_in,
        "llm_output_chars": chars_out,
        "code_attempts_by_stage": attempts,
        "code_succeeded_by_stage": succeeded,
    }


def question_of(turn: Path) -> str:
    query = turn / "query.txt"
    if query.exists():
        return query.read_text()
    manifest = turn / "manifest.json"
    if manifest.exists():
        return json.loads(manifest.read_text()).get("user_query") or ""
    return ""


# ---------------------------------------------------------------------------
# Reading the reported values out of an answer
# ---------------------------------------------------------------------------

EXTRACT_SYSTEM = (
    "You read the final answer a data-analysis assistant gave to a question and copy out "
    "the values it reports. You never judge, correct, compute or guess a value."
)

EXTRACT_RULES = """\
Rules:
1. Copy each value as the answer states it, with its sign and every digit it gives. Do not
   round it, and do not compute a value the answer does not state, not even from other
   values it does state (no subtracting two reported means to get a difference it does
   not report).
2. If the answer gives a value in another unit or scale than the one listed (m3 instead of
   1e15 m3, a fraction instead of a percentage, a month name instead of its number),
   convert it to the listed unit.
3. If the answer gives more than one value for the same quantity, take the one it presents
   as its result.
4. Use model and source names exactly as they appear in the item descriptions.
5. If the answer does not report the value, give null. Never guess.
6. For every value that is not null, give as evidence the words of the answer it comes
   from, at most 20 words.

Reply with one JSON object and nothing else, of the form
{"values": {"<id>": {"value": <value or null>, "evidence": "<quote, or empty>"}}}
with every id listed above present exactly once."""


def extraction_prompt(question: str, answer: str, items: List[Dict[str, Any]]) -> str:
    lines = [f"- id: {item['id']} | {item['description']} | unit: {item['unit'] or 'n/a'} "
             f"| give: {JSON_TYPE_HINTS[item['type']]}" for item in items]
    return (f"## The question the assistant was given\n\n{question.strip()}\n\n"
            f"## The assistant's answer\n\n{answer.strip()}\n\n"
            f"## Values to copy out of the answer\n\n" + "\n".join(lines) +
            "\n\n" + EXTRACT_RULES)


def reply_text(reply: Any) -> str:
    content = getattr(reply, "content", reply)
    if isinstance(content, list):  # some clients return content blocks
        content = "".join(part.get("text", "") if isinstance(part, dict) else str(part)
                          for part in content)
    return str(content)


def parse_reply(text: str, ids: List[str]) -> Dict[str, Dict[str, Any]]:
    """The reply's value and evidence for each id; an id it leaves out is null."""
    fenced = re.search(r"```(?:json)?\s*(.*?)```", text, re.S)
    if fenced:
        text = fenced.group(1)
    start, end = text.find("{"), text.rfind("}")
    if start < 0 or end < start:
        raise ValueError("the reply holds no JSON object")
    data = json.loads(text[start:end + 1])
    values = data.get("values", data) if isinstance(data, dict) else None
    if not isinstance(values, dict) or not any(i in values for i in ids):
        raise ValueError("the reply names none of the requested ids")
    parsed = {}
    for item_id in ids:
        entry = values.get(item_id)
        if isinstance(entry, dict) and "value" in entry:
            parsed[item_id] = {"value": entry.get("value"),
                               "evidence": str(entry.get("evidence") or "")}
        else:
            parsed[item_id] = {"value": entry, "evidence": ""}
    return parsed


class Extractor:
    """The language model that reads answers, built only when an answer needs reading."""

    def __init__(self, cfg: Dict[str, Any], config_path: Path):
        from src.models.llm_providers import build_llm

        self.llm, self.settings = build_llm(cfg.get("llm") or {},
                                            hint=f"{config_path} (extractor.llm)")
        self.retries = int(cfg.get("retries", 2))

    @property
    def name(self) -> str:
        return f"{self.settings['provider']}:{self.settings['model']}"

    def __call__(self, question: str, answer: str,
                 items: List[Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
        from langchain_core.messages import HumanMessage, SystemMessage

        messages = [SystemMessage(content=EXTRACT_SYSTEM),
                    HumanMessage(content=extraction_prompt(question, answer, items))]
        ids = [item["id"] for item in items]
        error: Optional[Exception] = None
        for _ in range(self.retries + 1):
            try:
                return parse_reply(reply_text(self.llm.invoke(messages)), ids)
            except Exception as exc:  # a failed call and an unreadable reply alike
                error = exc
        raise RuntimeError(f"no readable reply after {self.retries + 1} attempts: {error}")


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------

def to_number(value: Any) -> Optional[float]:
    """A reported value as a float, or None when it is not a finite number."""
    if isinstance(value, bool) or value is None:
        return None
    if isinstance(value, (int, float)):
        return float(value) if math.isfinite(value) else None
    text = str(value).strip().replace("−", "-").replace(",", "").rstrip("%").strip()
    text = re.sub(r"\s*[x×]\s*10\s*\^\s*", "e", text)  # 2.2 x 10^17
    try:
        number = float(text)
    except ValueError:
        return None
    return number if math.isfinite(number) else None


def round_to_precision(value: float, decimals: Optional[int] = None,
                       sig_figs: Optional[int] = None) -> float:
    """Round half away from zero, as a person rounding to that precision would."""
    if value == 0:
        return 0.0
    if sig_figs is not None:
        decimals = sig_figs - 1 - math.floor(math.log10(abs(value)))
    try:
        step = Decimal(1).scaleb(-int(decimals))
        return float(Decimal(repr(value)).quantize(step, rounding=ROUND_HALF_UP))
    except InvalidOperation:
        return round(value, int(decimals))


def wrap_longitude(value: float) -> float:
    return (value + 180.0) % 360.0 - 180.0


def normalize_exact(value: Any, kind: str) -> Any:
    """A value in the form two exact-match items are compared in; None if unusable."""
    if value is None:
        return None
    if kind == "boolean":
        if isinstance(value, bool):
            return value
        return {"true": True, "yes": True, "false": False, "no": False}.get(str(value).strip().lower())
    if kind == "integer":
        number = to_number(value)
        return int(number) if number is not None and number == int(number) else None
    if kind == "list":
        if isinstance(value, str):
            value = [part for part in re.split(r"\s*(?:,|;|<|>|→|->)\s*", value) if part]
        if not isinstance(value, list):
            return None
        return [normalize_scalar(element) for element in value]
    return normalize_scalar(value)


def normalize_scalar(value: Any) -> Any:
    """Numbers compare as numbers; names and dates ignore case, spaces and punctuation."""
    number = to_number(value)
    if number is not None:
        return number
    return re.sub(r"[^0-9a-z]", "", str(value).lower())


def score_item(item: Dict[str, Any], extracted: Dict[str, Any], tolerance: float,
               scoring: Dict[str, Any]) -> Dict[str, Any]:
    reported = extracted.get("value")
    result = {
        "id": item["id"],
        "description": item["description"],
        "unit": item["unit"],
        "type": item["type"],
        "reference": item["reference"],
        "reported": reported,
        "evidence": extracted.get("evidence", ""),
    }

    if item["type"] != "number":
        expected = normalize_exact(item["reference"], item["type"])
        given = normalize_exact(reported, item["type"])
        result["status"] = ("missing" if reported is None
                            else "unreadable" if given is None else "reported")
        result["correct"] = given is not None and given == expected
        return result

    reference = float(item["reference"])
    value = to_number(reported)
    if value is None:
        missing = scoring.get("missing_relative_error")
        result.update(status="missing" if reported is None else "unreadable", correct=False,
                      relative_error=None if missing is None else float(missing),
                      tolerance=tolerance)
        return result

    if item["transform"] == "longitude":
        value, reference = wrap_longitude(value), wrap_longitude(reference)
    precision = {"decimals": item["decimals"], "sig_figs": item["sig_figs"]}
    rounded = round_to_precision(value, **precision)
    exact = math.isclose(rounded, round_to_precision(reference, **precision),
                         rel_tol=1e-9, abs_tol=1e-12)
    graded = rounded if scoring.get("round_to_required_precision", True) else value
    # A reference of zero has no relative error; such an item is right only when exact.
    relative_error = abs(graded - reference) / abs(reference) if reference != 0 else None
    result.update(
        status="reported",
        graded_value=graded,
        exact_match=exact,
        relative_error=relative_error,
        tolerance=tolerance,
        correct=exact or (relative_error is not None and relative_error < tolerance),
    )
    return result


def matches_any(task: str, item_id: str, patterns) -> bool:
    key = f"{task}/{item_id}"
    return any(fnmatch.fnmatchcase(key, str(pattern)) for pattern in patterns or [])


def tolerance_for(task: str, item_id: str, scoring: Dict[str, Any]) -> float:
    key = f"{task}/{item_id}"
    for pattern, tolerance in (scoring.get("tolerance_overrides") or {}).items():
        if fnmatch.fnmatchcase(key, str(pattern)):
            return float(tolerance)
    return float(scoring.get("relative_tolerance", 0.05))


def graded_items(task: str, items: List[Dict[str, Any]],
                 scoring: Dict[str, Any]) -> List[Dict[str, Any]]:
    """The task's items that the scoring settings keep."""
    keep_exact = scoring.get("include_non_numeric", True)
    return [item for item in items
            if (keep_exact or item["type"] == "number")
            and not matches_any(task, item["id"], scoring.get("exclude_items"))]


def mean(values: List[Optional[float]]) -> Optional[float]:
    values = [v for v in values if v is not None]
    return statistics.fmean(values) if values else None


def stdev(values: List[Optional[float]]) -> Optional[float]:
    """Sample standard deviation (n - 1); None below two values."""
    values = [v for v in values if v is not None]
    return statistics.stdev(values) if len(values) > 1 else None


def mean_std(values: List[Optional[float]]) -> Dict[str, Any]:
    return {"mean": mean(values), "std": stdev(values),
            "n": sum(v is not None for v in values)}


def run_scores(scored: List[Dict[str, Any]]) -> Dict[str, Any]:
    errors = [item["relative_error"] for item in scored
              if item["type"] == "number" and item["relative_error"] is not None]
    n_correct = sum(item["correct"] for item in scored)
    n_answered = sum(item["status"] == "reported" for item in scored)
    return {
        "accuracy": n_correct / len(scored) if scored else None,
        "answer_rate": n_answered / len(scored) if scored else None,
        "mean_relative_error": statistics.fmean(errors) if errors else None,
        "median_relative_error": statistics.median(errors) if errors else None,
        "n_items": len(scored),
        "n_correct": n_correct,
        "n_answered": n_answered,
        "n_numeric": sum(item["type"] == "number" for item in scored),
        "n_relative_error": len(errors),
    }


# ---------------------------------------------------------------------------
# The evaluation
# ---------------------------------------------------------------------------

def previous_extractions(output_path: Path) -> Dict[str, Dict[str, Any]]:
    """What an earlier evaluation of this folder read from each run, by session folder."""
    if not output_path.exists():
        return {}
    try:
        previous = json.loads(output_path.read_text())
    except (json.JSONDecodeError, OSError):
        return {}
    runs = [run for task_runs in (previous.get("runs") or {}).values() for run in task_runs]
    # Files written before runs moved to the top level kept them under tasks.<task>.runs.
    runs += [run for task in (previous.get("tasks") or {}).values() for run in task.get("runs", [])]
    return {run["run"]: run["extraction"] for run in runs if run.get("extraction")}


def read_answers(jobs: List[Dict[str, Any]], items: Dict[str, List[Dict[str, Any]]],
                 extractor_cfg: Dict[str, Any], config_path: Path, output_path: Path) -> None:
    """Attach to every run with an answer the values read from it, reusing earlier reads."""
    reuse = previous_extractions(output_path) if extractor_cfg.get("reuse_extraction", True) else {}
    model = extractor_model_name(extractor_cfg)
    pending = []
    for job in jobs:
        if job["answer"] is None:
            continue
        job["sha256"] = sha256_of(job["answer"])
        cached = reuse.get(job["session"].name)
        if (cached and cached.get("answer_sha256") == job["sha256"]
                and cached.get("model") == model
                and all(item["id"] in cached.get("values", {}) for item in items[job["task"]])):
            job["extraction"] = cached
        else:
            pending.append(job)
    if not pending:
        return

    extractor = Extractor(extractor_cfg, config_path)
    print(f"[evaluate] Reading {len(pending)} answer(s) with {extractor.name}"
          f" ({len(jobs) - len(pending)} reused or without an answer)")

    def read(job):
        try:
            values = extractor(question_of(job["turn"]), job["answer"].read_text(),
                               items[job["task"]])
            return job, {"model": extractor.name, "answer_sha256": job["sha256"],
                         "values": values}, None
        except Exception as exc:
            return job, None, str(exc)

    with ThreadPoolExecutor(max_workers=max(1, int(extractor_cfg.get("max_workers", 4)))) as pool:
        for job, extraction, error in pool.map(read, pending):
            print(f"[evaluate]   {job['session'].name}: "
                  f"{'read' if extraction else 'FAILED - ' + error}")
            job["extraction"], job["error"] = extraction, error


def score_run(job: Dict[str, Any], task: str, kept: List[Dict[str, Any]],
              scoring: Dict[str, Any], result_dir: Path) -> Dict[str, Any]:
    """One run: its scores, its efficiency, and every graded item."""
    extraction = job.get("extraction")
    if job["answer"] is None:
        status, values = "no_answer", {}       # crashed or still running: every item missing
    elif extraction is None:
        status, values = "extraction_failed", None
    else:
        status, values = "scored", extraction["values"]
    scored = ([score_item(item, values.get(item["id"], {}),
                          tolerance_for(task, item["id"], scoring), scoring) for item in kept]
              if values is not None else [])

    run: Dict[str, Any] = {"run": job["session"].name, "task": task, "status": status}
    if status == "extraction_failed":
        run["error"] = job.get("error")
    run["answer_file"] = str(job["answer"].relative_to(result_dir)) if job["answer"] else None
    run["scores"] = run_scores(scored) if values is not None else None
    run["efficiency"] = run_efficiency(job["turn"])
    run["items"] = scored
    if extraction:
        run["extraction"] = extraction
    return run


def task_summary(runs: List[Dict[str, Any]], n_items: int) -> Dict[str, Any]:
    """Mean and standard deviation of every metric over the task's runs."""
    scored = [run for run in runs if run["scores"]]
    completed = [run for run in runs if run["efficiency"]["turn_completed"]]
    return {
        "n_items": n_items,
        "n_runs": len(runs),
        "n_runs_scored": len(scored),
        "n_runs_completed": len(completed),
        "scores": {key: mean_std([run["scores"][key] for run in scored]) for key in SCORE_KEYS},
        "efficiency": {key: mean_std([run["efficiency"][key] for run in completed])
                       for key in EFFICIENCY_KEYS},
    }


def form_repeats(runs_by_task: Dict[str, List[Dict[str, Any]]],
                 levels: Dict[str, List[str]]):
    """
    Repeat k is the k-th run of every task, oldest first.
    """
    if not runs_by_task:
        return [], []
    n_repeats = min(len(runs) for runs in runs_by_task.values())
    repeats = []
    for k in range(n_repeats):
        members = [runs[k] for runs in runs_by_task.values()]
        repeats.append({
            "repeat": k + 1,
            "runs": {run["task"]: run["run"] for run in members},
            "scores": {key: mean([run["scores"][key] if run["scores"] else None
                                  for run in members]) for key in SCORE_KEYS},
            "efficiency": {key: mean([run["efficiency"][key]
                                      if run["efficiency"]["turn_completed"] else None
                                      for run in members]) for key in EFFICIENCY_KEYS},
            "by_difficulty": {level: {key: mean([run["scores"][key] if run["scores"] else None
                                                 for run in members if run["task"] in tasks])
                                      for key in LEVEL_KEYS}
                              for level, tasks in levels.items()},
        })
    unpaired = [run["run"] for runs in runs_by_task.values() for run in runs[n_repeats:]]
    return repeats, unpaired


def evaluate(config: Dict[str, Any], config_path: Path, result_dir: Path) -> Dict[str, Any]:
    scoring = config.get("scoring") or {}
    extractor_cfg = config.get("extractor") or {}
    answer_name = config.get("answer_file", "solution.md")
    output_path = result_dir / config.get("output_file", "evaluation.json")
    spec_path = resolve_path(config.get("grading_spec", "cases/reference_ans/grading_spec.yaml"))
    spec = yaml.safe_load(spec_path.read_text()) or {}
    reference_dir = resolve_path(config.get("reference_dir", "cases/reference_ans"))

    tasks = [str(t) for t in (config.get("tasks") or spec)]
    items = load_items(spec, reference_dir, tasks)
    levels = difficulty_levels(config, [str(t) for t in spec])
    level_of = {task: level for level, members in levels.items() for task in members}
    found = find_runs(result_dir)

    jobs = []
    for task in tasks:
        for session in found.get(task, []):
            turn = last_turn(session)
            answer = turn / answer_name if turn else None
            jobs.append({"task": task, "session": session, "turn": turn,
                         "answer": answer if answer and answer.exists() else None})
    read_answers(jobs, items, extractor_cfg, config_path, output_path)

    runs_by_task: Dict[str, List[Dict[str, Any]]] = {}
    per_task: Dict[str, Dict[str, Any]] = {}
    for task in tasks:
        task_jobs = [job for job in jobs if job["task"] == task]
        if not task_jobs:
            continue
        kept = graded_items(task, items[task], scoring)
        runs_by_task[task] = [score_run(job, task, kept, scoring, result_dir) for job in task_jobs]
        per_task[task] = {"difficulty": level_of.get(task),
                          **task_summary(runs_by_task[task], len(kept))}
    repeats, unpaired = form_repeats(runs_by_task, levels)

    def overall(section: str, keys) -> Dict[str, Any]:
        return {key: {"mean": mean([per_task[t][section][key]["mean"] for t in per_task]),
                      "std_across_repeats": stdev([rep[section][key] for rep in repeats])}
                for key in keys}

    by_difficulty = {}
    for level, members in levels.items():
        scored = [task for task in members if task in per_task]
        by_difficulty[level] = {
            "tasks": members,
            "missing_tasks": [task for task in members if task not in per_task],
            **{key: {"mean": mean([per_task[t]["scores"][key]["mean"] for t in scored]),
                     "std_across_repeats": stdev([rep["by_difficulty"][level][key]
                                                  for rep in repeats])}
               for key in LEVEL_KEYS},
        }

    all_runs = [run for runs in runs_by_task.values() for run in runs]
    completed = [run["efficiency"] for run in all_runs if run["efficiency"]["turn_completed"]]
    summary = {
        "n_tasks": len(per_task),
        "n_runs": len(all_runs),
        "n_repeats": len(repeats),
        "missing_tasks": [task for task in tasks if task not in per_task],
        "failed_runs": [run["run"] for run in all_runs if run["status"] == "extraction_failed"],
        "unpaired_runs": unpaired,
        "tasks_without_difficulty": [task for task in per_task if levels and task not in level_of],
        "scores": overall("scores", SCORE_KEYS),
        "efficiency": overall("efficiency", EFFICIENCY_KEYS),
        "totals": {"runs_completed": len(completed),
                   **{key: round(sum(e[key] for e in completed if e[key] is not None), 1)
                      for key in EFFICIENCY_KEYS}},
    }
    return {
        "evaluated_at": datetime.now().isoformat(timespec="seconds"),
        "result_dir": str(result_dir),
        "settings": {
            "config": str(config_path),
            "grading_spec": str(spec_path),
            "grading_spec_sha256": sha256_of(spec_path),
            "reference_dir": str(reference_dir),
            "answer_file": answer_name,
            "scoring": scoring,
            "extractor_model": extractor_model_name(extractor_cfg),
            "difficulty_levels": levels,
        },
        "definitions": DEFINITIONS,
        "summary": summary,
        "by_difficulty": by_difficulty,
        "per_task": per_task,
        "repeats": repeats,
        "runs": runs_by_task,
    }


def extractor_model_name(extractor_cfg: Dict[str, Any]) -> str:
    """The provider:model an extraction is recorded under, without building a client."""
    from src.models.llm_providers import resolve_settings

    settings = resolve_settings(extractor_cfg.get("llm") or {})
    return f"{settings['provider']}:{settings['model']}"


def print_summary(report: Dict[str, Any]) -> None:
    def fmt(value, spread=None, digits=3):
        if value is None:
            return "-"
        text = f"{value:.{digits}f}"
        return text + (f" ± {spread:.{digits}f}" if spread is not None else "")

    def pct(value):
        return "-" if value is None else f"{100 * value:.2f}%"

    def minutes(value):
        return "-" if value is None else f"{value / 60:.1f}"

    header = (f"{'task':<6} {'runs':>4} {'items':>5}  {'answer rate':<15} {'accuracy':<15} "
              f"{'mean rel. err':>13} {'min/run':>8} {'LLM calls':>9} {'code tries':>10}")
    print("\n" + header + "\n" + "-" * len(header))
    for task, row in report["per_task"].items():
        s, e = row["scores"], row["efficiency"]
        print(f"{task:<6} {row['n_runs']:>4} {row['n_items']:>5}  "
              f"{fmt(s['answer_rate']['mean'], s['answer_rate']['std']):<15} "
              f"{fmt(s['accuracy']['mean'], s['accuracy']['std']):<15} "
              f"{pct(s['mean_relative_error']['mean']):>13} "
              f"{minutes(e['elapsed_seconds']['mean']):>8} "
              f"{fmt(e['llm_calls']['mean'], digits=1):>9} "
              f"{fmt(e['code_attempts']['mean'], digits=1):>10}")
    summary = report["summary"]
    s, e = summary["scores"], summary["efficiency"]
    print("-" * len(header))
    print(f"{'mean':<6} {'':>4} {'':>5}  "
          f"{fmt(s['answer_rate']['mean'], s['answer_rate']['std_across_repeats']):<15} "
          f"{fmt(s['accuracy']['mean'], s['accuracy']['std_across_repeats']):<15} "
          f"{pct(s['mean_relative_error']['mean']):>13} "
          f"{minutes(e['elapsed_seconds']['mean']):>8} "
          f"{fmt(e['llm_calls']['mean'], digits=1):>9} "
          f"{fmt(e['code_attempts']['mean'], digits=1):>10}")
    if report["by_difficulty"]:
        print(f"\n{'level':<8} {'tasks':<14} {'answer rate':<15} {'accuracy':<15}")
        for level, row in report["by_difficulty"].items():
            missing = f"  (no run: {' '.join(row['missing_tasks'])})" if row["missing_tasks"] else ""
            print(f"{level:<8} {' '.join(row['tasks']):<14} "
                  f"{fmt(row['answer_rate']['mean'], row['answer_rate']['std_across_repeats']):<15} "
                  f"{fmt(row['accuracy']['mean'], row['accuracy']['std_across_repeats']):<15}"
                  f"{missing}")
        if summary["tasks_without_difficulty"]:
            print(f"Tasks in no difficulty level: {', '.join(summary['tasks_without_difficulty'])}")

    totals = summary["totals"]
    print(f"\n± is the std over a task's runs in task rows, and across the {summary['n_repeats']} "
          f"repeat(s) in the mean and level rows. {totals['runs_completed']} completed run(s) took "
          f"{totals['elapsed_seconds'] / 3600:.2f} h in all, "
          f"{totals['llm_seconds'] / 3600:.2f} h of it waiting for the model.")
    if summary["missing_tasks"]:
        print(f"No run found for: {', '.join(summary['missing_tasks'])}")
    if summary["unpaired_runs"]:
        print(f"Runs in no complete repeat: {', '.join(summary['unpaired_runs'])}")
    if summary["failed_runs"]:
        print(f"Answers that could not be read, left out of the scores: "
              f"{', '.join(summary['failed_runs'])}")


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Score a results folder against the reference answers.")
    parser.add_argument("--config", type=Path, default=CONFIG_PATH,
                        help=f"configuration file (default: {CONFIG_PATH.relative_to(REPO_ROOT)})")
    parser.add_argument("--result-dir", type=Path, default=None,
                        help="results folder to score; overrides result_dir in the configuration")
    args = parser.parse_args(argv)

    config_path = args.config.resolve()
    config = yaml.safe_load(config_path.read_text()) or {}
    result_dir = (args.result_dir.resolve() if args.result_dir
                  else resolve_path(config.get("result_dir", "")))
    if not result_dir.is_dir():
        raise SystemExit(f"Results folder not found: {result_dir}")

    try:
        report = evaluate(config, config_path, result_dir)
    except SpecError as exc:
        raise SystemExit(f"Grading spec problem: {exc}")

    output_path = result_dir / config.get("output_file", "evaluation.json")
    output_path.write_text(json.dumps(report, indent=2, ensure_ascii=False))
    print_summary(report)
    print(f"\nWrote {output_path}")
    return 1 if report["summary"]["failed_runs"] else 0


if __name__ == "__main__":
    sys.exit(main())
