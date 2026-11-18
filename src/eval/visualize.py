"""
Draw the evaluation figures from the score files in the results folders.

    python -m src.eval.visualize
    python -m src.eval.visualize --only a_portrait_llms e_loss_sources

Every results folder that configs/visualize.yaml names holds what the scorers write into it:

  - evaluation.json (src/eval/evaluate.py): the answer rate (RC) and accuracy (ACC) of
    every run, the difficulty level of every task, and which run belongs to which repeat.
  - visualization_scores.json and report_scores.json: the 0-3 completeness and quality
    of each run's figures and report, divided here by 3.
  - loss_attribution.json: the cause category of every numerical item the runs lost.

A score is computed as in the summary of evaluation.json: the mean over tasks of the
per-task means over runs, with the sample standard deviation across repeats of the
per-repeat means as its SD (`sd_over: tasks` takes it across the per-task means instead,
for folders with a single repeat). A difficulty level does the same over its own tasks.

Each figure is saved as <output_dir>/<name>.<format> and again as <name>_notext.<format>,
in which every text element is hidden and everything else, the layout included, is the
same.
"""

import argparse
import json
import re
import statistics
import sys
import textwrap
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import yaml
from matplotlib import colors as mcolors
from matplotlib import patheffects
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from matplotlib.text import Text
from matplotlib.transforms import blended_transform_factory, offset_copy

REPO_ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = REPO_ROOT / "configs" / "visualize.yaml"

# The batch runner names each session folder <case_id>_<model>_<timestamp>.
RUN_STAMP = re.compile(r"_(\d{8}_\d{6})$")

DEFAULT_FILES = {
    "evaluation": "evaluation.json",
    "visualization": "visualization_scores.json",
    "report": "report_scores.json",
    "loss": "loss_attribution.json",
}

METRICS = {
    "rc": {"source": "evaluation", "field": "answer_rate", "full": 1.0,
           "label": "RC", "group": "Numerical", "name": "RC"},
    "acc": {"source": "evaluation", "field": "accuracy", "full": 1.0,
            "label": "ACC", "group": "Numerical", "name": "ACC"},
    "figure_completeness": {"source": "visualization", "field": "completeness", "full": 3.0,
                            "label": "Compl.", "group": "Figures",
                            "name": "Figure completeness"},
    "figure_quality": {"source": "visualization", "field": "quality", "full": 3.0,
                       "label": "Qual.", "group": "Figures", "name": "Figure quality"},
    "report_completeness": {"source": "report", "field": "completeness", "full": 3.0,
                            "label": "Compl.", "group": "Report",
                            "name": "Report completeness"},
    "report_quality": {"source": "report", "field": "quality", "full": 3.0,
                       "label": "Qual.", "group": "Report", "name": "Report quality"},
}

# Portrait colours: by default cividis on a continuous, evenly spaced 0-1 scale.
# Configured bounds give discrete steps instead, and the steps from
# PORTRAIT_LIGHTEN_FROM up can be lightened by raising their OKLab lightness (hue and
# chroma kept) by the configured amounts.
PORTRAIT_CMAP = "cividis"
PORTRAIT_BOUNDS = None
PORTRAIT_LIGHTEN_FROM = 0.7
PORTRAIT_LIGHTEN: List[float] = []
# Plot-area sizes in inches, shared by the figures of a kind so that they line up.
LINE_AXES = [2.48, 2.02]
PORTRAIT_CELL = [0.48, 0.36]

SERIES_COLORS = ["#0072B2", "#E69F00", "#CC79A7", "#56B4E9", "#000000"]
MARKERS = ["o", "s", "^", "D", "v", "P"]

INK = "#1A1A1A"          # text
MUTED = "#6B6B6B"        # secondary text, separators
GRID = "#E3E3E3"
SURFACE = "#FFFFFF"
MISSING = "#EDEDED"      # a portrait cell without a score

# Loss categories of loss_attribution.json, without a red-green pair: every pair of the
# five clears the same separation as SERIES_COLORS. Category 4, acceptable numerical
# difference, is the neutral grey and textured.
LOSS_STYLE = {
    "1": {"color": "#E69F00", "hatch": None},
    "2": {"color": "#56B4E9", "hatch": None},
    "3": {"color": "#0072B2", "hatch": None},
    "4": {"color": "#DDDDDD", "hatch": "////"},
    "5": {"color": "#CC79A7", "hatch": None},
    "unattributed": {"color": SURFACE, "hatch": "xxxx"},
}
HATCH_COLOR = "#8C8C8C"
# Used when a folder has no loss_attribution.json to read the names from.
LOSS_NAMES = {
    "1": "No result produced",
    "2": "Method / convention error",
    "3": "Data handling error",
    "4": "Acceptable numerical difference",
    "5": "Other",
    "unattributed": "Not attributed",
}

LLM_STAGES = [
    {"key": "knowledge", "label": "Knowledge retrieval",
     "match": r"^knowledge_", "color": "#99DDFF", "hatch": None},
    {"key": "planning", "label": "Planning", "match": r"^planning", "color": "#EE8866",
     "hatch": None},
    {"key": "code_generation", "label": "Code generation", "match": r"^code_generation",
     "color": "#77AADD", "hatch": None},
    {"key": "code_repair", "label": "Code repair", "match": r"^fix_", "color": "#DDDDDD",
     "hatch": "////"},
    {"key": "solution", "label": "Answer & report", "match": r"^solution", "color": "#44BB99",
     "hatch": None},
]
OTHER_STAGE = {"key": "other", "label": "Other", "color": SURFACE, "hatch": "xxxx"}
# RunBundle.record_exchange heads every saved model call with its duration.
EXCHANGE_SECONDS = re.compile(r"^# elapsed: ([\d.]+) s$", re.M)
TURN_FOLDER = re.compile(r"^turn(\d+)$")

_warned = set()


def warn(message: str, prefix: str = "Warning") -> None:
    if message not in _warned:
        _warned.add(message)
        print(f"{prefix}: {message}", file=sys.stderr)


def resolve_path(value: Any) -> Path:
    path = Path(str(value)).expanduser()
    return path if path.is_absolute() else (REPO_ROOT / path).resolve()


def mean(values) -> Optional[float]:
    values = [v for v in values if v is not None]
    return statistics.fmean(values) if values else None


def stdev(values) -> Optional[float]:
    """Sample standard deviation (n - 1); None below two values."""
    values = [v for v in values if v is not None]
    return statistics.stdev(values) if len(values) > 1 else None


# ---------------------------------------------------------------------------
# Reading a results folder
# ---------------------------------------------------------------------------
def repeat_numbers(evaluation: Dict[str, Any]) -> Dict[str, int]:
    """Run name -> repeat number, from evaluation.json or, in older files, as evaluate.py
    forms them: the k-th run of every task by timestamp, complete repeats only."""
    if evaluation.get("repeats") is not None:
        return {run: rep["repeat"] for rep in evaluation["repeats"]
                for run in rep["runs"].values()}
    by_task = {task: sorted((run["run"] for run in runs),
                            key=lambda name: (RUN_STAMP.search(name) or [None, name])[1])
               for task, runs in evaluation["runs"].items()}
    n_repeats = min((len(names) for names in by_task.values()), default=0)
    return {names[k]: k + 1 for names in by_task.values() for k in range(n_repeats)}


def category_items(info: Dict[str, Any]) -> Dict[str, int]:
    """Lost items by category in one task of loss_attribution.json, which gives them as
    by_category {category: {items, events}} or, in single-run files, by_category_items."""
    if info.get("by_category") is not None:
        return {str(c): v.get("items", 0) for c, v in info["by_category"].items() if v.get("items")}
    return {str(c): n for c, n in (info.get("by_category_items") or {}).items() if n}


def scored_entries(per_task_value: Any) -> List[Dict[str, Any]]:
    """The scored runs of one task in a visualization or report score file: a single
    mapping when one run was scored, a list of them (or {"runs": [...]}) for repeats."""
    if isinstance(per_task_value, list):
        return per_task_value
    if isinstance(per_task_value, dict) and isinstance(per_task_value.get("runs"), list):
        return per_task_value["runs"]
    return [per_task_value] if isinstance(per_task_value, dict) else []


class Folder:
    """The score files of one results folder."""

    def __init__(self, entry: Dict[str, Any], fallback_levels: Dict[str, List[str]]):
        self.label = entry["label"]
        self.dir = resolve_path(entry["dir"])
        files = {**DEFAULT_FILES, **(entry.get("files") or {})}

        path = self.dir / files["evaluation"]
        if not path.is_file():
            raise SystemExit(f"{self.label}: {path} not found")
        evaluation = json.loads(path.read_text())
        self.levels = ((evaluation.get("settings") or {}).get("difficulty_levels")
                       or fallback_levels)
        self.repeat_of = repeat_numbers(evaluation)
        self.n_repeats = len(set(self.repeat_of.values()))
        self.runs = {run["run"] for runs in evaluation["runs"].values() for run in runs}
        # The task of every counted run, and the LLM time the evaluation recorded for it.
        self.task_of = {run["run"]: task for task, runs in evaluation["runs"].items()
                        for run in runs}
        self.recorded_llm_seconds = {run["run"]: (run.get("efficiency") or {}).get("llm_seconds")
                                     for runs in evaluation["runs"].values() for run in runs}
        self._llm_time_read = False
        # n_lost[task][run] = the graded items the run got wrong
        self.n_lost = {task: {run["run"]: run["scores"]["n_items"] - run["scores"]["n_correct"]
                              for run in runs
                              if run.get("scores") and run["scores"].get("n_items") is not None}
                       for task, runs in evaluation["runs"].items()}

        # scores[source][field][task][run] = the run's score
        self.scores: Dict[str, Dict[str, Dict[str, Dict[str, float]]]] = {}
        for field in ("accuracy", "answer_rate"):
            self.scores.setdefault("evaluation", {})[field] = {
                task: {run["run"]: run["scores"][field] for run in runs
                       if run.get("scores") and run["scores"].get(field) is not None}
                for task, runs in evaluation["runs"].items()}
        for source in ("visualization", "report"):
            data = self._read_optional(self.dir / files[source],
                                       f"{source} scores are left out")
            other_runs = set()
            for field in ("completeness", "quality"):
                by_task = {}
                for task, value in ((data or {}).get("per_task") or {}).items():
                    runs = {}
                    for i, scored in enumerate(scored_entries(value)):
                        if scored.get("run") and scored["run"] not in self.runs:
                            other_runs.add(scored["run"])
                        elif scored.get(field) is not None:
                            runs[scored.get("run") or f"{task}#{i}"] = scored[field]
                    by_task[task] = runs
                self.scores.setdefault(source, {})[field] = by_task
            if other_runs:
                warn(f"{self.label}: {files[source]} also scores {len(other_runs)} runs that "
                     f"{path.name} does not; left out", prefix="Note")
        self._loss_path = self.dir / files["loss"]
        self._loss: Optional[Dict[str, Any]] = None
        self._attributed: Optional[Dict[str, Dict[Optional[str], Dict[str, int]]]] = None

        summary_acc = (((evaluation.get("summary") or {}).get("scores") or {})
                       .get("accuracy") or {}).get("mean")
        own_acc = self.stat("evaluation", "accuracy")[0]
        if summary_acc is not None and own_acc is not None and abs(summary_acc - own_acc) > 1e-9:
            warn(f"{self.label}: accuracy {own_acc:.4f} recomputed from the runs differs from "
                 f"the summary of {path.name} ({summary_acc:.4f})")

    def _read_optional(self, path: Path, consequence: str) -> Optional[Dict[str, Any]]:
        if not path.is_file():
            warn(f"{self.label}: no {path.name} in {self.dir}; its {consequence}")
            return None
        return json.loads(path.read_text())

    @property
    def loss(self) -> Optional[Dict[str, Any]]:
        """loss_attribution.json, read when a figure first needs it."""
        if self._attributed is None:
            self._loss = self._read_optional(self._loss_path,
                                             "lost items are drawn as unattributed")
            self._attributed = self._read_attributed(self._loss or {})
        return self._loss

    @property
    def attributed(self) -> Dict[str, Dict[Optional[str], Dict[str, int]]]:
        self.loss
        return self._attributed

    def _read_attributed(self, loss: Dict[str, Any]) -> Dict[str, Dict[Optional[str], Dict[str, int]]]:
        """
        Lost items by category, as attributed[task][run][category].
        """
        items = loss.get("items") or []
        if items and all(item.get("run") for item in items):
            attributed: Dict[str, Dict[Optional[str], Dict[str, int]]] = {}
            for item in items:
                if item["run"] in self.runs:
                    counts = attributed.setdefault(item["task"], {}).setdefault(item["run"], {})
                    category = str(item["category"])
                    counts[category] = counts.get(category, 0) + 1
            return attributed
        return {task: {None: category_items(info or {})}
                for task, info in ((loss.get("summary") or {}).get("by_task") or {}).items()}

    def stat(self, source: str, field: str, tasks: Optional[List[str]] = None,
             sd_over: str = "repeats") -> Tuple[Optional[float], Optional[float]]:
        """Mean over tasks of the per-task means, and its SD (see the module docstring)."""
        by_task = self.scores.get(source, {}).get(field, {})
        wanted = list(tasks) if tasks is not None else sorted(
            set(by_task) | set(self.scores["evaluation"]["accuracy"]))
        missing = [task for task in wanted if not by_task.get(task)]
        if missing and by_task:
            warn(f"{self.label}: no {source} {field} score for {', '.join(missing)}; "
                 f"left out of its means")
        wanted = [task for task in wanted if by_task.get(task)]
        if not wanted:
            return None, None
        task_means = [mean(by_task[task].values()) for task in wanted]
        if sd_over == "tasks":
            return mean(task_means), stdev(task_means)
        per_repeat: Dict[int, Dict[str, float]] = {}
        for task in wanted:
            for run, value in by_task[task].items():
                if run in self.repeat_of:
                    per_repeat.setdefault(self.repeat_of[run], {})[task] = value
        # Only a repeat with a score for every task enters the SD.
        repeat_means = [mean(values.values()) for values in per_repeat.values()
                        if len(values) == len(wanted)]
        return mean(task_means), stdev(repeat_means)

    def loss_split(self, categories: List[str]) -> Dict[str, float]:
        """
        1 - ACC split by cause: each run's lost fraction 1 - accuracy, divided in the
        proportions of its lost items in loss_attribution.json, averaged over the task's
        runs and then over tasks, as ACC is. A run whose loss the file does not
        attribute, or a folder without the file, counts as unattributed, so the parts
        always add up to 1 - ACC.
        """
        accuracy = {task: runs for task, runs in self.scores["evaluation"]["accuracy"].items()
                    if runs}
        split = {key: 0.0 for key in [*categories, "unattributed"]}
        unattributed_runs, mismatched = 0, []
        for task, runs in accuracy.items():
            attributed = self.attributed.get(task, {})
            task_lost = sum(self.n_lost.get(task, {}).values())
            if None in attributed and sum(attributed[None].values()) not in (0, task_lost):
                mismatched.append(task)
            for run, acc in runs.items():
                lost = 1.0 - acc
                counts = attributed.get(run, attributed.get(None, {}))
                total = sum(counts.values())
                weight = 1.0 / len(runs) / len(accuracy)
                if run in attributed and total != self.n_lost.get(task, {}).get(run, total):
                    mismatched.append(run)
                if not total:
                    split["unattributed"] += lost * weight
                    unattributed_runs += lost > 0
                    continue
                for category, n in counts.items():
                    if category not in split:
                        warn(f"{self.label}: loss category {category} is not among "
                             f"{categories}; counted as unattributed")
                        category = "unattributed"
                    split[category] += lost * n / total * weight
        if unattributed_runs and self.loss:
            warn(f"{self.label}: loss_attribution.json does not attribute the loss of "
                 f"{unattributed_runs} runs; drawn as unattributed")
        if mismatched:
            warn(f"{self.label}: loss_attribution.json counts other lost items than the "
                 f"evaluation for {', '.join(mismatched)}; drawn with the evaluation's "
                 f"height and the file's proportions")
        return split

    def loss_events(self, categories: List[str]) -> Dict[str, float]:
        """Loss events per run, by category. An event is one decision or defect of a run,
        however many items it cost; events that own no item (co-causes only) are left
        out, as in the file's summary."""
        loss = self.loss or {}
        split = {key: 0.0 for key in [*categories, "unattributed"]}
        accuracy = self.scores["evaluation"]["accuracy"]
        events = loss.get("events") or []
        if events and all(event.get("run") for event in events):
            # A run counts when the file attributes its loss, or when it lost nothing.
            n_runs = sum(1 for task, runs in accuracy.items() for run in runs
                         if run in self.attributed.get(task, {})
                         or not self.n_lost.get(task, {}).get(run))
            counted = [(str(event["category"]), 1) for event in events
                       if event["run"] in self.runs and event.get("n_items", 1)]
        else:
            n_runs = sum(len(runs) for runs in accuracy.values())
            counted = [(str(category), info.get("events") or 0) for category, info in
                       (((loss.get("summary") or {}).get("by_category")) or {}).items()]
        if not counted or not n_runs:
            warn(f"{self.label}: no loss events to count; its bar is left empty")
            return split
        for category, n in counted:
            split[category if category in split else "unattributed"] += n / n_runs
        return split

    def read_llm_time(self) -> None:
        """
        Seconds of language-model calls of every counted run, by pipeline stage.
        """
        if self._llm_time_read:
            return
        self._llm_time_read = True
        keys = [stage["key"] for stage in LLM_STAGES] + [OTHER_STAGE["key"], "total"]
        by_stage = {key: {} for key in keys}
        missing, mismatched, other = [], [], set()
        for run, task in self.task_of.items():
            session = self.dir / run
            turns = ([p for p in session.iterdir() if p.is_dir() and TURN_FOLDER.match(p.name)]
                     if session.is_dir() else [])
            prompts = (max(turns, key=lambda p: int(TURN_FOLDER.match(p.name).group(1)))
                       / "prompts") if turns else None
            if prompts is None or not prompts.is_dir():
                missing.append(run)
                continue
            seconds = {key: 0.0 for key in keys}
            for path in sorted(prompts.glob("*.txt")):
                found = EXCHANGE_SECONDS.search(path.read_text(errors="replace"))
                name = re.sub(r"^\d+_", "", path.stem)
                stage = next((s["key"] for s in LLM_STAGES if re.search(s["match"], name)), None)
                if stage is None:
                    stage = OTHER_STAGE["key"]
                    other.add(name)
                value = float(found.group(1)) if found else 0.0
                seconds[stage] += value
                seconds["total"] += value
            recorded = self.recorded_llm_seconds.get(run)
            if recorded is not None and abs(recorded - seconds["total"]) > 0.5:
                mismatched.append(run)
            for key, value in seconds.items():
                by_stage[key].setdefault(task, {})[run] = value
        self.scores["llm_time"] = by_stage
        if missing:
            warn(f"{self.label}: no prompts/ folder for {len(missing)} runs; left out of the "
                 f"LLM time")
        if other:
            warn(f"{self.label}: LLM calls of no known stage ({', '.join(sorted(other))}); "
                 f"drawn as {OTHER_STAGE['label']}")
        if mismatched:
            warn(f"{self.label}: LLM time read from prompts/ differs from llm_seconds in the "
                 f"evaluation for {', '.join(mismatched)}")

    def loss_names(self) -> Dict[str, str]:
        names = dict(LOSS_NAMES)
        for key, category in ((self.loss or {}).get("categories") or {}).items():
            if category.get("name"):
                names[str(key)] = re.sub(r"\s*\([^)]*\)$", "", category["name"])
        return names


def resolve_entries(config: Dict[str, Any], refs: List[Any]) -> List[Tuple[str, Dict[str, Any]]]:
    """A figure's list of entries, as (key, entry): a name from `entries`, or inline."""
    named = config.get("entries") or {}
    resolved = []
    for i, ref in enumerate(refs):
        if isinstance(ref, str):
            if ref not in named:
                raise SystemExit(f"Entry {ref!r} is not defined under `entries`")
            resolved.append((ref, named[ref]))
        else:
            resolved.append((ref.get("name") or f"inline:{ref['label']}", ref))
    return resolved


# ---------------------------------------------------------------------------
# Drawing
# ---------------------------------------------------------------------------
def apply_style(config: Dict[str, Any]) -> None:
    font = config.get("font") or {}
    family = font.get("family") or ["DejaVu Sans"]
    size = float(font.get("size", 8))
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ([family] if isinstance(family, str) else list(family)) + ["DejaVu Sans"],
        "font.size": size,
        "axes.titlesize": size + 1,
        "axes.labelsize": size,
        "xtick.labelsize": size,
        "ytick.labelsize": size,
        "legend.fontsize": size,
        "text.color": INK,
        "axes.labelcolor": INK,
        "xtick.color": INK,
        "ytick.color": INK,
        "axes.edgecolor": MUTED,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.facecolor": SURFACE,
        "figure.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "legend.frameon": False,
        "hatch.linewidth": 0.5,
        # Text stays text, so that a layout program can edit it.
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })


def sequential_cmap(name: str):
    if name == "batlow":
        try:
            from cmcrameri import cm as crameri
        except ImportError:
            raise SystemExit("The batlow colour map needs the cmcrameri package "
                             "(pip install cmcrameri), or set portrait_colors.cmap: cividis")
        return crameri.batlow.copy()
    return plt.get_cmap(name).copy()


# Linear sRGB -> LMS and cube-rooted LMS -> OKLab (Ottosson 2020).
_TO_LMS = np.array([[0.4122214708, 0.5363325363, 0.0514459929],
                    [0.2119034982, 0.6806995451, 0.1073969566],
                    [0.0883024619, 0.2817188376, 0.6299787005]])
_TO_OKLAB = np.array([[0.2104542553, 0.7936177850, -0.0040720468],
                      [1.9779984951, -2.4285922050, 0.4505937099],
                      [0.0259040371, 0.7827717662, -0.8086757660]])


def lighten(color: Any, lift: float) -> str:
    """The colour with its OKLab lightness raised by `lift`, hue and chroma kept as far
    as sRGB allows."""
    srgb = np.array(mcolors.to_rgb(color))
    linear = np.where(srgb <= 0.04045, srgb / 12.92, ((srgb + 0.055) / 1.055) ** 2.4)
    lab = _TO_OKLAB @ np.cbrt(_TO_LMS @ linear)
    lab[0] += lift
    linear = np.clip(np.linalg.solve(_TO_LMS, np.linalg.solve(_TO_OKLAB, lab) ** 3), 0, 1)
    srgb = np.where(linear <= 0.0031308, 12.92 * linear, 1.055 * linear ** (1 / 2.4) - 0.055)
    return mcolors.to_hex(np.clip(srgb, 0, 1))


def portrait_colors(config: Dict[str, Any]):
    """Colour map, norm and colour-bar ticks shared by every portrait."""
    colors = config.get("portrait_colors") or {}
    base = sequential_cmap(colors.get("cmap", PORTRAIT_CMAP))
    bounds = colors.get("bounds", PORTRAIT_BOUNDS)
    if not bounds:
        base.set_bad(MISSING)
        return base, mcolors.Normalize(0.0, 1.0), [0, 0.25, 0.5, 0.75, 1.0]
    n = len(bounds) - 1
    steps = [mcolors.to_hex(base(i / (n - 1))) for i in range(n)]
    lighter = [i for i in range(n) if bounds[i] >= colors.get("lighten_from", PORTRAIT_LIGHTEN_FROM)]
    lifts = colors.get("lighten", PORTRAIT_LIGHTEN) or []
    if lifts and len(lifts) != len(lighter):
        warn(f"portrait_colors.lighten has {len(lifts)} values for {len(lighter)} steps; "
             f"the extra steps are left as they are")
    for i, lift in zip(lighter, lifts):
        steps[i] = lighten(steps[i], lift)
    cmap = mcolors.ListedColormap(steps)
    cmap.set_bad(MISSING)
    # clip: a score of exactly 1 takes the top step.
    return cmap, mcolors.BoundaryNorm(bounds, n, clip=True), list(bounds)


def text_color_on(rgba) -> str:
    """Dark or light ink, whichever reads on the cell colour."""
    r, g, b = (c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4 for c in rgba[:3])
    return INK if 0.2126 * r + 0.7152 * g + 0.0722 * b > 0.25 else SURFACE


def title_of(spec: Dict[str, Any]) -> Optional[str]:
    return spec.get("title") or None


def new_figure(size: List[float]) -> Tuple[plt.Figure, plt.Axes]:
    # Generous to start with; fit_axes shrinks it to the plot area asked for.
    fig = plt.figure(figsize=(size[0] + 4.0, size[1] + 3.0), layout="constrained")
    return fig, fig.add_subplot()


def fit_axes(fig: plt.Figure, ax: plt.Axes, size: List[float]) -> None:
    """Resize the figure until the plot area is `size` inches (width, height). The layout
    engine keeps every margin, so figures given one size share their plot area."""
    for _ in range(6):
        fig.canvas.draw()
        box = ax.get_position()
        width, height = fig.get_size_inches()
        dw, dh = size[0] - box.width * width, size[1] - box.height * height
        if abs(dw) < 0.003 and abs(dh) < 0.003:
            return
        fig.set_size_inches(width + dw, height + dh)


def level_name(level: str, spec: Dict[str, Any]) -> str:
    return (spec.get("level_labels") or {}).get(level) or re.sub(r"^level\s*(\d+)$", r"L\1", level)


def draw_grid(spec: Dict[str, Any], rows: List[Folder], columns: List[Dict[str, Any]],
              config: Dict[str, Any]) -> plt.Figure:
    """
    Entries against scores, coloured by the normalized score. Each column is
    {label, group, full, stat}: stat(folder) gives (mean, SD) in the column's own units.
    """
    annotate = spec.get("annotate", "mean_sd")
    values = np.full((len(rows), len(columns)), np.nan)
    spreads = np.full_like(values, np.nan)
    for i, folder in enumerate(rows):
        lacks_sd = False
        for j, column in enumerate(columns):
            m, sd = column["stat"](folder)
            if m is not None:
                values[i, j] = m / column["full"]
                if sd is not None:
                    spreads[i, j] = sd / column["full"]
                lacks_sd |= sd is None
        if lacks_sd and annotate == "mean_sd":
            warn(f"{folder.label}: no SD for some cells ({folder.n_repeats} complete "
                 f"repeat(s)); they show the mean only")

    cell_w, cell_h = spec.get("cell_size") or PORTRAIT_CELL
    size = [len(columns) * cell_w, len(rows) * cell_h]
    fig, ax = new_figure(size)
    cmap, norm, ticks = portrait_colors(config)
    mesh = ax.pcolormesh(np.ma.masked_invalid(values), cmap=cmap, norm=norm,
                         edgecolors=SURFACE, linewidth=1.5)
    ax.set_xlim(0, len(columns))
    ax.set_ylim(len(rows), 0)
    ax.set_xticks(np.arange(len(columns)) + 0.5, [c["label"] for c in columns])
    ax.set_yticks(np.arange(len(rows)) + 0.5, [folder.label for folder in rows])
    ax.xaxis.tick_top()
    ax.tick_params(length=0, pad=3)
    for spine in ax.spines.values():
        spine.set_visible(False)

    size_pt = plt.rcParams["font.size"]
    small = size_pt - 1
    for i in range(len(rows)):
        for j in range(len(columns)):
            value = values[i, j]
            if np.isnan(value):
                # Texture, not text, so that the text-free version still shows the gap.
                ax.add_patch(Rectangle((j + 0.04, i + 0.06), 0.92, 0.88, fill=False,
                                       hatch="////", edgecolor="#BDBDBD", linewidth=0))
                continue
            ink = text_color_on(cmap(norm(value)))
            if annotate == "mean_sd" and not np.isnan(spreads[i, j]):
                # The y axis points down: "bottom" sits above the cell centre.
                ax.text(j + 0.5, i + 0.5, f"{value:.2f}", ha="center", va="bottom", color=ink)
                ax.text(j + 0.5, i + 0.5, f"±{spreads[i, j]:.2f}", ha="center", va="top",
                        color=ink, fontsize=small)
            elif annotate in ("mean", "mean_sd"):
                ax.text(j + 0.5, i + 0.5, f"{value:.2f}", ha="center", va="center", color=ink)

    # Each run of columns sharing a group gets its name above their labels, over a rule.
    runs: List[List[Any]] = []
    for j, column in enumerate(columns):
        if runs and runs[-1][0] == column.get("group"):
            runs[-1][2] = j
        else:
            runs.append([column.get("group"), j, j])
    above = blended_transform_factory(ax.transData, ax.transAxes)
    for group, start, end in runs:
        if not group:
            continue
        ax.add_artist(Line2D([start + 0.08, end + 0.92], [1, 1], color=MUTED, linewidth=0.6,
                             transform=offset_copy(above, fig=fig, y=size_pt + 6, units="points"),
                             clip_on=False))
        ax.text((start + end + 1) / 2, 1, group, ha="center", va="bottom",
                transform=offset_copy(above, fig=fig, y=size_pt + 8, units="points"))

    cbar = fig.colorbar(mesh, ax=ax, location="bottom", aspect=30, pad=0.03, spacing="uniform")
    cbar.set_ticks(ticks, labels=[f"{t:g}" for t in ticks])
    cbar.set_label("Normalized score")
    cbar.outline.set_visible(False)
    cbar.ax.tick_params(length=2, width=0.6, labelsize=small)
    title = title_of(spec)
    if title:
        # A given y turns off matplotlib's own placement, which clears the column labels
        # but not the group names above them.
        grouped = any(group for group, _, _ in runs)
        ax.set_title(title, loc="left", y=1.0,
                     pad=(2.2 * size_pt + 14) if grouped else (1.2 * size_pt + 8))
    fit_axes(fig, ax, size)
    return fig


def draw_portrait(spec: Dict[str, Any], rows: List[Folder], config: Dict[str, Any]) -> plt.Figure:
    sd_over = spec.get("sd_over", config.get("sd_over", "repeats"))
    columns = []
    for key in spec.get("columns") or list(METRICS):
        metric = METRICS[key]
        columns.append({"label": metric["label"], "group": metric.get("group"),
                        "full": metric["full"],
                        "stat": lambda folder, m=metric: folder.stat(m["source"], m["field"],
                                                                     sd_over=sd_over)})
    return draw_grid(spec, rows, columns, config)


def series_style(entry: Dict[str, Any], i: int) -> Tuple[str, str]:
    return (entry.get("color") or SERIES_COLORS[i % len(SERIES_COLORS)],
            entry.get("marker") or MARKERS[i % len(MARKERS)])


def level_label(level: str, tasks: List[str], spec: Dict[str, Any]) -> str:
    name = level_name(level, spec)
    return f"{name}\n{', '.join(tasks)}" if spec.get("show_tasks", True) and tasks else name


def draw_level_lines(spec: Dict[str, Any], series: List[Tuple[Dict[str, Any], Folder]],
                     config: Dict[str, Any]) -> plt.Figure:
    metric = METRICS[spec.get("metric", "acc")]
    sd_over = spec.get("sd_over", config.get("sd_over", "repeats"))
    band = spec.get("band", "none")
    levels = spec.get("levels") or list(series[0][1].levels)
    for _, folder in series[1:]:
        if folder.levels != series[0][1].levels:
            warn(f"{folder.label}: difficulty levels differ from {series[0][1].label}'s; "
                 f"each line uses its own folder's")
    x = np.arange(len(levels), dtype=float)
    dodge = float(spec.get("dodge", 0.05))
    offsets = (np.arange(len(series)) - (len(series) - 1) / 2) * dodge

    size = spec.get("axes_size") or LINE_AXES
    fig, ax = new_figure(size)
    ax.yaxis.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    for i, (entry, folder) in enumerate(series):
        color, marker = series_style(entry, i)
        stats = [folder.stat(metric["source"], metric["field"], folder.levels.get(level, []),
                             sd_over) for level in levels]
        m = np.array([np.nan if s[0] is None else s[0] / metric["full"] for s in stats])
        sd = np.array([np.nan if s[1] is None else s[1] / metric["full"] for s in stats])
        xs = x + offsets[i]
        if band == "sd":
            if np.isnan(sd[~np.isnan(m)]).all():
                warn(f"{folder.label}: no SD by level ({folder.n_repeats} complete repeat(s), "
                     f"sd_over {sd_over}); drawn without a band")
            else:
                # A score outside 0-1 is impossible, so the band stops there.
                ax.fill_between(xs, np.clip(m - sd, 0, 1), np.clip(m + sd, 0, 1),
                                color=color, alpha=0.18, linewidth=0, zorder=1)
        ax.plot(xs, m, color=color, linewidth=1.5, marker=marker, markersize=5,
                markeredgecolor=SURFACE, markeredgewidth=0.8, label=folder.label,
                clip_on=False, zorder=3)

    first_levels = series[0][1].levels
    ax.set_xticks(x, [level_label(level, first_levels.get(level, []), spec) for level in levels])
    ax.set_xlim(-0.35, len(levels) - 0.65)
    ax.set_ylim(0, 1)
    ax.set_yticks(np.linspace(0, 1, 6))
    ax.set_xlabel(spec.get("xlabel", "Difficulty level"))
    ax.set_ylabel(spec.get("ylabel", metric["name"]))
    ax.tick_params(axis="x", length=0)
    fig.legend(loc="outside right upper", handlelength=2.2)
    title = title_of(spec)
    if title:
        ax.set_title(title, loc="left")
    fit_axes(fig, ax, size)
    return fig


def draw_loss_bars(spec: Dict[str, Any], groups: List[Tuple[Optional[str], List[Folder]]],
                   config: Dict[str, Any]) -> plt.Figure:
    categories = [str(c) for c in spec.get("categories") or ["1", "2", "3", "4", "5"]]
    events = spec.get("measure", "items") == "events"
    folders = [folder for _, members in groups for folder in members]
    names = next((f.loss_names() for f in folders if f.loss), dict(LOSS_NAMES))
    width, gap = 0.62, 0.6

    positions, spans, x = [], [], 0.0
    for label, members in groups:
        start = x
        for _ in members:
            positions.append(x)
            x += 1.0
        spans.append((label, start, x - 1.0))
        x += gap

    size = spec.get("axes_size") or [0.6 * (positions[-1] + 1) + 0.2, LINE_AXES[1]]
    fig, ax = new_figure(size)
    ax.yaxis.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    # Parts higher than this get their value written on them, or beside the bar when the
    # value does not fit in the part.
    label_parts = spec.get("label_parts")
    used, totals, part_labels = set(), [], []
    for position, folder in zip(positions, folders):
        split = folder.loss_events(categories) if events else folder.loss_split(categories)
        bottom = 0.0
        for key in [*categories, "unattributed"]:
            height = split.get(key, 0.0)
            if height <= 0:
                continue
            used.add(key)
            style = LOSS_STYLE.get(key, LOSS_STYLE["5"])
            # Fill and texture first; the surface-coloured outline drawn over it leaves
            # the gap between stacked parts (matplotlib draws hatches in the edge colour).
            ax.bar(position, height, width, bottom=bottom, color=style["color"],
                   hatch=style["hatch"], edgecolor=HATCH_COLOR, linewidth=0, zorder=2)
            ax.bar(position, height, width, bottom=bottom, fill=False,
                   edgecolor=SURFACE, linewidth=0.8, zorder=3)
            if label_parts is not None and height > float(label_parts):
                part_labels.append((position, bottom, height, style))
            bottom += height
        totals.append((position, bottom))

    top = max([total for _, total in totals] + [0.05])
    ax.set_ylim(0, top * (1.28 if any(label for label, _, _ in spans) else 1.15))
    label_size = plt.rcParams["font.size"] - 1
    # The height of a line of label text, in data units.
    line_height = (label_size / 72) * ax.get_ylim()[1] / size[1]
    side_floor = {}
    for position, bottom, height, style in part_labels:
        centre = bottom + height / 2
        if height < 1.15 * line_height:
            # Beside the bar, stacked in the order of the parts where they would collide.
            centre = max(centre, side_floor.get(position, centre))
            side_floor[position] = centre + 1.05 * line_height
            ax.text(position + width / 2 + 0.03, centre, f"{height:.2f}", ha="left",
                    va="center", color=INK, zorder=4, fontsize=label_size)
            continue
        # A patch of the fill behind the value keeps it clear of the hatching.
        fill = mcolors.to_rgba(style["color"])
        ax.text(position, centre, f"{height:.2f}", ha="center", va="center",
                color=text_color_on(fill), zorder=4, fontsize=label_size,
                bbox=dict(boxstyle="square,pad=0.15", facecolor=fill,
                          edgecolor="none") if style["hatch"] else None)
    if spec.get("show_total", True):
        for position, total in totals:
            ax.text(position, total + top * 0.025, f"{total:.1f}" if events else f"{total:.3f}",
                    ha="center", va="bottom", color=MUTED, fontsize=plt.rcParams["font.size"] - 1)
    ax.set_xticks(positions, [folder.label for folder in folders])
    rotation = float(spec.get("xtick_rotation", 0))
    if rotation:
        plt.setp(ax.get_xticklabels(), rotation=rotation, ha="right", rotation_mode="anchor")
    ax.tick_params(axis="x", length=0)
    ax.set_xlim(positions[0] - 0.6, positions[-1] + 0.6)
    ax.set_ylabel(spec.get("ylabel", "Loss events per run" if events
                           else "1 − ACC (fraction of items lost)"))

    above = blended_transform_factory(ax.transData, ax.transAxes)
    for i, (label, start, end) in enumerate(spans):
        if label:
            ax.text((start + end) / 2, 0.98, label, transform=above, ha="center", va="top",
                    color=MUTED)
        if i:
            ax.axvline(start - (1.0 + gap) / 2, color=GRID, linewidth=0.8, zorder=1)

    handles = [Patch(facecolor=LOSS_STYLE.get(key, LOSS_STYLE["5"])["color"],
                     hatch=LOSS_STYLE.get(key, LOSS_STYLE["5"])["hatch"],
                     edgecolor=HATCH_COLOR, linewidth=0.4 if key == "unattributed" else 0,
                     label=f"{key}  {names.get(key, key)}" if key != "unattributed"
                     else names["unattributed"])
               for key in [*categories, "unattributed"]
               if key in categories or key in used]
    # Top of the stack first, as the parts are seen.
    fig.legend(handles=handles[::-1], loc="outside right upper", handlelength=1.4,
               handleheight=1.1, title=spec.get("legend_title", "Cause of loss"),
               title_fontsize=plt.rcParams["font.size"], alignment="left")
    title = title_of(spec)
    if title:
        ax.set_title(title, loc="left")
    fit_axes(fig, ax, size)
    return fig


def draw_llm_time_bars(spec: Dict[str, Any], folders: List[Folder], config: Dict[str, Any]) -> plt.Figure:
    """
    One bar per entry: the mean language-model time of a run, as every score is averaged
    (over a task's runs, then over tasks), split by the pipeline stage of the calls. A
    whisker gives the SD of the total (sd_over, as for scores), and the total is written
    above it.
    """
    sd_over = spec.get("sd_over", config.get("sd_over", "repeats"))
    width = float(spec.get("bar_width", 0.62))
    positions = list(range(len(folders)))
    size = spec.get("axes_size") or [0.6 * len(folders) + 0.2, LINE_AXES[1]]
    fig, ax = new_figure(size)
    ax.yaxis.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)

    stages = [*LLM_STAGES, OTHER_STAGE]
    used, tops = set(), []
    for position, folder in zip(positions, folders):
        folder.read_llm_time()
        bottom = 0.0
        for stage in stages:
            height = folder.stat("llm_time", stage["key"], sd_over=sd_over)[0] or 0.0
            if height <= 0:
                continue
            used.add(stage["key"])
            ax.bar(position, height, width, bottom=bottom, color=stage["color"],
                   hatch=stage["hatch"], edgecolor=HATCH_COLOR, linewidth=0, zorder=2)
            ax.bar(position, height, width, bottom=bottom, fill=False,
                   edgecolor=SURFACE, linewidth=0.8, zorder=3)
            bottom += height
        total, sd = folder.stat("llm_time", "total", sd_over=sd_over)
        if total is None:
            warn(f"{folder.label}: no LLM time to draw; its bar is left empty")
            tops.append((position, None, None))
            continue
        if abs(total - bottom) > 0.05:
            warn(f"{folder.label}: LLM time stages add up to {bottom:.1f} s, not {total:.1f} s")
        if sd is not None:
            # A surface-coloured halo keeps the whisker legible where it crosses the parts
            # (the line only: on the caps it would draw them double).
            bars = ax.errorbar(position, total, yerr=sd, fmt="none", ecolor=INK, elinewidth=0.8,
                               capsize=2.5, capthick=0.8, zorder=4)
            halo = [patheffects.withStroke(linewidth=2.4, foreground=SURFACE)]
            for artist in bars.lines[2]:
                artist.set_path_effects(halo)
        elif sd_over != "tasks":
            warn(f"{folder.label}: no SD of the LLM time ({folder.n_repeats} complete "
                 f"repeat(s), sd_over {sd_over}); drawn without a whisker")
        tops.append((position, total, sd))

    highest = max([t + (s or 0) for _, t, s in tops if t is not None] + [1.0])
    ax.set_ylim(0, highest * 1.15)
    if spec.get("show_total", True):
        # On the top of the bar, beside the whisker, so that it reads as the bar's height.
        for position, total, sd in tops:
            if total is not None:
                ax.text(position + 0.04, total + highest * 0.012, f"{total:.0f}",
                        ha="left", va="bottom", color=INK,
                        fontsize=plt.rcParams["font.size"] - 1)
    wrap = spec.get("wrap_labels")
    labels = [textwrap.fill(folder.label, int(wrap)) if wrap else folder.label
              for folder in folders]
    ax.set_xticks(positions, labels)
    rotation = float(spec.get("xtick_rotation", 0))
    if rotation:
        plt.setp(ax.get_xticklabels(), rotation=rotation, ha="right", rotation_mode="anchor")
    ax.tick_params(axis="x", length=0)
    ax.set_xlim(positions[0] - 0.6, positions[-1] + 0.6)
    ax.set_ylabel(spec.get("ylabel", "LLM time per run (s)"))

    handles = [Patch(facecolor=stage["color"], hatch=stage["hatch"], edgecolor=HATCH_COLOR,
                     linewidth=0.4 if stage is OTHER_STAGE else 0, label=stage["label"])
               for stage in stages if stage is not OTHER_STAGE or stage["key"] in used]
    # Top of the stack first, as the parts are seen; then what the whisker is.
    handles = handles[::-1]
    sd_label = {"tasks": "± SD over tasks", "repeats": "± SD over repeats"}.get(sd_over)
    if sd_label:
        handles.append(Line2D([], [], color=INK, linewidth=0.8, marker="_", markersize=6,
                              markeredgewidth=0.8, label=f"Total, {sd_label}"))
    fig.legend(handles=handles, loc="outside right upper", handlelength=1.4,
               handleheight=1.1, title=spec.get("legend_title", "LLM call"),
               title_fontsize=plt.rcParams["font.size"], alignment="left")
    title = title_of(spec)
    if title:
        ax.set_title(title, loc="left")
    fit_axes(fig, ax, size)
    return fig


def save(fig: plt.Figure, out_dir: Path, name: str, formats: List[str], dpi: int) -> List[Path]:
    """Write the figure with all its text, then again with every text element hidden.
    The layout is fixed before the first write, so the two differ only in the text."""
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.canvas.draw()
    fig.set_layout_engine("none")
    written = []
    for fmt in formats:
        written.append(out_dir / f"{name}.{fmt}")
        fig.savefig(written[-1], dpi=dpi)
    for ax in fig.axes:
        ax.tick_params(which="both", labelbottom=False, labeltop=False,
                       labelleft=False, labelright=False)
    for text in fig.findobj(Text):
        text.set_visible(False)
    for fmt in formats:
        written.append(out_dir / f"{name}_notext.{fmt}")
        fig.savefig(written[-1], dpi=dpi)
    plt.close(fig)
    return written


def render(config: Dict[str, Any], only: Optional[List[str]], out_dir: Path) -> List[Path]:
    apply_style(config)
    fallback_levels = config.get("difficulty_levels") or {}
    formats = config.get("formats") or ["pdf", "png"]
    dpi = int(config.get("dpi", 600))
    figures = config.get("figures") or {}
    unknown = [name for name in only or [] if name not in figures]
    if unknown:
        raise SystemExit(f"No figure named {', '.join(unknown)} in the configuration")

    cache: Dict[str, Folder] = {}

    def folder_of(entry: Dict[str, Any]) -> Folder:
        cache_key = json.dumps([entry["label"], entry["dir"], entry.get("files")], sort_keys=True)
        if cache_key not in cache:
            cache[cache_key] = Folder(entry, fallback_levels)
        return cache[cache_key]

    written = []
    for name, spec in figures.items():
        if only and name not in only:
            continue
        kind = spec.get("kind")
        if kind == "portrait":
            rows = [folder_of(e) for _, e in resolve_entries(config, spec["rows"])]
            fig = draw_portrait(spec, rows, config)
        elif kind == "level_lines":
            series = [(e, folder_of(e)) for _, e in resolve_entries(config, spec["series"])]
            fig = draw_level_lines(spec, series, config)
        elif kind == "loss_bars":
            groups = spec.get("groups") or [{"label": None, "bars": spec["bars"]}]
            members = [(group.get("label"),
                        [folder_of(e) for _, e in resolve_entries(config, group["bars"])])
                       for group in groups]
            fig = draw_loss_bars(spec, members, config)
        elif kind == "llm_time_bars":
            bars = [folder_of(e) for _, e in resolve_entries(config, spec["bars"])]
            fig = draw_llm_time_bars(spec, bars, config)
        else:
            raise SystemExit(f"Figure {name}: unknown kind {kind!r}")
        written += save(fig, out_dir, name, formats, dpi)
        print(f"{name}: {', '.join(p.name for p in written[-2 * len(formats):])}")
    return written


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Draw the evaluation figures.")
    parser.add_argument("--config", type=Path, default=CONFIG_PATH,
                        help=f"configuration file (default: {CONFIG_PATH.relative_to(REPO_ROOT)})")
    parser.add_argument("--only", nargs="+", metavar="FIGURE",
                        help="draw only these figures, by their name in the configuration")
    parser.add_argument("--output-dir", type=Path, default=None,
                        help="overrides output_dir in the configuration")
    args = parser.parse_args(argv)

    config = yaml.safe_load(args.config.resolve().read_text()) or {}
    out_dir = (args.output_dir.resolve() if args.output_dir
               else resolve_path(config.get("output_dir", "results/figures")))
    render(config, args.only, out_dir)
    print(f"\nWrote {out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
