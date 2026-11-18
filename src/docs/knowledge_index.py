"""
The catalogue of knowledge fragments in docs/knowledge, and the checks that keep it sound.

Each fragment is one Markdown file with YAML front matter (the format is described in
docs/knowledge/README.md). Two things are built from those headers:

  INDEX.md    one line per fragment - id, severity, variables and summary - which is
              what a retrieval step shows the model when it chooses fragments, so that
              only the chosen ones are loaded in full.
  checks      the header fields exist, the id equals the file name, the type matches
              the folder, and every name a fragment points at actually exists. 

Run from anywhere:

    python -m src.docs.knowledge_index            # check, then rewrite INDEX.md
    python -m src.docs.knowledge_index --check    # check only; exit 1 on any problem
"""

import argparse
import re
import sys
from pathlib import Path
from typing import Any, Dict, List

import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
KNOWLEDGE_ROOT = REPO_ROOT / "docs" / "knowledge"
TOOLS_ROOT = REPO_ROOT / "src" / "tools"
DATASETS_CONFIG = REPO_ROOT / "configs" / "datasets.yaml"

# Folder -> the `type` every fragment inside it must declare.
FOLDER_TYPES = {
    "pitfalls": "pitfall",
    "concepts": "concept",
    "methods": "method",
    "datasets": "dataset",
}
# Index sections, in the order a planner should read them: silent errors first.
SECTION_TITLES = {
    "pitfall": "Pitfalls: silent errors",
    "concept": "Concepts: what a variable or attribute means",
    "method": "Methods: how to compute or compare",
    "dataset": "Datasets: models, experiments, grids",
}
SEVERITY_ORDER = {"must": 0, "should": 1, "info": 2}

REQUIRED_FIELDS = ("id", "title", "type", "summary", "keywords", "severity", "sources")
SOURCE_FIELDS = ("citation", "locator", "url")
# Sections the retrieval step may cut out and inject on their own.
REQUIRED_SECTIONS = ("## Key points", "## Implications for this project", "## Source excerpts")

_FRONT_MATTER = re.compile(r"^---\n(.*?)\n---\n(.*)$", re.S)
_BODY_REFERENCE = re.compile(r"`([a-z0-9]+(?:-[a-z0-9]+){2,})`")
_REGISTERED_TOOL = re.compile(r'@register_tools\(\s*name="([A-Za-z_0-9]+)"')


def load_fragments(root: Path = KNOWLEDGE_ROOT) -> List[Dict[str, Any]]:
    """
    Every fragment under `root`, as {front-matter fields..., "path", "body"}.
    """
    fragments = []
    for folder in FOLDER_TYPES:
        for path in sorted((root / folder).glob("*.md")):
            match = _FRONT_MATTER.match(path.read_text(encoding="utf-8"))
            if not match:
                raise ValueError(f"{path}: no YAML front matter between '---' lines")
            try:
                header = yaml.safe_load(match.group(1)) or {}
            except yaml.YAMLError as exc:
                # The usual cause is an unquoted title or summary containing ': '.
                raise ValueError(f"{path}: front matter is not valid YAML ({exc})") from exc
            header["path"] = path
            header["body"] = match.group(2)
            fragments.append(header)
    return fragments


def registered_tools() -> set:
    names = set()
    for path in TOOLS_ROOT.glob("*.py"):
        names.update(_REGISTERED_TOOL.findall(path.read_text(encoding="utf-8")))
    return names


def registry_datasets() -> set:
    return set((yaml.safe_load(DATASETS_CONFIG.read_text(encoding="utf-8")) or {}).get("datasets", {}))


def check_fragments(fragments: List[Dict[str, Any]], root: Path = KNOWLEDGE_ROOT) -> List[str]:
    """Every problem found, one line each; an empty list means the catalogue is sound."""
    problems = []
    tools, datasets = registered_tools(), registry_datasets()

    seen = {}
    for fragment in fragments:
        fid = fragment.get("id")
        if fid in seen:
            problems.append(f"{fid}: duplicate id in {fragment['path']} and {seen[fid]}")
        seen[fid] = fragment["path"]
    ids = set(seen)

    for fragment in fragments:
        path, fid = fragment["path"], fragment.get("id", "?")
        where = path.relative_to(root)

        for key in REQUIRED_FIELDS:
            if not fragment.get(key):
                problems.append(f"{where}: missing '{key}'")
        if fid != path.stem:
            problems.append(f"{where}: id '{fid}' does not match the file name")
        expected_type = FOLDER_TYPES[path.parent.name]
        if fragment.get("type") != expected_type:
            problems.append(f"{where}: type '{fragment.get('type')}', but the folder holds '{expected_type}'")
        if fragment.get("severity") not in SEVERITY_ORDER:
            problems.append(f"{where}: severity '{fragment.get('severity')}' is not one of {list(SEVERITY_ORDER)}")

        for source in fragment.get("sources") or []:
            missing = [key for key in SOURCE_FIELDS if not source.get(key)]
            if missing:
                problems.append(f"{where}: a source lacks {missing}")

        for related in fragment.get("related") or []:
            if related not in ids:
                problems.append(f"{where}: related '{related}' is not a fragment")
        for tool in fragment.get("tools") or []:
            if tool not in tools:
                problems.append(f"{where}: tool '{tool}' is not registered in src/tools")
        for dataset in fragment.get("datasets") or []:
            if dataset not in datasets:
                problems.append(f"{where}: dataset '{dataset}' is not in configs/datasets.yaml")

        body = fragment["body"]
        for reference in _BODY_REFERENCE.findall(body):
            if reference not in ids:
                problems.append(f"{where}: body refers to `{reference}`, which is not a fragment")
        for section in REQUIRED_SECTIONS:
            if section not in body:
                problems.append(f"{where}: missing section '{section}'")
    return problems


def render_index(fragments: List[Dict[str, Any]]) -> str:
    lines = [
        "# Index of knowledge fragments",
        "",
        "One line per fragment. The `summary` column is the text a retrieval step can show to "
        "the model when it picks fragments; the full file is loaded only for the fragments it "
        "picks.",
        "",
        "Generated by `python -m src.docs.knowledge_index` from the front matter of each "
        "fragment. Do not edit by hand; edit the fragment and regenerate.",
        "",
        "Severity: **must** = ignoring it produces a wrong number; **should** = ignoring it "
        "produces a misleading comparison or interpretation; **info** = background for "
        "interpretation.",
        "",
    ]
    for ftype, title in SECTION_TITLES.items():
        members = sorted(
            (f for f in fragments if f.get("type") == ftype),
            key=lambda f: (SEVERITY_ORDER.get(f.get("severity"), 9), f["id"]),
        )
        if not members:
            continue
        lines += [f"## {title}", "", "| id | severity | variables | summary |", "|---|---|---|---|"]
        for f in members:
            link = f["path"].relative_to(KNOWLEDGE_ROOT).as_posix()
            variables = ", ".join(f.get("variables") or []) or "–"
            summary = " ".join(str(f["summary"]).split()).replace("|", "/")
            lines.append(f"| [{f['id']}]({link}) | {f['severity']} | {variables} | {summary} |")
        lines.append("")
    return "\n".join(lines)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    parser.add_argument("--check", action="store_true",
                        help="only check the fragments; do not rewrite INDEX.md")
    args = parser.parse_args(argv)

    try:
        fragments = load_fragments()
    except ValueError as exc:
        print(f"[error] {exc}")
        return 1

    problems = check_fragments(fragments)
    for problem in problems:
        print(f"[problem] {problem}")
    print(f"{len(fragments)} fragments, {len(problems)} problem(s)")
    if problems:
        # An index built from a broken catalogue would advertise fragments that
        # point at nothing, so it is not rewritten until the problems are fixed.
        return 1
    if not args.check:
        index_path = KNOWLEDGE_ROOT / "INDEX.md"
        index_path.write_text(render_index(fragments), encoding="utf-8")
        print(f"wrote {index_path.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
