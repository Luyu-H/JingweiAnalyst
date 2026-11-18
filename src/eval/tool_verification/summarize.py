"""
Aggregate a verification run into report.md and summary.json.

    python -m src.eval.tool_verification.summarize results/tool_verification/<run>
"""
import glob
import json
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..")))

CATEGORY_TITLES = {
    "cat1_analytic": "1 analytic",
    "cat2_reference": "2 reference impl.",
    "cat3_products": "3 published products",
    "cat4_magnitude": "4 magnitudes",
    "cat5_consistency": "5 consistency",
    "cat6_montecarlo": "6 Monte Carlo",
    "cat7_safeguards": "7 safeguards",
}


def load(run_dir):
    records = []
    for path in sorted(glob.glob(os.path.join(run_dir, "cat*.jsonl"))):
        with open(path) as handle:
            records.extend(json.loads(line) for line in handle if line.strip())
    return records


def _fmt(value):
    if isinstance(value, float):
        return f"{value:.4g}"
    text = json.dumps(value, ensure_ascii=False) if isinstance(value, (dict, list)) else str(value)
    return text if len(text) < 160 else text[:157] + "..."


def write_report(run_dir):
    from src.tools.registry import TOOLS, auto_import_tools
    auto_import_tools()
    records = load(run_dir)
    tools = sorted(set(TOOLS) | {r["tool"] for r in records})
    cats = [c for c in CATEGORY_TITLES if any(r["category"] == c for r in records)]

    def count(tool, cat=None, status=None):
        return sum(1 for r in records if r["tool"] == tool and (cat is None or r["category"] == cat)
                   and (status is None or r["status"] == status))

    lines = ["# Tool verification report", "", f"Run directory: `{run_dir}`", ""]
    totals = {s: sum(r["status"] == s for r in records) for s in ("pass", "fail", "info", "error")}
    lines += [f"**{len(records)} metrics**: {totals['pass']} pass, {totals['fail']} fail, {totals['error']} error, "
              f"{totals['info']} info (reported against a stated range, no threshold).", ""]

    header = "| tool | " + " | ".join(CATEGORY_TITLES[c] for c in cats) + " | verdict |"
    lines += ["## Per tool (pass/fail/error, info in brackets)", "", header,
              "|" + "---|" * (len(cats) + 2)]
    summary = {}
    for tool in tools:
        cells = []
        for c in cats:
            p, f, e, i = (count(tool, c, s) for s in ("pass", "fail", "error", "info"))
            cells.append("" if p + f + e + i == 0 else f"{p}/{f}/{e}" + (f" ({i})" if i else ""))
        fails, errors = count(tool, status="fail"), count(tool, status="error")
        tested = count(tool) > 0
        verdict = "not tested" if not tested else ("ERROR" if errors else ("FAIL" if fails else "pass"))
        summary[tool] = {"verdict": verdict, "pass": count(tool, status="pass"), "fail": fails, "error": errors,
                         "info": count(tool, status="info")}
        registered = "" if tool in TOOLS else " *(not a registered tool)*"
        lines.append(f"| `{tool}`{registered} | " + " | ".join(cells) + f" | {verdict} |")

    problems = [r for r in records if r["status"] in ("fail", "error")]
    lines += ["", f"## Failures and errors ({len(problems)})", ""]
    for r in problems:
        lines.append(f"- **{r['tool']}** [{r['category']}] {r['case']} — {r['metric']}: "
                     f"value {_fmt(r.get('value'))}, threshold {r.get('op', '')} {_fmt(r.get('threshold'))}"
                     + (f". {r['note'][:400]}" if r.get("note") else ""))
    infos = [r for r in records if r["status"] == "info"]
    lines += ["", f"## Reported values ({len(infos)})", ""]
    for r in infos:
        lines.append(f"- **{r['tool']}** [{r['category']}] {r['metric']}: {_fmt(r.get('value'))}"
                     + (f" — {r['note'][:300]}" if r.get("note") else ""))

    path = os.path.join(run_dir, "report.md")
    with open(path, "w") as handle:
        handle.write("\n".join(lines) + "\n")
    with open(os.path.join(run_dir, "summary.json"), "w") as handle:
        json.dump({"totals": totals, "tools": summary}, handle, indent=2)
    return path


if __name__ == "__main__":
    print(write_report(sys.argv[1]))
