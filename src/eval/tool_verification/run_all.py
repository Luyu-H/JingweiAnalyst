"""
Run every verification category, one process per category, then write the report.

    python -u -m src.eval.tool_verification.run_all [--out results/tool_verification/<run>] [--only cat1 cat3]
"""
import argparse
import os
import subprocess
import sys
import time

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
sys.path.insert(0, REPO)

from src.eval.tool_verification.run_category import CATEGORIES  # noqa: E402
from src.eval.tool_verification.summarize import write_report  # noqa: E402


def log(path, message):
    line = f"{time.strftime('%Y-%m-%d %H:%M:%S')} {message}"
    print(line, flush=True)
    with open(path, "a") as handle:
        handle.write(line + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default=os.path.join(REPO, "results", "tool_verification",
                                                      time.strftime("%Y%m%d_%H%M%S")))
    parser.add_argument("--only", nargs="*", help="category name prefixes, e.g. cat1 cat3")
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)
    status = os.path.join(args.out, "status.log")
    categories = [c for c in CATEGORIES if not args.only or any(c.startswith(o) for o in args.only)]
    log(status, f"START {categories} -> {args.out}")

    for category in categories:
        for attempt in (1, 2):
            started = time.time()
            jsonl = os.path.join(args.out, f"{category}.jsonl")
            if os.path.exists(jsonl):
                os.remove(jsonl)            # a retry starts the category afresh
            with open(os.path.join(args.out, f"{category}.log"), "w") as out:
                proc = subprocess.run(
                    [sys.executable, "-u", "-m", "src.eval.tool_verification.run_category", category, "--out", args.out],
                    cwd=REPO, stdout=out, stderr=subprocess.STDOUT, env={**os.environ, "PYTHONPATH": REPO},
                )
            rc = proc.returncode
            log(status, f"END {category} attempt={attempt} rc={rc} ({time.time() - started:.0f}s)")
            if rc >= 0:                     # finished, whatever the test outcomes
                break
            log(status, f"{category} died on signal {-rc}; retrying once")

    report = write_report(args.out)
    log(status, f"DONE report={report}")


if __name__ == "__main__":
    main()
