"""
Run one verification category in this process.

    python -u -m src.eval.tool_verification.run_category cat1_analytic --out results/tool_verification/<run>
    python -u -m src.eval.tool_verification.run_category cat2_reference --out ... --only compute_pdo_index
"""
import argparse
import importlib
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..")))

from src.eval.tool_verification.common import run_cases  # noqa: E402

CATEGORIES = [
    "cat1_analytic",
    "cat2_reference",
    "cat3_products",
    "cat4_magnitude",
    "cat5_consistency",
    "cat6_montecarlo",
    "cat7_safeguards",
]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("category", choices=CATEGORIES)
    parser.add_argument("--out", required=True)
    parser.add_argument("--only", nargs="*", help="substrings of 'tool/case' to run")
    args = parser.parse_args()

    module_name = f"src.eval.tool_verification.{args.category}"
    importlib.import_module(module_name)
    run_cases(args.category, module_name, args.out, only=args.only)


if __name__ == "__main__":
    main()
