from __future__ import annotations

import argparse
import os
import sys

# Allow running from repo root without installing
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from analytics.aggregate import run_and_export_ullmann


def main():
    parser = argparse.ArgumentParser(description="Analyze reaction datasets and export analytics.")
    parser.add_argument("--reaction-type", default="Ullmann", help="Reaction type to analyze (Milestone 1 supports Ullmann)")
    parser.add_argument("--out-dir", default=REPO_ROOT, help="Output base directory (repo root by default)")
    args = parser.parse_args()

    rt = args.reaction_type
    if rt.lower().startswith("ullmann"):
        out = run_and_export_ullmann(args.out_dir)
        print(out)
    else:
        print(f"Unsupported reaction type for Milestone 1: {rt}")
        sys.exit(2)


if __name__ == "__main__":
    main()
