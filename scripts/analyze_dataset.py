#!/usr/bin/env python3
from __future__ import annotations
"""
Dataset analytics: compute simple frequency priors for reagents by reaction type.

Outputs JSON summaries that EnhancedRecommendationEngine can load as priors.

Currently optimized for Ullmann: extracts bases, solvents, and ligands from
CSV files under data/reaction_dataset/*.csv and writes to data/analytics/Ullmann/latest.json.

Usage examples:
  python scripts/analyze_dataset.py --reaction-type Ullmann
  python scripts/analyze_dataset.py --reaction-type Ullmann --top-limit 5 --co-limit 20
"""

import argparse
import csv
import glob
import json
import os
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from typing import Dict, Iterable, List, Tuple

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))


def _parse_list_field(val: str) -> List[str]:
    """Parse fields that may contain a JSON-ish list string like ["A", "B"].
    Falls back to splitting by comma. Returns trimmed non-empty items.
    """
    if val is None:
        return []
    s = str(val).strip()
    if not s:
        return []
    # Try JSON first
    if s.startswith('[') and s.endswith(']'):
        try:
            arr = json.loads(s)
            out = []
            for it in arr:
                t = str(it).strip()
                if t:
                    out.append(t)
            return out
        except Exception:
            pass
    # Fallback: strip brackets/quotes and split by comma
    s2 = s.strip('[]').replace('"', '').replace("'", '')
    parts = [p.strip() for p in s2.split(',') if p.strip()]
    if parts:
        return parts
    return [s] if s else []


def _canon_text(x: str) -> str:
    return (x or '').strip()


def _canon_solvent(x: str) -> str:
    s = (x or '').strip()
    # Normalize a few common synonyms/abbreviations
    low = s.lower().replace(' ', '')
    mapping = {
        'n-methyl-2-pyrrolidone': 'N-Methyl-2-pyrrolidone',
        'nmp': 'N-Methyl-2-pyrrolidone',
        'dimethylformamide': 'Dimethylformamide',
        'dmf': 'DMF',
        'dimethylsulfoxide': 'Dimethyl sulfoxide',
        'dmso': 'DMSO',
        '1,4-dioxane': '1,4-Dioxane',
        'dioxane': '1,4-Dioxane',
        'toluene': 'Toluene',
        'thf': 'THF',
    }
    if low in mapping:
        return mapping[low]
    # If field appears as CAS-only or unresolved, keep original
    return _canon_text(s)


def _canon_base(x: str) -> str:
    s = (x or '').strip()
    # Favor the formula inside parentheses if present, otherwise normalize common bases
    if '(' in s and ')' in s:
        inner = s[s.rfind('(')+1:s.rfind(')')].strip()
        if inner:
            s = inner
    low = s.lower().replace(' ', '')
    mapping = {
        'k2co3': 'K2CO3',
        'cesiumcarbonate': 'Cs2CO3',
        'cs2co3': 'Cs2CO3',
        'koh': 'KOH',
        'naoh': 'NaOH',
        'triethylamine': 'Et3N',
        'et3n': 'Et3N',
        'k3po4': 'K3PO4',
        'ko tbu': 'KOtBu', 'kotbu': 'KOtBu', 'tbuok': 'KOtBu', 'po(tbu)k': 'KOtBu',
        'naotbu': 'NaOtBu',
        'na2co3': 'Na2CO3',
    }
    if low in mapping:
        return mapping[low]
    return _canon_text(s)


def _canon_ligand(x: str) -> str:
    # Keep ligand names as-is, just trim
    return _canon_text(x)


def _extract_solvents(row: Dict[str, str]) -> List[str]:
    raw = row.get('Solvent') or row.get('SOLName') or ''
    items = _parse_list_field(raw)
    if not items and raw:
        items = [str(raw)]
    return [_canon_solvent(it) for it in items if it]


def _extract_bases(row: Dict[str, str]) -> List[str]:
    # Try Base column first; then ReagentRaw or RGTName and filter base tokens
    items: List[str] = []
    for key in ('Base', 'ReagentRaw', 'RGTName'):
        raw = row.get(key) or ''
        parts = _parse_list_field(raw)
        for p in parts:
            items.append(p)
    # Filter to typical base names
    base_like = []
    tokens = (
        'k2co3','cs2co3','k3po4','kotbu','tbuok','ko tbu','naotbu','na2co3','koh','naoh','triethylamine','et3n','carbonate','phosphate','hydroxide','t‑butoxide','tert-butoxide','t-butoxide'
    )
    for it in items:
        low = it.lower().replace(' ', '')
        if any(tok in low for tok in tokens):
            base_like.append(_canon_base(it))
    # If Base column gave nothing but there is a dedicated Base field
    if not base_like and row.get('Base'):
        for p in _parse_list_field(row['Base']):
            base_like.append(_canon_base(p))
    return [b for b in base_like if b]


def _extract_ligands(row: Dict[str, str]) -> List[str]:
    raw = row.get('Ligand') or ''
    items = _parse_list_field(raw)
    return [_canon_ligand(it) for it in items if it]


def analyze_reaction_type(reaction_type: str, dataset_glob: str) -> Dict:
    files = sorted(glob.glob(os.path.join(ROOT, dataset_glob)))
    if not files:
        raise FileNotFoundError(f"No dataset files found under pattern: {dataset_glob}")

    solvents = Counter()
    bases = Counter()
    ligands = Counter()
    co_base_solvent = Counter()  # simple co-occurrence of base-solvent per row
    rows_total = 0
    files_considered: List[str] = []

    for path in files:
        try:
            with open(path, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f)
                matched_any = False
                for row in reader:
                    rtype = (row.get('ReactionType') or '').strip()
                    if not rtype or rtype.lower() != reaction_type.lower():
                        continue
                    matched_any = True
                    rows_total += 1
                    s_list = _extract_solvents(row)
                    b_list = _extract_bases(row)
                    l_list = _extract_ligands(row)
                    if s_list:
                        solvents.update(s_list)
                    if b_list:
                        bases.update(b_list)
                    if l_list:
                        ligands.update(l_list)
                    # co-occurrence: each base with each solvent in the same row
                    for b in set(b_list):
                        for s in set(s_list):
                            co_base_solvent[(b, s)] += 1
                if matched_any:
                    files_considered.append(os.path.relpath(path, ROOT))
        except Exception as e:
            # Skip bad file, continue
            continue

    summary = {
        'reaction_type': reaction_type,
        'generated_at': datetime.now(timezone.utc).isoformat(),
        'source_files': files_considered,
        'totals': {
            'rows': rows_total,
            'unique_solvents': len(solvents),
            'unique_bases': len(bases),
            'unique_ligands': len(ligands),
        },
        'top': {
            'solvents': [],
            'bases': [],
            'ligands': [],
            'co_base_solvent': [],
        }
    }

    return summary, solvents, bases, ligands, co_base_solvent


def _top_list(counter: Counter, total_rows: int, limit: int) -> List[Dict[str, object]]:
    out: List[Dict[str, object]] = []
    if total_rows <= 0:
        return out
    for name, cnt in counter.most_common(limit if limit > 0 else None):
        pct = float(cnt) / float(total_rows)
        out.append({'name': name, 'count': int(cnt), 'pct': round(pct, 4)})
    return out


def main(argv: List[str]) -> int:
    p = argparse.ArgumentParser(description='Compute reagent priors from datasets.')
    p.add_argument('--reaction-type', required=True, help='Reaction type to analyze (e.g., Ullmann)')
    p.add_argument('--dataset-glob', default='data/reaction_dataset/*.csv', help='Glob for dataset CSVs relative to repo root')
    p.add_argument('--top-limit', type=int, default=10, help='Top-N entries per category')
    p.add_argument('--co-limit', type=int, default=30, help='Top-N co-occurrence pairs')
    args = p.parse_args(argv)

    reaction_type = args.reaction_type.strip()

    summary, solvents, bases, ligands, co = analyze_reaction_type(reaction_type, args.dataset_glob)

    # Populate top lists
    rows_total = summary['totals']['rows']
    summary['top']['solvents'] = _top_list(solvents, rows_total, args.top_limit)
    summary['top']['bases'] = _top_list(bases, rows_total, args.top_limit)
    summary['top']['ligands'] = _top_list(ligands, rows_total, args.top_limit)
    # Co-occurrence pairs
    co_items: List[Tuple[Tuple[str, str], int]] = co.most_common(args.co_limit if args.co_limit > 0 else None)
    summary['top']['co_base_solvent'] = [
        {'base': b, 'solvent': s, 'count': int(cnt)} for (b, s), cnt in co_items
    ]

    # Write to data/analytics/<reaction_type>/latest.json and timestamped file
    out_dir = os.path.join(ROOT, 'data', 'analytics', reaction_type)
    os.makedirs(out_dir, exist_ok=True)
    latest_path = os.path.join(out_dir, 'latest.json')
    ts = datetime.now().strftime('%Y%m%d-%H%M%S')
    stamp_path = os.path.join(out_dir, f'summary_{ts}.json')
    with open(latest_path, 'w', encoding='utf-8') as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    with open(stamp_path, 'w', encoding='utf-8') as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)

    print(f"Wrote analytics summary to: {os.path.relpath(latest_path, ROOT)}")
    print(f"Rows analyzed: {rows_total}")
    print('Top solvents:', [x['name'] for x in summary['top']['solvents'][:5]])
    print('Top bases:', [x['name'] for x in summary['top']['bases'][:5]])
    print('Top ligands:', [x['name'] for x in summary['top']['ligands'][:5]])
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
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
    parser.add_argument("--reaction-type", default="Ullmann", help="Reaction type to analyze (currently supports Ullmann)")
    parser.add_argument("--out-dir", default=REPO_ROOT, help="Output base directory (repo root by default)")
    parser.add_argument("--no-json", action="store_true", help="Do not write JSON outputs")
    parser.add_argument("--no-csv", action="store_true", help="Do not write CSV outputs")
    parser.add_argument("--top-limit", type=int, default=None, help="Limit number of top items per category in outputs")
    parser.add_argument("--co-limit", type=int, default=None, help="Limit number of co-occurrence rows per pair list")
    args = parser.parse_args()

    rt = args.reaction_type
    if rt.lower().startswith("ullmann"):
        out = run_and_export_ullmann(
            args.out_dir,
            write_json=not args.no_json,
            write_csv=not args.no_csv,
            top_limit=args.top_limit,
            co_limit=args.co_limit,
        )
        print(out)
    else:
        print(f"Unsupported reaction type for Milestone 1: {rt}")
        sys.exit(2)


if __name__ == "__main__":
    main()
