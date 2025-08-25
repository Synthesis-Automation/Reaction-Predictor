#!/usr/bin/env python3
"""
Auto-generate an evaluation gold CSV from data/reaction_dataset/*.
Produces columns: ReactionID, ReactionType, ReactantSMILES, ProductSMILES, Ligand, Base, Solvent
- Handles CSV/TSV
- Maps flexible fields (Ligand list-like, Base via Reagent/Role, Solvent via Solvent/SOLName)
- Limits number of rows (default 200), keeping only rows with at least one of Ligand/Base/Solvent
"""
from __future__ import annotations

import os, sys, csv, json
from typing import List, Dict, Optional

_ROOT = os.path.abspath(os.path.dirname(__file__) + os.sep + "..")


def _ensure_dir(p: str):
    os.makedirs(p, exist_ok=True)


def _first_from_listlike(s: Optional[str]) -> str:
    if not s:
        return ''
    txt = str(s).strip()
    if not txt:
        return ''
    if txt.startswith('[') and txt.endswith(']'):
        try:
            arr = json.loads(txt)
            for x in arr:
                v = str(x).strip()
                if v:
                    return v
        except Exception:
            pass
    parts = [p.strip() for p in txt.strip('[]').replace('"','').replace("'",'').split(',') if p.strip()]
    return parts[0] if parts else txt


def _derive_base(row: Dict) -> str:
    base = row.get('Base') or ''
    if base:
        return _first_from_listlike(base)
    role = (row.get('ReagentRole') or '').lower()
    if 'base' in role:
        return _first_from_listlike(row.get('RGTName') or row.get('Reagent') or '')
    txt = ' '.join([str(row.get('Reagent') or ''), str(row.get('RGTName') or '')]).lower()
    for tok in ['k2co3','cs2co3','k3po4','kotbu','naotbu','koh','triethylamine','et3n','dipea','dbu']:
        if tok in txt:
            return tok
    return ''


def _scan_dataset(data_dir: str) -> List[Dict]:
    pool: List[Dict] = []
    for fname in os.listdir(data_dir):
        if not fname.lower().endswith(('.csv', '.tsv')):
            continue
        path = os.path.join(data_dir, fname)
        delim = '\t' if fname.lower().endswith('.tsv') else ','
        try:
            with open(path, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f, delimiter=delim)
                for row in reader:
                    pool.append(row)
        except Exception:
            continue
    return pool


def make_gold(limit: int = 200, out_path: str = os.path.join(_ROOT, 'data', 'eval_gold_auto.csv')) -> str:
    data_dir = os.path.join(_ROOT, 'data', 'reaction_dataset')
    pool = _scan_dataset(data_dir)
    rows: List[Dict] = []
    for r in pool:
        lig = _first_from_listlike(r.get('Ligand') or '')
        base = _derive_base(r)
        solv = _first_from_listlike(r.get('Solvent') or r.get('SOLName') or '')
        # keep rows that have at least one label to evaluate
        if not (lig or base or solv):
            continue
        rows.append({
            'ReactionID': r.get('ReactionID') or '',
            'ReactionType': r.get('ReactionType') or '',
            'ReactantSMILES': r.get('ReactantSMILES') or '',
            'ProductSMILES': r.get('ProductSMILES') or '',
            'Ligand': lig,
            'Base': base,
            'Solvent': solv,
        })
        if len(rows) >= limit:
            break
    out_dir = os.path.dirname(out_path)
    _ensure_dir(out_dir)
    with open(out_path, 'w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['ReactionID','ReactionType','ReactantSMILES','ProductSMILES','Ligand','Base','Solvent'])
        w.writeheader()
        w.writerows(rows)
    return out_path


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default=os.path.join(_ROOT, 'data', 'eval_gold_auto.csv'))
    ap.add_argument('--limit', type=int, default=200)
    args = ap.parse_args()
    path = make_gold(args.limit, args.out)
    print(path)
