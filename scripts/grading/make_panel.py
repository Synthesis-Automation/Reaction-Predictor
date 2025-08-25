#!/usr/bin/env python3
"""
Create a 20-reaction grading panel stratified across reaction families.
Outputs CSV and HTML under reports/grading/.
"""
from __future__ import annotations

import os, sys, csv, json, random, time
from typing import List, Dict

_ROOT = os.path.abspath(os.path.dirname(__file__) + os.sep + ".." + os.sep + "..")
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from enhanced_recommendation_engine import create_recommendation_engine

random.seed(42)

FAMILIES = [
    ('Ullmann', 5),
    ('Amide Formation', 5),
    ('Cross-Coupling', 5),
    ('Other', 5),
]


def _ensure_dir(p: str):
    os.makedirs(p, exist_ok=True)


def _load_pool(data_dir: str) -> List[Dict]:
    pool: List[Dict] = []
    for fname in os.listdir(data_dir):
        if not fname.lower().endswith(('.csv', '.tsv')):
            continue
        path = os.path.join(data_dir, fname)
        try:
            with open(path, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f, delimiter='\t') if fname.lower().endswith('.tsv') else csv.DictReader(f)
                for row in reader:
                    pool.append(row)
        except Exception:
            continue
    return pool


def _rxn_smiles(row: Dict) -> str:
    return f"{row.get('ReactantSMILES') or ''}>>{row.get('ProductSMILES') or ''}"


def _family_of(row: Dict) -> str:
    rt = (row.get('ReactionType') or '').strip()
    if not rt:
        return 'Other'
    low = rt.lower()
    if 'ullmann' in low:
        return 'Ullmann'
    if 'amide' in low or 'amidation' in low:
        return 'Amide Formation'
    if any(x in low for x in ['buchwald', 'suzuki', 'heck', 'sonogashira', 'stille', 'negishi', 'chan-lam', 'cross-coupling']):
        return 'Cross-Coupling'
    return 'Other'


def make_panel() -> Dict:
    data_dir = os.path.join(_ROOT, 'data', 'reaction_dataset')
    pool = _load_pool(data_dir)
    by_family: Dict[str, List[Dict]] = { k: [] for k,_ in FAMILIES }
    for row in pool:
        fam = _family_of(row)
        if fam not in by_family:
            continue
        by_family[fam].append(row)

    panel: List[Dict] = []
    for fam, n in FAMILIES:
        cand = by_family.get(fam, [])
        if len(cand) <= n:
            panel.extend(cand)
        else:
            panel.extend(random.sample(cand, n))

    # Build rubric entries with top-3 per role
    eng = create_recommendation_engine()
    out_rows: List[Dict] = []
    for r in panel:
        rxn = _rxn_smiles(r)
        rt = r.get('ReactionType') or 'Auto-detect'
        recs = eng.get_recommendations(rxn, rt)
        ligs = [x.get('ligand') for x in (recs.get('ligand_recommendations') or [])[:3]]
        bases = [x.get('base') for x in (recs.get('base_recommendations') or [])[:3]]
        solvs = [x.get('solvent') for x in (recs.get('solvent_recommendations') or [])[:3]]
        out_rows.append({
            'reaction_smiles': rxn,
            'reaction_type': rt,
            'ligand_top3': ', '.join([x for x in ligs if x]),
            'base_top3': ', '.join([x for x in bases if x]),
            'solvent_top3': ', '.join([x for x in solvs if x]),
            'grade_ligand': '',
            'grade_base': '',
            'grade_solvent': '',
            'notes': ''
        })

    out_dir = os.path.join(_ROOT, 'reports', 'grading')
    _ensure_dir(out_dir)
    csv_path = os.path.join(out_dir, 'grading_panel.csv')
    html_path = os.path.join(out_dir, 'grading_panel.html')

    with open(csv_path, 'w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()) if out_rows else [])
        w.writeheader()
        w.writerows(out_rows)

    # Simple HTML
    html = ['<html><head><meta charset="utf-8"><title>Grading Panel</title></head><body>']
    html.append('<h2>20-Reaction Grading Panel</h2>')
    html.append('<table border="1" cellspacing="0" cellpadding="4">')
    if out_rows:
        # header
        html.append('<tr>' + ''.join(f'<th>{k}</th>' for k in out_rows[0].keys()) + '</tr>')
        for r in out_rows:
            html.append('<tr>' + ''.join(f'<td>{(r.get(k) or "")}</td>' for k in out_rows[0].keys()) + '</tr>')
    html.append('</table></body></html>')
    with open(html_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(html))

    return { 'csv': csv_path, 'html': html_path, 'count': len(out_rows), 'generated_at': time.strftime('%Y-%m-%dT%H:%M:%S') }


if __name__ == '__main__':
    print(json.dumps(make_panel(), ensure_ascii=False))
