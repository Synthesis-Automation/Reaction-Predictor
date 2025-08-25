#!/usr/bin/env python3
"""
Evaluate agent (ligand/solvent/base) recommendations for a set of reactions.
Outputs summary JSON and per-reaction CSV under reports/.

Phase 2 scope: agents only, Top-K metrics with bootstrap 95% CIs.
Flexible schema mapper supports updated dataset formats (Ullmann/Buchwald/Amide).
"""
from __future__ import annotations

import os, sys, json, csv, time
from typing import List, Dict, Tuple, Optional
import random

_ROOT = os.path.abspath(os.path.dirname(__file__) + os.sep + "..")
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from enhanced_recommendation_engine import create_recommendation_engine


def _ensure_dir(p: str):
    os.makedirs(p, exist_ok=True)


def _load_gold(path: str) -> List[Dict]:
        """Load gold CSV/TSV with flexible columns.
        Accepts:
            - ReactionSMILES or (ReactantSMILES, ProductSMILES)
            - Ligand as a scalar or list-like string; uses first item
            - Base from columns: Base, or derived from (Reagent/ReagentRole/RGTName)
            - Solvent from Solvent or SOLName
        """
        out: List[Dict] = []
        delim = '\t' if path.lower().endswith('.tsv') else ','
        with open(path, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f, delimiter=delim)
                for row in reader:
                        out.append(row)
        return out


def _canon(s: str) -> str:
    return (s or '').strip().casefold()


def _topk_contains(cands: List[str], gold: str, k: int) -> bool:
    g = _canon(gold)
    if not g:
        return False
    return any(_canon(x) == g for x in cands[:k])


def _extract_names(items: List[Dict], key: str) -> List[str]:
    names = []
    for it in items or []:
        nm = it.get(key)
        if nm:
            names.append(str(nm))
    return names


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
    # fallback split by comma
    parts = [p.strip() for p in txt.strip('[]').replace('"','').replace("'",'').split(',') if p.strip()]
    return parts[0] if parts else txt


def _derive_base(row: Dict) -> str:
    base = row.get('Base') or ''
    if base:
        return _first_from_listlike(base)
    # Derive from reagent columns if present
    # If roles included, pick where role contains 'base'
    role = (row.get('ReagentRole') or '').lower()
    if 'base' in role:
        return _first_from_listlike(row.get('RGTName') or row.get('Reagent') or '')
    # Fallback: scan tokens for common base patterns
    txt = ' '.join([str(row.get('Reagent') or ''), str(row.get('RGTName') or '')]).lower()
    for tok in ['k2co3','cs2co3','k3po4','kotbu','naotbu','koh','triethylamine','et3n','dipea','dbu']:
        if tok in txt:
            return tok
    return ''


def _rxn_smiles_from_row(row: Dict) -> str:
    rs = row.get('ReactionSMILES')
    if rs:
        return str(rs)
    return f"{row.get('ReactantSMILES') or ''}>>{row.get('ProductSMILES') or ''}"


def _bootstrap_ci(values: List[int], total: int, n_iter: int = 1000, alpha: float = 0.05) -> Tuple[float,float]:
    if total <= 0:
        return (0.0, 0.0)
    acc = sum(values)/total
    if total < 2:
        return (acc, acc)
    random.seed(1337)
    props: List[float] = []
    for _ in range(n_iter):
        sample = [values[random.randrange(0, len(values))] for _ in range(len(values))]
        props.append(sum(sample)/len(sample) if sample else 0.0)
    props.sort()
    lo = props[int((alpha/2)*len(props))]
    hi = props[int((1-alpha/2)*len(props))-1]
    return (lo, hi)


def evaluate(gold_csv: str, k: int = 3) -> Dict:
    gold = _load_gold(gold_csv)
    eng = create_recommendation_engine()
    results: List[Dict] = []
    role_hits = { 'ligand': 0, 'base': 0, 'solvent': 0 }
    role_vecs: Dict[str, List[int]] = { 'ligand': [], 'base': [], 'solvent': [] }
    total = 0

    for row in gold:
        total += 1
        rxn = _rxn_smiles_from_row(row)
        rtype = row.get('ReactionType') or 'Auto-detect'
        recs = eng.get_recommendations(rxn, rtype)
        ligs = _extract_names(recs.get('ligand_recommendations') or [], 'ligand')
        bases = _extract_names(recs.get('base_recommendations') or [], 'base')
        solvs = _extract_names(recs.get('solvent_recommendations') or [], 'solvent')

        lig_gold = _first_from_listlike(row.get('Ligand') or '')
        base_gold = _derive_base(row)
        solv_gold = _first_from_listlike(row.get('Solvent') or row.get('SOLName') or '')

        lig_hit = _topk_contains(ligs, lig_gold, k)
        base_hit = _topk_contains(bases, base_gold, k)
        solv_hit = _topk_contains(solvs, solv_gold, k)

        role_hits['ligand'] += 1 if lig_hit else 0
        role_hits['base'] += 1 if base_hit else 0
        role_hits['solvent'] += 1 if solv_hit else 0
        role_vecs['ligand'].append(1 if lig_hit else 0)
        role_vecs['base'].append(1 if base_hit else 0)
        role_vecs['solvent'].append(1 if solv_hit else 0)

        results.append({
            'reaction_id': row.get('ReactionID'),
            'reaction_type': rtype,
            'ligand_gold': lig_gold,
            'base_gold': base_gold,
            'solvent_gold': solv_gold,
            'ligand_topk': ligs[:k],
            'base_topk': bases[:k],
            'solvent_topk': solvs[:k],
            'ligand_hit': lig_hit,
            'base_hit': base_hit,
            'solvent_hit': solv_hit,
        })

    # Bootstrap CIs per role
    cis = { r: _bootstrap_ci(role_vecs[r], total) for r in role_vecs }

    summary = {
        'k': k,
        'total': total,
        'topk_accuracy': { r: (role_hits[r] / total if total else 0.0) for r in role_hits },
        'topk_accuracy_ci95': { r: {'lo': cis[r][0], 'hi': cis[r][1]} for r in cis },
        'generated_at': time.strftime('%Y-%m-%dT%H:%M:%S')
    }

    out_dir = os.path.join(_ROOT, 'reports')
    _ensure_dir(out_dir)
    with open(os.path.join(out_dir, 'evaluation_summary.json'), 'w', encoding='utf-8') as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    with open(os.path.join(out_dir, 'evaluation_detailed.csv'), 'w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(results[0].keys()) if results else [])
        w.writeheader()
        w.writerows(results)

    return summary


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--gold-csv', required=True, help='Path to gold CSV with Ligand/Base/Solvent columns')
    ap.add_argument('--k', type=int, default=3)
    args = ap.parse_args()
    s = evaluate(args.gold_csv, args.k)
    print(json.dumps(s, ensure_ascii=False))
