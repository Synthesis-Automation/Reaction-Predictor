#!/usr/bin/env python3
"""
Quick smoke test using the default Python environment (no special env needed).
- Calls the enhanced engine directly (same logic as CLI) and builds the export payload
- Prints a compact JSON: { meta.providers, top_conditions_len }
- Exits 0 on success (>=1 top condition), non-zero otherwise
"""
from __future__ import annotations

import sys, json
from typing import Optional

# Local imports (repo-root is added by modules themselves)
from enhanced_recommendation_engine import create_recommendation_engine
from prediction_export import build_export_payload


def run(smiles: str, selected_type: str) -> int:
    eng = create_recommendation_engine()
    recs = eng.get_recommendations(smiles, selected_type)
    result = {
        'analysis_type': recs.get('analysis_type'),
        'reaction_smiles': smiles,
        'reaction_type': selected_type,
        'status': recs.get('status', 'success' if 'error' not in recs else 'failed'),
        'recommendations': {
            'reaction_type': recs.get('reaction_type'),
            'ligand_recommendations': recs.get('ligand_recommendations', []),
            'solvent_recommendations': recs.get('solvent_recommendations', []),
            'base_recommendations': recs.get('base_recommendations', []),
            'combined_conditions': recs.get('combined_conditions', []),
            'dataset_info': recs.get('dataset_info', {}),
        },
    }
    if recs.get('providers'):
        result['providers'] = recs.get('providers')
    export = build_export_payload(result)
    compact = {
        'meta': { 'providers': (export.get('meta') or {}).get('providers') },
        'top_conditions_len': len(export.get('top_conditions') or [])
    }
    print(json.dumps(compact, ensure_ascii=False))
    return 0 if compact['top_conditions_len'] >= 1 else 2


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--smiles', default='Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1')
    ap.add_argument('--type', default='C-N Coupling - Ullmann')
    args = ap.parse_args()
    raise SystemExit(run(args.smiles, args.type))
