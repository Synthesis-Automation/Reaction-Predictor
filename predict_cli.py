#!/usr/bin/env python3
"""
Command-line predictor that shares the same enhanced recommendation engine and export builder as the GUI.

Usage (PowerShell example):
  python predict_cli.py '{"reaction_smiles":"Brc1cc...","selected_reaction_type":"C-N Coupling - Ullmann"}'

Or read from stdin:
  Get-Content input.json | python predict_cli.py
"""
from __future__ import annotations

import sys
import os
import json

# Ensure project root on path
_HERE = os.path.dirname(__file__)
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

from enhanced_recommendation_engine import create_recommendation_engine
from prediction_export import build_export_payload


def _load_input() -> dict:
    # Try argv first
    if len(sys.argv) > 1 and sys.argv[1].strip():
        arg = sys.argv[1]
        try:
            return json.loads(arg)
        except Exception:
            pass
    # Fallback to stdin
    try:
        data = sys.stdin.read()
        if data.strip():
            return json.loads(data)
    except Exception:
        pass
    return {}


def main() -> int:
    payload = _load_input()
    reaction_smiles = payload.get('reaction_smiles') or ''
    selected_type = payload.get('selected_reaction_type') or 'Auto-detect'

    engine = create_recommendation_engine()
    recs = engine.get_recommendations(reaction_smiles, selected_type)

    # Normalize to GUI-like result to reuse export builder
    result = {
        'analysis_type': recs.get('analysis_type'),
        'reaction_smiles': reaction_smiles,
        'reaction_type': selected_type,
        'status': recs.get('status', 'success' if 'error' not in recs else 'failed'),
        'recommendations': {
            'reaction_type': recs.get('reaction_type'),
            'ligand_recommendations': recs.get('ligand_recommendations', []),
            'solvent_recommendations': recs.get('solvent_recommendations', []),
            'base_recommendations': recs.get('base_recommendations', []),
            'combined_conditions': recs.get('combined_conditions', []),
            'property_based_alternatives': recs.get('property_based_alternatives', {}),
            'dataset_info': recs.get('dataset_info', {}),
            'related_reactions': recs.get('related_reactions', []),
        }
    }

    export = build_export_payload(result)
    sys.stdout.write(json.dumps(export, ensure_ascii=False))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
