from __future__ import annotations

from prediction_export import build_export_payload


def test_export_payload_simplified_schema():
    # Minimal fake result resembling GUI/CLI result structure
    result = {
        'analysis_type': 'enhanced',
        'reaction_smiles': 'Clc1ccncc1.NCC>>CCNc1ccncc1',
        'reaction_type': 'C-N Coupling - Ullmann',
        'status': 'ok',
        'recommendations': {
            'reaction_type': 'Ullmann',
            'dataset_info': {'ligands_available': 10, 'solvents_available': 5, 'reaction_types_supported': ['Ullmann'], 'analytics_loaded': False},
            'combined_conditions': [
                {
                    'ligand': 'L-Proline',
                    'solvent': 'DMF',
                    'solvent_abbreviation': 'DMF',
                    'ligand_compatibility': 0.8,
                    'solvent_compatibility': 0.9,
                    'combined_score': 0.85,
                    'recommendation_confidence': 'High',
                    'typical_conditions': {'temperature': '120 °C', 'time': '6 h', 'atmosphere': 'inert', 'base': 'K3PO4', 'catalyst_loading': '10 mol% Cu'},
                }
            ],
        },
    }

    export = build_export_payload(result)

    # Schema checks
    assert 'recommendations' not in export  # simplified output should not include this block
    for key in ('meta', 'input', 'detection', 'dataset', 'top_conditions', 'related_reactions'):
        assert key in export

    # top_conditions should be a short list with chemicals/conditions
    assert isinstance(export['top_conditions'], list)
    if export['top_conditions']:
        tc0 = export['top_conditions'][0]
        assert 'chemicals' in tc0 and 'conditions' in tc0
