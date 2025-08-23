"""
Shared export builder for prediction results.
Generates the structured JSON payload used by both GUI and CLI.
"""
from __future__ import annotations

import os
import time
import json
import csv
from typing import Optional


def build_export_payload(result: dict, related_reactions: Optional[list] = None) -> dict:
    """Convert internal result into a clean, stable JSON payload.

    Schema:
      {
        meta: { generated_at, analysis_type, status },
        input: { reaction_smiles, selected_reaction_type },
        detection: { reaction_type },
        dataset: { ligands_available, solvents_available, reaction_types_supported },
        recommendations: {
          combined: [ { ligand, ligand_compatibility, solvent, abbreviation, solvent_compatibility, combined_score, synergy_bonus, recommendation_confidence, typical_conditions, suggested_base } ],
          ligands: [ { ligand, compatibility_score, applications, reaction_suitability } ],
          solvents: [ { solvent, abbreviation, compatibility_score, applications, reaction_suitability } ],
          alternatives: { budget_friendly_ligands, low_boiling_solvents, green_solvents }
        },
        top_conditions: [ for top 3, a single 'chemicals' list with all required components (starting materials, metal precursor, ligand, base, solvent) and 'conditions' with time/temperature ],
        related_reactions: [ { reaction_smiles, yield, catalyst, ligand, solvent, temperature, time, similarity, reaction_id, reference } ]
      }
    """

    recs = result.get('recommendations', {}) or {}

    combined = []
    for c in recs.get('combined_conditions', []) or []:
        combined.append({
            'ligand': c.get('ligand'),
            'ligand_compatibility': c.get('ligand_compatibility'),
            'solvent': c.get('solvent'),
            'solvent_abbreviation': c.get('solvent_abbreviation'),
            'solvent_compatibility': c.get('solvent_compatibility'),
            'combined_score': c.get('combined_score'),
            'recommendation_confidence': c.get('recommendation_confidence'),
            'typical_conditions': c.get('typical_conditions', {}),
            'synergy_bonus': c.get('synergy_bonus', 0),
            'suggested_base': c.get('suggested_base'),
        })

    ligands = []
    for l in recs.get('ligand_recommendations', []) or []:
        ligands.append({
            'ligand': l.get('ligand'),
            'compatibility_score': l.get('compatibility_score'),
            'applications': l.get('applications'),
            'reaction_suitability': l.get('reaction_suitability'),
        })

    solvents = []
    for s in recs.get('solvent_recommendations', []) or []:
        solvents.append({
            'solvent': s.get('solvent'),
            'abbreviation': s.get('abbreviation'),
            'compatibility_score': s.get('compatibility_score'),
            'applications': s.get('applications'),
            'reaction_suitability': s.get('reaction_suitability'),
        })

    alternatives = {}
    alt = recs.get('property_based_alternatives', {}) or {}
    for key in ('budget_friendly_ligands', 'low_boiling_solvents', 'green_solvents'):
        if key in alt:
            alternatives[key] = alt[key]

    related = related_reactions or recs.get('related_reactions', []) or []

    # Helpers for top_conditions construction
    def _split_smiles(sm: str):
        try:
            lhs, rhs = (sm or '').split('>>', 1)
            reactants = [t for t in lhs.split('.') if t]
            product = rhs
            return reactants, product
        except Exception:
            return [], sm

    def _solvent_cas(name: Optional[str]):
        try:
            from reagents.solvent import create_solvent_dataframe  # type: ignore
            df = create_solvent_dataframe()
            if name and 'Solvent' in df.columns and 'CAS Number' in df.columns:
                row = df[df['Solvent'] == name]
                if not row.empty:
                    return row.iloc[0].get('CAS Number')
        except Exception:
            pass
        return None

    rxn_smiles = result.get('reaction_smiles', '')
    reactants, _ = _split_smiles(rxn_smiles)
    detected_type = recs.get('reaction_type') or result.get('reaction_type') or ''

    def _default_metal_precursor(rt: str):
        if isinstance(rt, str) and rt.lower() == 'ullmann':
            return {'name': 'CuI', 'cas': '7681-65-4', 'smiles': None, 'equivalents': None}
        else:
            return {'name': 'Pd(OAc)2', 'cas': '3375-31-3', 'smiles': None, 'equivalents': None}

    # Lightweight CAS lookup from data/cas_dictionary.csv for ligands/bases
    def _load_cas_map():
        cas_by_token = {}
        cas_by_name = {}
        try:
            data_dir = os.path.join(os.path.dirname(__file__), 'data')
            path = os.path.join(data_dir, 'cas_dictionary.csv')
            with open(path, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    cas = (row.get('CAS') or '').strip()
                    name = (row.get('Name') or '').strip()
                    token = (row.get('Token') or '').strip()
                    if cas:
                        if token:
                            cas_by_token[token.lower()] = cas
                        if name:
                            cas_by_name[name.lower()] = cas
        except Exception:
            pass
        aliases = {
            'k2co3': 'K2CO3',
            'cs2co3': 'Cs2CO3',
            'k3po4': 'K3PO4',
            'kotbu': 'KOtBu',
            'potassium tert-butoxide': 'KOtBu',
            'naotbu': 'NaOtBu',
            'dipea': 'DIPEA',
            'dbu': 'DBU',
            'pyridine': 'Pyridine',
            'xphos': 'XPhos',
            'sphos': 'SPhos',
            'ruphos': 'RuPhos',
            'brettphos': 'BrettPhos',
            'tbuxphos': 'tBuXPhos',
            'johnphos': 'JohnPhos',
            'xantphos': 'XantPhos',
            'dppe': 'DPPE',
            'dppf': 'DPPF',
            'binap': 'BINAP',
        }
        return cas_by_token, cas_by_name, aliases

    cas_by_token, cas_by_name, alias_map = _load_cas_map()

    def _lookup_cas(name: Optional[str]):
        if not name:
            return None
        key = name.strip().lower()
        if key in cas_by_token:
            return cas_by_token[key]
        if key in cas_by_name:
            return cas_by_name[key]
        if key in alias_map:
            alias = alias_map[key]
            ak = alias.lower()
            if ak in cas_by_token:
                return cas_by_token[ak]
            if ak in cas_by_name:
                return cas_by_name[ak]
        return None

    top_conditions = []
    for c in (recs.get('combined_conditions', []) or [])[:3]:
        conditions = c.get('typical_conditions', {}) or {}
        base_name = c.get('suggested_base') or conditions.get('base')
        chemicals = []
        # Starting materials
        for smi in reactants:
            chemicals.append({'name': None, 'cas': None, 'smiles': smi, 'equivalents': None, 'role': 'starting_material'})
        # Metal precursor
        mp = _default_metal_precursor(detected_type)
        chemicals.append({**mp, 'role': 'metal_precursor'})
        # Ligand
        lig_name = c.get('ligand')
        chemicals.append({'name': lig_name, 'cas': _lookup_cas(lig_name), 'smiles': None, 'equivalents': None, 'role': 'ligand'})
        # Base
        if base_name:
            chemicals.append({'name': base_name, 'cas': _lookup_cas(base_name), 'smiles': None, 'equivalents': 2.0, 'role': 'base'})
        # Solvent
        chemicals.append({'name': c.get('solvent'), 'abbreviation': c.get('solvent_abbreviation'), 'cas': _solvent_cas(c.get('solvent')), 'smiles': None, 'equivalents': None, 'role': 'solvent'})

        top_conditions.append({
            'reaction': {'smiles': rxn_smiles},
            'chemicals': chemicals,
            'conditions': {
                'temperature': conditions.get('temperature'),
                'time': conditions.get('time'),
                'atmosphere': conditions.get('atmosphere'),
            }
        })

    # Optional: attach a compact analytics snippet (Milestone 3)
    def _load_analytics_snippet(rt: str):
        try:
            if not isinstance(rt, str) or rt.strip().lower() != 'ullmann':
                return None
            base = os.path.join(os.path.dirname(__file__), 'data', 'analytics', 'Ullmann')
            latest = os.path.join(base, 'latest.json')
            if not os.path.exists(latest):
                return None
            with open(latest, 'r', encoding='utf-8') as f:
                summ = json.load(f)
            top = (summ.get('top') or {})
            co = (summ.get('cooccurrence') or {})
            def _simple(items, n=3):
                out = []
                for it in (items or [])[:n]:
                    out.append({
                        'name': it.get('name'),
                        'pct': it.get('pct'),
                        'count': it.get('count')
                    })
                return out
            def _best(pair_list):
                if not pair_list:
                    return None
                first = pair_list[0]
                return {
                    'a': first.get('a'),
                    'b': first.get('b'),
                    'pct': first.get('pct'),
                    'count': first.get('count')
                }
            return {
                'source': 'Ullmann',
                'top': {
                    'ligands': _simple(top.get('ligands')),
                    'solvents': _simple(top.get('solvents')),
                    'bases': _simple(top.get('bases')),
                },
                'cooccurrence': {
                    'best_ligand_solvent': _best(co.get('ligand_solvent')),
                    'best_base_solvent': _best(co.get('base_solvent')),
                }
            }
        except Exception:
            return None

    analytics_snippet = _load_analytics_snippet(detected_type)

    payload = {
        'meta': {
            'generated_at': time.strftime('%Y-%m-%dT%H:%M:%S'),
            'analysis_type': result.get('analysis_type', 'unknown'),
            'status': result.get('status'),
        },
        'input': {
            'reaction_smiles': result.get('reaction_smiles'),
            'selected_reaction_type': result.get('reaction_type') or result.get('selected_reaction_type'),
        },
        'detection': {
            'reaction_type': recs.get('reaction_type'),
        },
        'dataset': recs.get('dataset_info', {}),
        'recommendations': {
            'combined': combined,
            'ligands': ligands,
            'solvents': solvents,
            'alternatives': alternatives,
        },
        'top_conditions': top_conditions,
        'related_reactions': related,
    }
    if analytics_snippet:
        payload['analytics'] = analytics_snippet
    return payload
