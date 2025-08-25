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
        top_conditions: [ for top 3, a single 'chemicals' list with all required components (starting materials, metal precursor, ligand, base, solvent) and 'conditions' with time/temperature ],
    analytics: { source, top: { ligands[], solvents[], bases[] }, cooccurrence: { best_ligand_solvent, best_base_solvent }, numeric_stats: { temperature_c, time_h, yield_pct }, typical_catalyst_loading },
        related_reactions: [ { reaction_smiles, yield, catalyst, ligand, solvent, temperature, time, similarity, reaction_id, reference } ]
      }
    """

    recs = result.get('recommendations', {}) or {}

    # We no longer include the detailed 'recommendations' block in the export.

    related = related_reactions or recs.get('related_reactions', []) or []
    # Normalize any name|CAS tokens in related reactions (catalyst/ligand/solvent)
    def _name_only(tok: Optional[str]) -> Optional[str]:
        try:
            if tok is None:
                return None
            txt = str(tok)
            if '|' in txt:
                left, right = txt.split('|', 1)
                left = left.strip()
                right = right.strip()
                return left or right or None
            return txt.strip() or None
        except Exception:
            return str(tok).strip() if tok else None

    def _sanitize_field(val: Optional[str]) -> Optional[str]:
        if val is None:
            return None
        s = str(val).strip()
        if not s:
            return None
        # Try JSON array first
        if s.startswith('[') and s.endswith(']'):
            try:
                arr = json.loads(s)
                items = [_name_only(x) for x in arr if str(x).strip()]
                items = [x for x in items if x]
                return ", ".join(items) if items else None
            except Exception:
                pass
        # Fallback: split by commas
        parts = [p.strip() for p in s.split(',') if p.strip()]
        if len(parts) > 1:
            cleaned = [_name_only(p) for p in parts]
            cleaned = [x for x in cleaned if x]
            return ", ".join(cleaned) if cleaned else _name_only(s)
        return _name_only(s)

    if isinstance(related, list) and related:
        cleaned_related = []
        for r in related:
            try:
                if not isinstance(r, dict):
                    cleaned_related.append(r)
                    continue
                rr = dict(r)
                for key in ('catalyst', 'ligand', 'solvent'):
                    if key in rr:
                        rr[key] = _sanitize_field(rr.get(key))
                cleaned_related.append(rr)
            except Exception:
                cleaned_related.append(r)
        related = cleaned_related

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
                # CSV expected for cas_dictionary; unchanged
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

    # Optional: attach a compact analytics snippet (expanded)
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
            nums = (summ.get('numeric_stats') or {})
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
            snippet = {
                'source': 'Ullmann',
                'top': {
                    'ligands': _simple(top.get('ligands')),
                    'solvents': _simple(top.get('solvents')),
                    'bases': _simple(top.get('bases')),
                },
                'cooccurrence': {
                    'best_ligand_solvent': _best(co.get('ligand_solvent')),
                    'best_base_solvent': _best(co.get('base_solvent')),
                },
                'numeric_stats': {
                    'temperature_c': nums.get('temperature_c'),
                    'time_h': nums.get('time_h'),
                    'yield_pct': nums.get('yield_pct'),
                }
            }
            # Try to add a typical catalyst_loading string if available from combined conditions
            try:
                for c in (recs.get('combined_conditions') or []):
                    tc = c.get('typical_conditions') or {}
                    if tc.get('catalyst_loading'):
                        snippet['typical_catalyst_loading'] = tc.get('catalyst_loading')
                        break
            except Exception:
                pass
            return snippet
        except Exception:
            return None

    analytics_snippet = _load_analytics_snippet(detected_type)

    # Build dataset block and embed analytics (if available)
    dataset_block = recs.get('dataset_info', {}) or {}
    if not isinstance(dataset_block, dict):
        dataset_block = {}
    if analytics_snippet:
        dataset_block = dict(dataset_block)
        dataset_block['analytics'] = analytics_snippet

    payload = {
        'meta': {
            'generated_at': time.strftime('%Y-%m-%dT%H:%M:%S'),
            'analysis_type': result.get('analysis_type', 'unknown'),
            'status': result.get('status'),
            'providers': result.get('providers') or None,
        },
        'input': {
            'reaction_smiles': result.get('reaction_smiles'),
            'selected_reaction_type': result.get('reaction_type') or result.get('selected_reaction_type'),
        },
        'detection': {
            'reaction_type': recs.get('reaction_type'),
    },
    'dataset': dataset_block,
        'top_conditions': top_conditions,
        'related_reactions': related,
    }

    # Include cross-dataset general suggestions when available (GUI parity)
    general = result.get('general_recommendations') or recs.get('general_recommendations')
    if isinstance(general, dict) and general:
        payload['general'] = {
            'ligands': general.get('ligands', [])[:20],
            'solvents': general.get('solvents', [])[:20],
            'bases': general.get('bases', [])[:20],
            'top_hits': general.get('top_hits', [])[:50],
        }
    return payload
