#!/usr/bin/env python3
"""
Enhanced recommendation engine that integrates ligand and solvent recommendation systems
for comprehensive reaction condition predictions.
"""

import os
import sys
import json
from typing import Dict, List, Optional, Tuple
import re
from collections import defaultdict

# Ensure project root (containing 'reagents' package) is on sys.path
_HERE = os.path.dirname(__file__)
_ROOT = os.path.abspath(_HERE)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

try:
    from reagents.ligand import (
        recommend_ligands_for_reaction, 
        get_reaction_specific_ligands,
        create_ligand_dataframe
    )
    from reagents.solvent import (
        recommend_solvents_for_reaction,
        get_reaction_specific_solvents, 
        create_solvent_dataframe
    )
    ENHANCED_REAGENTS_AVAILABLE = True
except ImportError as e:
    print(f"Warning: Enhanced reagent systems not available: {e}")
    ENHANCED_REAGENTS_AVAILABLE = False

# Try to import existing recommendation engine
try:
    from recommendation_engine import RecommendationEngine as BaseRecommendationEngine
    BASE_ENGINE_AVAILABLE = True
except ImportError:
    BASE_ENGINE_AVAILABLE = False
    print("Base recommendation engine not available, using enhanced-only mode")

class EnhancedRecommendationEngine:
    """Enhanced recommendation engine with integrated ligand and solvent recommendations"""
    
    def __init__(self):
        self.base_engine = None
        # Analytics priors configuration (Milestone 2)
        self._analytics_cfg = {
            'enabled': True,
            'apply_to': {
                'solvents': True,
                'bases': True,
                'ligands': True,  # Milestone 3: enable ligand priors
            },
            # Base weight (fallback) and per-category overrides for frequency priors
            'w_freq': 0.30,
            'w_freq_solvents': 0.45,
            'w_freq_bases': 0.40,
            'w_freq_ligands': 0.35,
            # Soft penalty: down-weight items not supported by analytics (pct < min_support_pct)
            'soft_penalty': True,
            'penalty_factor': 0.85,
            'penalty_factor_solvents': 0.85,
            'penalty_factor_bases': 0.85,
            'penalty_factor_ligands': 0.88,
            'min_support_pct': 0.01,  # ignore items <1% support
        }
        if BASE_ENGINE_AVAILABLE:
            try:
                self.base_engine = BaseRecommendationEngine()
            except:
                pass
    
    def analyze_reaction_type(self, reaction_smiles: str, suggested_type: str = None) -> str:
        """Analyze and determine the reaction type from SMILES"""
        
        # If user specified a type, try to map it
        if suggested_type and suggested_type not in ("Auto-detect", "Auto detect reaction type"):
            mapped_type = self._map_reaction_type(suggested_type)
            if mapped_type:
                return mapped_type
        
        # Auto-detect based on SMILES pattern analysis
        if ">>" not in reaction_smiles:
            return "General Organic Reaction"
        
        try:
            reactants, products = reaction_smiles.split(">>", 1)
            
            # Look for common patterns
            if self._is_cross_coupling_pattern(reactants, products):
                return "Cross-Coupling"
            elif self._is_hydrogenation_pattern(reactants, products):
                return "Hydrogenation"
            elif self._is_carbonylation_pattern(reactants, products):
                return "Carbonylation"
            elif self._is_ch_activation_pattern(reactants, products):
                return "C-H_Activation"
            else:
                return "Cross-Coupling"  # Default to cross-coupling for organometallic reactions
                
        except Exception:
            return "General Organic Reaction"
    
    def _map_reaction_type(self, gui_type: str) -> Optional[str]:
        """Map GUI reaction types to our enhanced system types"""
        # Strip trailing metal tags like " (Pd)" or " (Cu)" from GUI label
        base_gui = gui_type
        try:
            if gui_type and gui_type.endswith(')') and ' (' in gui_type:
                base_gui = gui_type[:gui_type.rfind(' (')]
        except Exception:
            base_gui = gui_type

        mapping = {
            # Couplings
            "Suzuki-Miyaura Coupling": "Cross-Coupling",
            "C-C Coupling - Suzuki-Miyaura": "Cross-Coupling",
            "Buchwald-Hartwig Amination": "Cross-Coupling",
            "C-N Coupling - Buchwald-Hartwig": "Cross-Coupling",
            "Heck Coupling": "Cross-Coupling",
            "C-C Coupling - Heck": "Cross-Coupling",
            "Sonogashira Coupling": "Cross-Coupling",
            "C-C Coupling - Sonogashira": "Cross-Coupling",
            "Stille Coupling": "Cross-Coupling",
            "C-C Coupling - Stille": "Cross-Coupling",
            "Negishi Coupling": "Cross-Coupling",
            "C-C Coupling - Negishi": "Cross-Coupling",
            # Chan-Lam oxidative C-N coupling
            "Chan-Lam Coupling": "Cross-Coupling",
            "C-N Oxidative Coupling - Chan-Lam": "Cross-Coupling",
            "C-N Coupling - Chan-Lam": "Cross-Coupling",
            # Ullmann variants
            "Ullmann Ether Synthesis": "Ullmann",
            "Ullmann Reaction": "Ullmann",
            "C-N Coupling - Ullmann": "Ullmann",
            "C-O Coupling - Ullmann Ether": "Ullmann",
            "C-O Coupling - Ullmann": "Ullmann",
            # Other categories
            "Hydrogenation": "Hydrogenation",
            "Carbonylation": "Carbonylation",
            "Oxidation": "C-H_Activation",
            "C-H Activation": "C-H_Activation"
        }
        return mapping.get(base_gui)
    
    def _is_cross_coupling_pattern(self, reactants: str, products: str) -> bool:
        """Check if reaction pattern matches cross-coupling"""
        # Look for halogens, boronic acids, etc.
        halogen_pattern = r'[Br]|[Cl]|[I]'
        boron_pattern = r'B\(O\)'
        nitrogen_pattern = r'N[^a-z]'  # Nitrogen not part of aromatic system
        
        has_halogen = bool(re.search(halogen_pattern, reactants))
        has_boron = bool(re.search(boron_pattern, reactants))
        has_nitrogen = bool(re.search(nitrogen_pattern, reactants))
        
        return (has_halogen and has_boron) or (has_halogen and has_nitrogen)
    
    def _is_hydrogenation_pattern(self, reactants: str, products: str) -> bool:
        """Check if reaction pattern matches hydrogenation"""
        # Look for reduction of double bonds, carbonyls, etc.
        # Simple heuristic: check for decrease in unsaturation
        reactant_double_bonds = reactants.count('=') + reactants.count('#') * 2
        product_double_bonds = products.count('=') + products.count('#') * 2
        
        return product_double_bonds < reactant_double_bonds
    
    def _is_carbonylation_pattern(self, reactants: str, products: str) -> bool:
        """Check if reaction pattern matches carbonylation"""
        # Look for CO insertion patterns
        reactant_carbonyls = reactants.count('C=O') + reactants.count('C(=O)')
        product_carbonyls = products.count('C=O') + products.count('C(=O)')
        
        return product_carbonyls > reactant_carbonyls
    
    def _is_ch_activation_pattern(self, reactants: str, products: str) -> bool:
        """Check if reaction pattern matches C-H activation"""
        # Look for aromatic substitution patterns
        aromatic_reactants = reactants.count('c') + reactants.count('C')
        aromatic_products = products.count('c') + products.count('C')
        
        # Simple heuristic: more aromatic complexity in products
        return aromatic_products > aromatic_reactants * 1.2
    
    def get_recommendations(self, reaction_smiles: str, reaction_type: str = "Auto-detect") -> Dict:
        """Get comprehensive recommendations including ligands and solvents"""
        
        try:
            # Determine actual reaction type
            actual_reaction_type = self.analyze_reaction_type(reaction_smiles, reaction_type)
            
            result = {
                'analysis_type': 'enhanced',
                'reaction_type': actual_reaction_type,
                'detected_from': reaction_type,
                'status': 'success'
            }
            
            if not ENHANCED_REAGENTS_AVAILABLE:
                result.update({
                    'error': 'Enhanced reagent systems not available',
                    'message': 'Please ensure ligand.py and solvent.py are properly installed'
                })
                return result
            
            # Get enhanced recommendations
            enhanced_recs = self._get_enhanced_recommendations(reaction_smiles, actual_reaction_type)
            result.update(enhanced_recs)

            # If user left type as Auto-detect, also compute a general, cross-dataset
            # similarity-based recommendation set as supplemental guidance.
            try:
                if (reaction_type or '').lower().startswith('auto'):
                    gen = self._get_general_similarity_recommendations(reaction_smiles)
                    if gen:
                        result['general_recommendations'] = gen
            except Exception as _e:
                # Keep enhanced path intact even if general similarity fails
                result.setdefault('debug_warnings', []).append(f"general_similarity_failed: {str(_e)}")
            
            # Try to get base engine recommendations as supplementary info
            # Suppress Buchwald analysis for Ullmann selections/detections
            if self.base_engine and (actual_reaction_type or "").lower() != "ullmann":
                try:
                    base_recs = self.base_engine.get_recommendations(reaction_smiles, reaction_type)
                    if base_recs and base_recs.get('analysis_type') == 'buchwald_hartwig':
                        result['buchwald_hartwig_analysis'] = base_recs
                        result['analysis_type'] = 'comprehensive'  # Both enhanced + base
                except Exception as e:
                    print(f"Base engine error: {e}")
            
            return result
            
        except Exception as e:
            return {
                'analysis_type': 'error',
                'error': str(e),
                'status': 'failed'
            }
    
    def _get_enhanced_recommendations(self, reaction_smiles: str, reaction_type: str) -> Dict:
        """Get enhanced ligand and solvent recommendations"""
        
        recommendations = {
            'ligand_recommendations': [],
            'solvent_recommendations': [],
            'base_recommendations': [],
            'combined_conditions': [],
            'property_based_alternatives': {},
            'reaction_specific_notes': ""
        }
        
        try:
            # Normalize scoring mapping for internal parsing, but pass original type to allow Ullmann boosts
            scoring_type = 'Cross-Coupling' if (reaction_type or '').lower() == 'ullmann' else reaction_type

            # Analytics priors (latest.json) take precedence when present; fallback to CSV harvest
            priors = self._load_analytics_summary(reaction_type)
            evidence_ligands = None
            evidence_solvents = None
            evidence_bases = None
            if priors:
                try:
                    evidence_solvents = self._extract_priors(priors, 'solvents')
                except Exception:
                    evidence_solvents = None
                try:
                    evidence_bases = self._extract_priors(priors, 'bases')
                except Exception:
                    evidence_bases = None
                try:
                    evidence_ligands = self._extract_priors(priors, 'ligands')
                except Exception:
                    evidence_ligands = None
            else:
                # Evidence-aware context: mine dataset for ligands/solvents/bases used in similar reactions (if available)
                evidence_ligands = self._harvest_evidence_ligands(reaction_type)
                evidence_solvents = self._harvest_evidence_solvents(reaction_type)
                evidence_bases = self._harvest_evidence_bases(reaction_type)

            # Get top ligands for this reaction type, with evidence-aware boost
            ligands = recommend_ligands_for_reaction(
                reaction_type=reaction_type,
                top_n=5,
                min_compatibility=0.4,
                evidence_ligands=evidence_ligands or None
            )
            # Drop placeholder/dummy items
            if ligands:
                ligands = [L for L in ligands if str(L.get('ligand') or '').strip().lower() not in ('none', 'n/a', '-')]
            # Apply analytics priors to ligands if configured
            if priors and self._analytics_cfg['enabled'] and self._analytics_cfg['apply_to']['ligands']:
                try:
                    ligands = self._apply_freq_priors_ligands(ligands, priors)
                except Exception:
                    pass
            recommendations['ligand_recommendations'] = ligands
            
            # Get top solvents for this reaction type (pass evidence for gentle boost)
            solvents = recommend_solvents_for_reaction(
                reaction_type=reaction_type,
                top_n=5,
                min_compatibility=0.4,
                evidence_solvents=evidence_solvents or None
            )
            # Apply analytics priors to solvents if configured
            if priors and self._analytics_cfg['enabled'] and self._analytics_cfg['apply_to']['solvents']:
                try:
                    solvents = self._apply_freq_priors_solvents(solvents, priors)
                except Exception:
                    pass
            recommendations['solvent_recommendations'] = solvents
            
            # Try to get base recommendations if available
            try:
                from reagents.base import recommend_bases_for_reaction  # type: ignore
                bases = recommend_bases_for_reaction(
                    reaction_type=reaction_type,
                    top_n=5,
                    min_compatibility=0.4
                )
                # Apply analytics priors to bases if configured
                if priors and self._analytics_cfg['enabled'] and self._analytics_cfg['apply_to']['bases']:
                    try:
                        bases = self._apply_freq_priors_bases(bases, priors)
                    except Exception:
                        pass
                else:
                    # Legacy evidence-aware gentle boost to bases
                    try:
                        if evidence_bases:
                            # Normalize weights
                            max_w = max(float(v) for v in evidence_bases.values()) if evidence_bases else 0.0
                            if max_w > 0:
                                def _canon_base(nm: str) -> str:
                                    s = (nm or '').strip()
                                    if '(' in s and ')' in s:
                                        inner = s[s.rfind('(')+1:s.rfind(')')].strip()
                                        if inner:
                                            s = inner
                                    return s.lower()

                                ev_map = { _canon_base(k): float(v) for k, v in evidence_bases.items() }
                                boosted = []
                                for b in bases:
                                    name = b.get('base') or ''
                                    key = _canon_base(str(name))
                                    w = ev_map.get(key, 0.0)
                                    adjusted = b.get('compatibility_score', 0.0)
                                    if w > 0:
                                        norm = w / max_w
                                        boost = 0.05 + 0.10 * float(norm)
                                        adjusted = max(0.0, min(1.0, adjusted + boost))
                                    nb = dict(b)
                                    nb['compatibility_score'] = round(adjusted, 3)
                                    boosted.append(nb)
                                boosted.sort(key=lambda x: x.get('compatibility_score', 0.0), reverse=True)
                                bases = boosted
                    except Exception:
                        pass
                recommendations['base_recommendations'] = bases
            except Exception:
                bases = []

            # Create combined condition recommendations
            combined_conditions = self._create_combined_conditions(ligands, solvents, reaction_type)
            # Attach first base as a suggestion if present
            if bases:
                suggested = bases[0].get('base')
                for cc in combined_conditions:
                    cc['suggested_base'] = suggested
            recommendations['combined_conditions'] = combined_conditions
            
            # Get property-based alternatives
            property_alternatives = self._get_property_alternatives(reaction_type)
            recommendations['property_based_alternatives'] = property_alternatives
            
            # Add reaction-specific notes
            recommendations['reaction_specific_notes'] = self._get_reaction_notes(reaction_type)
            
            # Add dataset statistics
            try:
                ligand_df = create_ligand_dataframe()
                solvent_df = create_solvent_dataframe()
                recommendations['dataset_info'] = {
                    'ligands_available': len(ligand_df),
                    'solvents_available': len(solvent_df),
                    'reaction_types_supported': ['Cross-Coupling', 'Hydrogenation', 'Metathesis', 'C-H_Activation', 'Carbonylation'],
                    'analytics_loaded': bool(priors)
                }
            except:
                pass
            
        except Exception as e:
            recommendations['error'] = f"Enhanced recommendation error: {str(e)}"
        
        return recommendations

    # ===== General cross-dataset similarity (for undefined/auto type) =====
    def _get_general_similarity_recommendations(self, reaction_smiles: str) -> Optional[Dict]:
        """Search all datasets for similar reactions (by SMILES) and aggregate
        ligands, solvents, and bases from the top hits.

        Returns a dict with:
        - top_hits: list of {ReactionID, ReactionType, similarity, CondKey, Ligand[], Base[], Solvent[]}
        - ligand_recommendations/solvent_recommendations/base_recommendations: ranked with scores 0..1
        """
        try:
            # Lazy import RDKit; if not available, gracefully skip
            try:
                from rdkit import Chem
                from rdkit.Chem import AllChem
                from rdkit import DataStructs
            except Exception:
                return None

            def _mols_from_mixture(smistr: str) -> List:
                mols = []
                if not smistr:
                    return mols
                for part in str(smistr).split('.'):
                    part = part.strip()
                    if not part:
                        continue
                    try:
                        m = Chem.MolFromSmiles(part)
                        if m:
                            mols.append(m)
                    except Exception:
                        continue
                return mols

            def _fp_from_mixture(smistr: str):
                mols = _mols_from_mixture(smistr)
                if not mols:
                    return None
                fps = [AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048) for m in mols]
                # Combine by bitwise OR
                combo = fps[0]
                for fv in fps[1:]:
                    combo |= fv
                return combo

            # Build query fingerprints for reactants and products
            react_smi, prod_smi = '', ''
            if '>>' in reaction_smiles:
                react_smi, prod_smi = reaction_smiles.split('>>', 1)
            else:
                react_smi = reaction_smiles
                prod_smi = ''
            qfp_r = _fp_from_mixture(react_smi)
            qfp_p = _fp_from_mixture(prod_smi) if prod_smi else None

            if not qfp_r and not qfp_p:
                return None

            data_dir = os.path.join(_ROOT, 'data', 'reaction_dataset')
            if not os.path.isdir(data_dir):
                return None

            # Collect similarities across all dataset files (CSV/TSV)
            candidates: List[Dict] = []
            import csv
            for fname in os.listdir(data_dir):
                if not (fname.lower().endswith('.csv') or fname.lower().endswith('.tsv')):
                    continue
                path = os.path.join(data_dir, fname)
                try:
                    with open(path, 'r', encoding='utf-8') as f:
                        # Auto-select delimiter by extension
                        if fname.lower().endswith('.tsv'):
                            reader = csv.DictReader(f, delimiter='\t')
                        else:
                            reader = csv.DictReader(f)
                        for row in reader:
                            rs = row.get('ReactantSMILES') or ''
                            ps = row.get('ProductSMILES') or ''
                            rfp = _fp_from_mixture(rs)
                            pfp = _fp_from_mixture(ps) if ps else None
                            sim_r = 0.0
                            sim_p = 0.0
                            try:
                                if qfp_r is not None and rfp is not None:
                                    sim_r = DataStructs.TanimotoSimilarity(qfp_r, rfp)
                                if qfp_p is not None and pfp is not None:
                                    sim_p = DataStructs.TanimotoSimilarity(qfp_p, pfp)
                            except Exception:
                                pass
                            # Weighted combination; favor product when available
                            if qfp_p is not None and pfp is not None:
                                sim = 0.6 * sim_p + 0.4 * sim_r
                            else:
                                sim = sim_r
                            if sim <= 0:
                                continue
                            # Flexible field harvesting for richer hit details
                            yield_val = (
                                row.get('Yield') or row.get('Yield_%') or row.get('Yield(%)') or
                                row.get('Yield%') or row.get('Yield %') or row.get('Yield_pct') or ''
                            )
                            temp_val = (
                                row.get('Temperature') or row.get('Temp') or row.get('Temperature_C') or
                                row.get('TempC') or row.get('Temp_C') or ''
                            )
                            time_val = (
                                row.get('Time') or row.get('Hours') or row.get('Time_h') or row.get('Duration') or ''
                            )
                            catalyst_val = (
                                row.get('Catalyst') or row.get('Cat') or row.get('CopperSource') or
                                row.get('PdSource') or row.get('Metal') or ''
                            )
                            reference_val = (
                                row.get('Reference') or row.get('DOI') or row.get('URL') or
                                row.get('Source') or row.get('JournalRef') or ''
                            )
                            candidates.append({
                                'ReactionID': row.get('ReactionID') or '',
                                'ReactionType': row.get('ReactionType') or '',
                                'CondKey': row.get('CondKey') or '',
                                'ReactantSMILES': rs,
                                'ProductSMILES': ps,
                                'Ligand': row.get('Ligand') or '',
                                # Support both legacy and new columns
                                'Reagent': row.get('Reagent') or row.get('ReagentRaw') or '',
                                'ReagentRole': row.get('ReagentRole') or '',
                                'RGTName': row.get('RGTName') or '',
                                'SOLName': row.get('SOLName') or '',
                                'Solvent': row.get('Solvent') or '',
                                'CoreDetail': row.get('CoreDetail') or '',
                                'CoreGeneric': row.get('CoreGeneric') or '',
                                'YieldPct': yield_val,
                                'Temperature': temp_val,
                                'Time': time_val,
                                'CatalystLike': catalyst_val,
                                'Reference': reference_val,
                                'DatasetFile': fname,
                                'similarity': float(sim)
                            })
                except Exception:
                    continue

            if not candidates:
                return None

            # Keep top K hits
            candidates.sort(key=lambda x: x['similarity'], reverse=True)
            top_hits = candidates[:50]

            # Aggregate ligands, solvents, bases weighted by similarity
            lig_counts: Dict[str, float] = defaultdict(float)
            solv_counts: Dict[str, float] = defaultdict(float)
            base_counts: Dict[str, float] = defaultdict(float)

            def _parse_listlike(val: str) -> List[str]:
                if not val:
                    return []
                s = str(val).strip()
                # Try JSON list first
                if s.startswith('[') and s.endswith(']'):
                    try:
                        arr = json.loads(s)
                        return [str(x).strip() for x in arr if str(x).strip()]
                    except Exception:
                        pass
                # Fallback: strip quotes/brackets and split by comma
                s2 = s.strip('[]').replace('"', '').replace("'", '')
                parts = [p.strip() for p in s2.split(',') if p.strip()]
                return parts

            # New datasets may encode tokens as "name|CAS"; prefer name, fallback to CAS if name missing
            def _name_only(tok: str) -> str:
                try:
                    if tok is None:
                        return ''
                    txt = str(tok)
                    if '|' in txt:
                        left, right = txt.split('|', 1)
                        left = left.strip()
                        right = right.strip()
                        return left or right or ''
                    return txt.strip()
                except Exception:
                    return str(tok).strip() if tok is not None else ''

            # Base alias normalization
            base_alias = {
                'potassium carbonate (k2co3)': 'K2CO3',
                'cesium carbonate (cs2co3)': 'Cs2CO3',
                'tripotassium phosphate (k3po4)': 'K3PO4',
                'potassium tert-butoxide (kotbu)': 'KOtBu',
                'sodium tert-butoxide (naotbu)': 'NaOtBu',
                'sodium carbonate (na2co3)': 'Na2CO3',
                'potassium hydroxide (koh)': 'KOH',
                'triethylamine': 'Et3N',
            }

            def _canon_base(nm: str) -> str:
                s = (nm or '').strip()
                low = s.lower()
                # Pull formula inside parentheses when present
                if '(' in s and ')' in s:
                    inner = s[s.rfind('(')+1:s.rfind(')')].strip()
                    if inner:
                        low = inner.lower()
                low = low.replace(' ', '')
                return base_alias.get(low, s)

            # Solvent abbreviation lookup via solvent dataframe if available
            solv_abbrev: Dict[str, str] = {}
            try:
                sdf = create_solvent_dataframe()
                for _, rowx in sdf.iterrows():
                    nm = str(rowx.get('name') or '').strip()
                    ab = str(rowx.get('abbreviation') or '').strip()
                    if nm:
                        solv_abbrev[nm] = ab
                    if ab and ab not in solv_abbrev:
                        solv_abbrev[ab] = ab
            except Exception:
                pass

            for hit in top_hits:
                w = float(hit.get('similarity') or 0.0)
                if w <= 0:
                    continue
                # Ligands
                for src in ('Ligand',):
                    for item in _parse_listlike(hit.get(src) or ''):
                        name = _name_only(item)
                        if name and name.lower() != 'none':
                            lig_counts[name] += w
                # Solvents (handle both columns)
                sol_items = _parse_listlike(hit.get('Solvent') or '')
                sol_items += _parse_listlike(hit.get('SOLName') or '')
                for sitem in sol_items:
                    name = _name_only(sitem)
                    if name and name.lower() != 'none':
                        solv_counts[name] += w
                # Bases (from reagent columns)
                # New unified reagent columns
                rraw = _parse_listlike(hit.get('Reagent') or '')
                rgtn = _parse_listlike(hit.get('RGTName') or '')
                for bitem in rraw + rgtn:
                    cb = _canon_base(_name_only(bitem))
                    if cb and cb.lower() not in ('none', 'unk'):
                        base_counts[cb] += w

            # Normalize to 0..1 by max weight
            def _rank_map(d: Dict[str, float]) -> List[Tuple[str, float]]:
                if not d:
                    return []
                mx = max(d.values()) if d else 1.0
                return [(k, (v / mx) if mx > 0 else 0.0) for k, v in sorted(d.items(), key=lambda kv: kv[1], reverse=True)]

            lig_rank = _rank_map(lig_counts)
            solv_rank = _rank_map(solv_counts)
            base_rank = _rank_map(base_counts)

            # Adapt to engine's output shapes
            # Utility: detect CAS-only tokens (e.g., 108-88-3)
            cas_re = re.compile(r"^\d{2,7}-\d{2}-\d$")
            def _filter_cas(rank_list: List[Tuple[str, float]]) -> List[Tuple[str, float]]:
                non_cas = [(n, s) for n, s in rank_list if not cas_re.match(str(n or '').strip())]
                return non_cas if non_cas else rank_list

            lig_rank = _filter_cas(lig_rank)
            solv_rank = _filter_cas(solv_rank)
            base_rank = _filter_cas(base_rank)

            lig_out = [{
                'ligand': name,
                'compatibility_score': round(score, 3)
            } for name, score in lig_rank[:5]]

            solv_out = [{
                'solvent': name,
                'abbreviation': (solv_abbrev.get(name) or name),
                'compatibility_score': round(score, 3)
            } for name, score in solv_rank[:5]]

            base_out = [{
                'base': name,
                'compatibility_score': round(score, 3)
            } for name, score in base_rank[:5]]

            # Provide a richer view of top hits (top 15)
            base_tokens = ['k2co3', 'cs2co3', 'k3po4', 'kotbu', 'naotbu', 'na2co3', 'koh', 'tbuok', 'ko-tbu', 'triethylamine', 'et3n']
            cas_only_re = re.compile(r'^\d{2,7}-\d{2}-\d$')
            def _dedup(seq: List[str]) -> List[str]:
                seen = set()
                out: List[str] = []
                for x in seq:
                    k = (x or '').strip().lower()
                    if not k or k in seen:
                        continue
                    seen.add(k)
                    out.append(x)
                return out

            def _extract_ligands_from_hit(h: Dict) -> List[str]:
                ligs = []
                for item in _parse_listlike(h.get('Ligand') or ''):
                    nm = _name_only(item)
                    if nm and nm.lower() not in ('none',) and not cas_only_re.match(nm):
                        ligs.append(nm)
                return _dedup(ligs)

            def _extract_solvents_from_hit(h: Dict) -> List[str]:
                sol_items = _parse_listlike(h.get('Solvent') or '') + _parse_listlike(h.get('SOLName') or '')
                sols = []
                for sitem in sol_items:
                    nm = _name_only(sitem)
                    if nm and nm.lower() not in ('none',) and not cas_only_re.match(nm):
                        sols.append(nm)
                return _dedup(sols)

            def _extract_bases_from_hit(h: Dict) -> List[str]:
                bases = []
                for bitem in _parse_listlike(h.get('Reagent') or '') + _parse_listlike(h.get('RGTName') or ''):
                    nm = _name_only(bitem)
                    low = nm.lower().replace(' ', '')
                    if cas_only_re.match(nm):
                        continue
                    if any(tok in low for tok in base_tokens):
                        bases.append(_canon_base(nm))
                return _dedup(bases)

            def _as_float(val: str) -> Optional[float]:
                try:
                    s = str(val).strip()
                    if not s:
                        return None
                    # strip percent sign or other non-numeric
                    s = s.replace('%', '').replace('\u200b', '').strip()
                    return float(s)
                except Exception:
                    return None

            hits_view = []
            for h in top_hits[:15]:
                y = _as_float(h.get('YieldPct') or '')
                if y is not None and y > 100:
                    y = min(y, 100.0)
                hits_view.append({
                    'reaction_id': h.get('ReactionID') or None,
                    'reaction_type': h.get('ReactionType') or None,
                    'dataset_file': h.get('DatasetFile') or None,
                    'cond_key': h.get('CondKey') or None,
                    'similarity': round(float(h.get('similarity') or 0.0), 3),
                    'reactant_smiles': h.get('ReactantSMILES') or '',
                    'product_smiles': h.get('ProductSMILES') or '',
                    'reaction_smiles': (f"{h.get('ReactantSMILES')}>>{h.get('ProductSMILES')}"
                                        if (h.get('ReactantSMILES') and h.get('ProductSMILES')) else None),
                    'ligands': _extract_ligands_from_hit(h),
                    'solvents': _extract_solvents_from_hit(h),
                    'bases': _extract_bases_from_hit(h),
                    'temperature': (h.get('Temperature') or None),
                    'time': (h.get('Time') or None),
                    'yield_pct': y,
                    'catalyst': (h.get('CatalystLike') or None),
                    'reference': (h.get('Reference') or None),
                    'core_detail': h.get('CoreDetail') or None,
                    'core_generic': h.get('CoreGeneric') or None,
                })

            # Create combined suggestions for convenience
            combined_conditions = self._create_combined_conditions(lig_out, solv_out, self.analyze_reaction_type(reaction_smiles, None)) if lig_out and solv_out else []
            if base_out and combined_conditions:
                suggested = base_out[0].get('base')
                for cc in combined_conditions:
                    cc['suggested_base'] = suggested

            return {
                'top_hits': hits_view,
                'ligand_recommendations': lig_out,
                'solvent_recommendations': solv_out,
                'base_recommendations': base_out,
                'combined_conditions': combined_conditions
            }
        except Exception as e:
            return {'error': f'general_similarity_error: {e}'}

    # ===== Milestone 2: analytics loading and priors application =====
    def _load_analytics_summary(self, reaction_type: str) -> Optional[dict]:
        """Load data/analytics/<reaction_type>/latest.json if present.

        For Milestone 2, we support Ullmann only.
        """
        try:
            rt = (reaction_type or '').strip()
            if rt.lower() != 'ullmann':
                return None
            base = os.path.join(_ROOT, 'data', 'analytics', 'Ullmann')
            latest = os.path.join(base, 'latest.json')
            if os.path.exists(latest):
                with open(latest, 'r', encoding='utf-8') as f:
                    data = json.load(f)
                    return data
        except Exception:
            return None
        return None

    def _extract_priors(self, summary: dict, kind: str) -> Optional[dict]:
        """Extract a mapping name -> pct for a top list kind (e.g., 'solvents', 'bases').
        Returns dict of canonical_name -> pct (0..1).
        """
        try:
            top = ((summary or {}).get('top') or {}).get(kind) or []
            if not top:
                return None
            pri = {}
            for item in top:
                nm = str(item.get('name') or '').strip()
                pct = float(item.get('pct') or 0.0)
                if not nm or pct <= 0:
                    continue
                pri[nm] = pct
            return pri or None
        except Exception:
            return None

    def _apply_freq_priors_solvents(self, solvents: List[Dict], summary: dict) -> List[Dict]:
        """Apply frequency-based priors to solvent compatibility scores.

        final = base * (1 + w * sqrt(pct)) with cap to 1.0
        """
        try:
            import math
            pri = self._extract_priors(summary, 'solvents') or {}
            if not pri:
                return solvents
            # Optional normalization helpers
            def _canon(name: str) -> str:
                s = (name or '').strip()
                # try to harmonize common cases
                s_low = s.lower().replace(' ', '')
                # special-case DMSO often appears as 'dms o' canonical in analytics
                if s_low in ('dmso', 'dms o', 'dimethylsulfoxide'):
                    return 'dms o'
                return s_low

            # Build lookup
            pri_map = { _canon(k): float(v) for k, v in pri.items() }
            w = float(self._analytics_cfg.get('w_freq_solvents', self._analytics_cfg.get('w_freq', 0.30)))
            min_pct = float(self._analytics_cfg.get('min_support_pct', 0.01) or 0.01)
            boosted: List[Dict] = []
            for s in solvents:
                base_score = float(s.get('compatibility_score', 0.0) or 0.0)
                key = _canon(str(s.get('solvent') or ''))
                pct = pri_map.get(key, 0.0)
                adj = base_score
                if base_score > 0:
                    if pct >= min_pct:
                        adj = min(1.0, base_score * (1.0 + w * math.sqrt(pct)))
                    elif self._analytics_cfg.get('soft_penalty', True):
                        pen = float(self._analytics_cfg.get('penalty_factor_solvents', self._analytics_cfg.get('penalty_factor', 0.85)))
                        adj = max(0.0, base_score * pen)
                ns = dict(s)
                ns['compatibility_score'] = round(adj, 3)
                boosted.append(ns)
            boosted.sort(key=lambda x: x.get('compatibility_score', 0.0), reverse=True)
            return boosted
        except Exception:
            return solvents

    def _apply_freq_priors_bases(self, bases: List[Dict], summary: dict) -> List[Dict]:
        """Apply frequency-based priors to base compatibility scores using analytics.
        """
        try:
            import math
            pri = self._extract_priors(summary, 'bases') or {}
            if not pri:
                return bases

            def _canon_base(nm: str) -> str:
                s = (nm or '').strip()
                # favor formula inside parentheses if present
                if '(' in s and ')' in s:
                    inner = s[s.rfind('(')+1:s.rfind(')')].strip()
                    if inner:
                        s = inner
                return s.lower().replace(' ', '')

            pri_map = { _canon_base(k): float(v) for k, v in pri.items() }
            w = float(self._analytics_cfg.get('w_freq_bases', self._analytics_cfg.get('w_freq', 0.30)))
            min_pct = float(self._analytics_cfg.get('min_support_pct', 0.01) or 0.01)
            out: List[Dict] = []
            for b in bases:
                base_score = float(b.get('compatibility_score', 0.0) or 0.0)
                key = _canon_base(str(b.get('base') or ''))
                pct = pri_map.get(key, 0.0)
                adj = base_score
                if base_score > 0:
                    if pct >= min_pct:
                        adj = min(1.0, base_score * (1.0 + w * math.sqrt(pct)))
                    elif self._analytics_cfg.get('soft_penalty', True):
                        pen = float(self._analytics_cfg.get('penalty_factor_bases', self._analytics_cfg.get('penalty_factor', 0.85)))
                        adj = max(0.0, base_score * pen)
                nb = dict(b)
                nb['compatibility_score'] = round(adj, 3)
                out.append(nb)
            out.sort(key=lambda x: x.get('compatibility_score', 0.0), reverse=True)
            return out
        except Exception:
            return bases

    def _apply_freq_priors_ligands(self, ligands: List[Dict], summary: dict) -> List[Dict]:
        """Apply frequency-based priors to ligand compatibility scores using analytics.

        Similar to solvents/bases. Uses ligand name key directly.
        """
        try:
            import math
            pri = self._extract_priors(summary, 'ligands') or {}
            if not pri:
                return ligands

            def _canon(nm: str) -> str:
                return (nm or '').strip().casefold()

            pri_map = { _canon(k): float(v) for k, v in pri.items() }
            w = float(self._analytics_cfg.get('w_freq_ligands', self._analytics_cfg.get('w_freq', 0.30)))
            min_pct = float(self._analytics_cfg.get('min_support_pct', 0.01) or 0.01)
            out: List[Dict] = []
            for L in ligands:
                base_score = float(L.get('compatibility_score', 0.0) or 0.0)
                key = _canon(str(L.get('ligand') or ''))
                pct = pri_map.get(key, 0.0)
                adj = base_score
                if base_score > 0:
                    if pct >= min_pct:
                        adj = min(1.0, base_score * (1.0 + w * math.sqrt(pct)))
                    elif self._analytics_cfg.get('soft_penalty', True):
                        pen = float(self._analytics_cfg.get('penalty_factor_ligands', self._analytics_cfg.get('penalty_factor', 0.85)))
                        adj = max(0.0, base_score * pen)
                nL = dict(L)
                nL['compatibility_score'] = round(adj, 3)
                out.append(nL)
            out.sort(key=lambda x: x.get('compatibility_score', 0.0), reverse=True)
            return out
        except Exception:
            return ligands

    def _harvest_evidence_ligands(self, reaction_type: str) -> dict:
        """Collect a small frequency map of ligands from built-in datasets for this reaction type.

        Keeps it lightweight and offline: scans CSVs in data/reaction_dataset for matching ReactionType.
        """
        evidence: dict[str, float] = {}
        try:
            data_dir = os.path.join(_ROOT, 'data', 'reaction_dataset')
            if not os.path.isdir(data_dir):
                return evidence
            # Consider a few known files; ignore huge processing
            # Helper for name|CAS collapsing
            def _name_only(tok: str) -> str:
                try:
                    if tok is None:
                        return ''
                    txt = str(tok)
                    if '|' in txt:
                        left, right = txt.split('|', 1)
                        left = left.strip()
                        right = right.strip()
                        return left or right or ''
                    return txt.strip()
                except Exception:
                    return str(tok).strip() if tok is not None else ''

            for fname in os.listdir(data_dir):
                if not (fname.lower().endswith('.csv') or fname.lower().endswith('.tsv')):
                    continue
                path = os.path.join(data_dir, fname)
                try:
                    import csv
                    with open(path, 'r', encoding='utf-8') as f:
                        reader = csv.DictReader(f, delimiter='\t') if fname.lower().endswith('.tsv') else csv.DictReader(f)
                        for row in reader:
                            rtype = (row.get('ReactionType') or '').strip()
                            if not rtype:
                                continue
                            # Simple match: same type token (case-insensitive)
                            if rtype.lower() != (reaction_type or '').lower():
                                continue
                            lig_raw = row.get('Ligand') or ''
                            if not lig_raw:
                                continue
                            # Ligand field often like ["L-Proline"] or a list string
                            # Strip brackets/quotes and split by comma when safe
                            items = []
                            s = str(lig_raw).strip()
                            # crude parse: remove [] and quotes
                            s = s.strip('[]')
                            s = s.replace('"', '').replace("'", '')
                            # split by comma only if present
                            parts = [p.strip() for p in s.split(',') if p.strip()]
                            items = parts if parts else ([s] if s else [])
                            for it in items:
                                if not it:
                                    continue
                                name = _name_only(it)
                                evidence[name] = evidence.get(name, 0) + 1
                except Exception:
                    # ignore a bad file
                    continue
            # retain top few
            if evidence:
                top = sorted(evidence.items(), key=lambda kv: kv[1], reverse=True)[:10]
                return {k: float(v) for k, v in top}
        except Exception:
            pass
        return evidence

    def _harvest_evidence_solvents(self, reaction_type: str) -> dict:
        """Collect frequency map of solvents from built-in datasets for this reaction type."""
        evidence: dict[str, float] = {}
        try:
            data_dir = os.path.join(_ROOT, 'data', 'reaction_dataset')
            if not os.path.isdir(data_dir):
                return evidence
            def _name_only(tok: str) -> str:
                try:
                    if tok is None:
                        return ''
                    txt = str(tok)
                    if '|' in txt:
                        left, right = txt.split('|', 1)
                        left = left.strip()
                        right = right.strip()
                        return left or right or ''
                    return txt.strip()
                except Exception:
                    return str(tok).strip() if tok is not None else ''
            for fname in os.listdir(data_dir):
                if not (fname.lower().endswith('.csv') or fname.lower().endswith('.tsv')):
                    continue
                path = os.path.join(data_dir, fname)
                try:
                    import csv
                    with open(path, 'r', encoding='utf-8') as f:
                        reader = csv.DictReader(f, delimiter='\t') if fname.lower().endswith('.tsv') else csv.DictReader(f)
                        for row in reader:
                            rtype = (row.get('ReactionType') or '').strip()
                            if not rtype or rtype.lower() != (reaction_type or '').lower():
                                continue
                            raw = row.get('Solvent') or row.get('SOLName') or ''
                            if not raw:
                                continue
                            s = str(raw).strip().strip('[]').replace('"', '').replace("'", '')
                            parts = [p.strip() for p in s.split(',') if p.strip()]
                            items = parts if parts else ([s] if s else [])
                            for it in items:
                                if not it:
                                    continue
                                name = _name_only(it)
                                evidence[name] = evidence.get(name, 0) + 1
                except Exception:
                    continue
            if evidence:
                top = sorted(evidence.items(), key=lambda kv: kv[1], reverse=True)[:10]
                return {k: float(v) for k, v in top}
        except Exception:
            pass
        return evidence

    def _harvest_evidence_bases(self, reaction_type: str) -> dict:
        """Collect frequency map of bases from built-in datasets for this reaction type.

        The datasets may store bases under different columns: ReagentRaw/ReagentRole or RGTName, etc.
        We'll look for common base tokens and increment counts by base name.
        """
        evidence: dict[str, float] = {}
        try:
            data_dir = os.path.join(_ROOT, 'data', 'reaction_dataset')
            if not os.path.isdir(data_dir):
                return evidence
            base_tokens = ['k2co3', 'cs2co3', 'k3po4', 'kotbu', 'naotbu', 'na2co3', 'koh', 'tbuok', 'ko-tbu', 'triethylamine', 'et3n']
            # Helper for name|CAS collapsing
            def _name_only(tok: str) -> str:
                try:
                    if tok is None:
                        return ''
                    txt = str(tok)
                    if '|' in txt:
                        left, right = txt.split('|', 1)
                        left = left.strip()
                        right = right.strip()
                        return left or right or ''
                    return txt.strip()
                except Exception:
                    return str(tok).strip() if tok is not None else ''
            def _maybe_add(text: str):
                if not text:
                    return
                s = str(text).strip().strip('[]').replace('"', '').replace("'", '')
                parts = [p.strip() for p in s.split(',') if p.strip()]
                items = parts if parts else ([s] if s else [])
                for it in items:
                    low = it.lower()
                    for tok in base_tokens:
                        if tok in low:
                            name = _name_only(it)
                            evidence[name] = evidence.get(name, 0) + 1
                            break

            for fname in os.listdir(data_dir):
                if not (fname.lower().endswith('.csv') or fname.lower().endswith('.tsv')):
                    continue
                path = os.path.join(data_dir, fname)
                try:
                    import csv
                    with open(path, 'r', encoding='utf-8') as f:
                        reader = csv.DictReader(f, delimiter='\t') if fname.lower().endswith('.tsv') else csv.DictReader(f)
                        for row in reader:
                            rtype = (row.get('ReactionType') or '').strip()
                            if not rtype or rtype.lower() != (reaction_type or '').lower():
                                continue
                            # check common columns
                            # New column name Reagent (with roles in ReagentRole)
                            _maybe_add(row.get('Reagent') or row.get('ReagentRaw') or '')
                            _maybe_add(row.get('RGTName') or '')
                            _maybe_add(row.get('Base') or '')
                except Exception:
                    continue
            if evidence:
                top = sorted(evidence.items(), key=lambda kv: kv[1], reverse=True)[:10]
                return {k: float(v) for k, v in top}
        except Exception:
            pass
        return evidence
    
    def _create_combined_conditions(self, ligands: List[Dict], solvents: List[Dict], reaction_type: str) -> List[Dict]:
        """Create optimized ligand-solvent combinations"""
        
        combined = []
        
        # Take top 3 ligands and top 3 solvents
        top_ligands = ligands[:3]
        top_solvents = solvents[:3]
        
        for i, ligand in enumerate(top_ligands):
            for j, solvent in enumerate(top_solvents):
                # Calculate combined score
                combined_score = (ligand['compatibility_score'] + solvent['compatibility_score']) / 2
                
                # Add synergy bonus for known good combinations
                synergy_bonus = self._calculate_synergy_bonus(
                    ligand['ligand'], solvent['solvent'], reaction_type
                )
                
                final_score = combined_score + synergy_bonus
                
                combined.append({
                    'rank': len(combined) + 1,
                    'ligand': ligand['ligand'],
                    'ligand_compatibility': ligand['compatibility_score'],
                    'solvent': solvent['solvent'],
                    'solvent_abbreviation': solvent['abbreviation'],
                    'solvent_compatibility': solvent['compatibility_score'],
                    'combined_score': round(final_score, 3),
                    'synergy_bonus': round(synergy_bonus, 3),
                    'recommendation_confidence': 'High' if final_score > 0.8 else 'Medium' if final_score > 0.6 else 'Low',
                    'typical_conditions': self._get_typical_conditions(ligand['ligand'], solvent['solvent'], reaction_type)
                })
        
        # Sort by combined score and return top 5
        combined.sort(key=lambda x: x['combined_score'], reverse=True)
        return combined[:5]
    
    def _calculate_synergy_bonus(self, ligand: str, solvent: str, reaction_type: str) -> float:
        """Calculate synergy bonus for specific ligand-solvent combinations"""
        
        # Known synergistic combinations
        synergies = {
            'Cross-Coupling': {
                ('SPhos', 'DMF'): 0.1,
                ('XPhos', 'THF'): 0.1,
                ('RuPhos', 'DMF'): 0.08,
                ('BINAP', 'Toluene'): 0.05,
                ('PPh3', 'THF'): 0.05
            },
            'Ullmann': {
                ('1,10-Phenanthroline', 'DMSO'): 0.10,
                ("2,2'-Bipyridine", 'DMSO'): 0.10,
                ('L-Proline', 'DMSO'): 0.08,
                ('Ethylenediamine', 'DMF'): 0.08,
                ('DMEDA', 'Toluene'): 0.06,
            },
            'Hydrogenation': {
                ('BINAP', 'Ethanol'): 0.15,
                ('Tol-BINAP', 'Methanol'): 0.15,
                ('PPh3', 'Ethanol'): 0.08,
                ('DPPF', 'Ethanol'): 0.08
            },
            'Metathesis': {
                ('IPr', 'Dichloromethane'): 0.12,
                ('IMes', 'Dichloromethane'): 0.12,
                ('SIPr', 'Toluene'): 0.1
            }
        }
        
        reaction_synergies = synergies.get(reaction_type, {})
        return reaction_synergies.get((ligand, solvent), 0.0)
    
    def _get_typical_conditions(self, ligand: str, solvent: str, reaction_type: str) -> Dict:
        """Get typical reaction conditions for ligand-solvent combination"""
        
        # Default conditions by reaction type
        conditions = {
            'Cross-Coupling': {
                'temperature': '80-120°C',
                'time': '4-24 hours',
                'atmosphere': 'Inert (N₂ or Ar)',
                'base': 'K₂CO₃ or Cs₂CO₃',
                'catalyst_loading': '1-5 mol%'
            },
            'Ullmann': {
                'temperature': '80-140°C',
                'time': '6-24 hours',
                'atmosphere': 'Inert (N₂ or Ar)',
                'base': 'K₂CO₃, Cs₂CO₃, K₃PO₄ or KOtBu',
                'catalyst_loading': '5-20 mol% Cu',
                'additives': 'Ligands: phen, bipy, L-proline, diamines'
            },
            'Hydrogenation': {
                'temperature': '20-80°C',
                'time': '2-16 hours',
                'atmosphere': 'H₂ (1-50 atm)',
                'catalyst_loading': '0.1-2 mol%',
                'additives': 'May require acid'
            },
            'Metathesis': {
                'temperature': '20-60°C',
                'time': '1-8 hours',
                'atmosphere': 'Inert (N₂ or Ar)',
                'catalyst_loading': '1-5 mol%',
                'additives': 'Avoid moisture'
            },
            'C-H_Activation': {
                'temperature': '100-160°C',
                'time': '6-48 hours',
                'atmosphere': 'Inert or air',
                'catalyst_loading': '5-10 mol%',
                'additives': 'May require oxidant'
            },
            'Carbonylation': {
                'temperature': '60-140°C',
                'time': '4-24 hours',
                'atmosphere': 'CO (1-20 atm)',
                'catalyst_loading': '1-5 mol%',
                'base': 'Organic base (Et₃N)'
            }
        }
        
        return conditions.get(reaction_type, {
            'temperature': '20-100°C',
            'time': '1-24 hours',
            'atmosphere': 'Inert',
            'catalyst_loading': '1-5 mol%'
        })
    
    def _get_property_alternatives(self, reaction_type: str) -> Dict:
        """Get property-based alternative recommendations"""
        
        alternatives = {}
        
        try:
            # Get budget-friendly options
            budget_ligands = get_reaction_specific_ligands(
                reaction_type=reaction_type,
                property_preferences={'price_category_max': 3}
            )
            if budget_ligands:
                alternatives['budget_friendly_ligands'] = budget_ligands[:3]
            
            # Get low BP solvents for easy removal
            low_bp_solvents = get_reaction_specific_solvents(
                reaction_type=reaction_type,
                property_preferences={'bp_max': 100}
            )
            if low_bp_solvents:
                alternatives['low_boiling_solvents'] = low_bp_solvents[:3]
            
            # Get green/sustainable options
            green_solvents = get_reaction_specific_solvents(
                reaction_type=reaction_type,
                property_preferences={'polarity_min': 3, 'bp_max': 150}
            )
            if green_solvents:
                alternatives['green_solvents'] = green_solvents[:3]
            
        except Exception as e:
            alternatives['error'] = f"Could not generate alternatives: {e}"
        
        return alternatives
    
    def _get_reaction_notes(self, reaction_type: str) -> str:
        """Get reaction-specific guidance notes"""
        
        notes = {
            'Cross-Coupling': """
💡 Cross-Coupling Optimization Tips:
• Use bulky phosphines (XPhos, SPhos) for challenging substrates
• Polar aprotic solvents (DMF, NMP) often give best results
• Consider base choice: K₂CO₃ for most substrates, Cs₂CO₃ for difficult cases
• Temperature typically 80-120°C depending on substrate reactivity
• Degassing is critical - use Schlenk techniques or glovebox
            """,
            'Ullmann': """
💡 Ullmann Coupling Optimization Tips:
• Copper sources: CuI, CuBr, Cu(OAc)₂, Cu₂O; often with simple ligands
• Ligands: diamines (e.g., ethylenediamine), amino acids (e.g., L-proline), phenanthroline
• Bases: K₂CO₃, Cs₂CO₃, K₃PO₄, KOtBu; water sometimes beneficial
• Solvents: DMSO, DMF, toluene, dioxane; 80–140°C typical
• For C–O/C–N: substrate electronics impact rates; consider stronger base for aryl chlorides
            """,
            'Hydrogenation': """
💡 Hydrogenation Optimization Tips:
• Bidentate ligands (BINAP, DuPhos) excellent for asymmetric reductions
• Protic solvents (alcohols) often enhance reactivity
• Start with low pressure (1-5 atm H₂) and increase if needed
• Temperature usually mild (20-80°C) to avoid over-reduction
• Check for catalyst poisoning from sulfur/nitrogen compounds
            """,
            'Metathesis': """
💡 Metathesis Optimization Tips:
• NHC ligands (IPr, IMes) provide high activity and stability
• Non-coordinating solvents (DCM, toluene) are preferred
• Strict exclusion of moisture and oxygen is essential
• Low catalyst loadings (1-5 mol%) usually sufficient
• Consider ring-closing vs cross-metathesis selectivity
            """,
            'C-H_Activation': """
💡 C-H Activation Optimization Tips:
• High temperatures (100-160°C) often required
• Polar solvents (DMSO, DMF) can facilitate C-H cleavage
• Consider directing groups for regioselectivity
• Oxidants may be required for catalytic turnover
• Screen different bases for optimal reactivity
            """,
            'Carbonylation': """
💡 Carbonylation Optimization Tips:
• CO pressure critical for good conversion (1-20 atm)
• Polar solvents (DMF, NMP) enhance CO solubility
• Phosphine ligands (PPh3, DPPF) commonly effective
• Base helps remove HX byproducts
• Monitor for catalyst degradation at high CO pressure
            """
        }
        
        return notes.get(reaction_type, "General organometallic reaction guidelines apply.")
    
    def get_available_recommenders(self) -> List[str]:
        """Get list of available recommendation systems"""
        recommenders = []
        
        if ENHANCED_REAGENTS_AVAILABLE:
            recommenders.extend(['Enhanced Ligand System', 'Enhanced Solvent System'])
        
        if self.base_engine:
            recommenders.extend(self.base_engine.get_available_recommenders())
        
        return recommenders

def create_recommendation_engine():
    """Factory function to create the enhanced recommendation engine"""
    return EnhancedRecommendationEngine()

# For backward compatibility
RecommendationEngine = EnhancedRecommendationEngine

if __name__ == "__main__":
    # Test the enhanced recommendation engine
    engine = create_recommendation_engine()
    
    # Test cross-coupling recommendation
    test_smiles = "c1ccc(Br)cc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1"
    recommendations = engine.get_recommendations(test_smiles, "Suzuki-Miyaura Coupling")
    
    print("Enhanced Recommendation Engine Test")
    print("=" * 50)
    print(f"Reaction: {test_smiles}")
    print(f"Analysis Type: {recommendations.get('analysis_type', 'unknown')}")
    print(f"Detected Reaction Type: {recommendations.get('reaction_type', 'unknown')}")
    
    if 'ligand_recommendations' in recommendations:
        print(f"\nTop 3 Ligands:")
        for rec in recommendations['ligand_recommendations'][:3]:
            print(f"  • {rec['ligand']} (Score: {rec['compatibility_score']})")
    
    if 'solvent_recommendations' in recommendations:
        print(f"\nTop 3 Solvents:")
        for rec in recommendations['solvent_recommendations'][:3]:
            print(f"  • {rec['solvent']} ({rec['abbreviation']}) (Score: {rec['compatibility_score']})")
    
    if 'combined_conditions' in recommendations:
        print(f"\nTop Combined Condition:")
        if recommendations['combined_conditions']:
            best = recommendations['combined_conditions'][0]
            print(f"  Ligand: {best['ligand']}")
            print(f"  Solvent: {best['solvent']} ({best['solvent_abbreviation']})")
            print(f"  Combined Score: {best['combined_score']}")
            print(f"  Confidence: {best['recommendation_confidence']}")
