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
            "Chan-Lam Coupling": "Cross-Coupling",
            "C-N Coupling - Chan-Lam": "Cross-Coupling",
            # Ullmann variants
            "Ullmann Ether Synthesis": "Ullmann",
            "Ullmann Reaction": "Ullmann",
            "C-N Coupling - Ullmann": "Ullmann",
            "C-O Coupling - Ullmann Ether": "Ullmann",
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

            # Evidence-aware context: mine dataset for ligands used in similar reactions (if available)
            evidence_ligands = self._harvest_evidence_ligands(reaction_type)

            # Get top ligands for this reaction type, with evidence-aware boost
            ligands = recommend_ligands_for_reaction(
                reaction_type=reaction_type,
                top_n=5,
                min_compatibility=0.4,
                evidence_ligands=evidence_ligands or None
            )
            recommendations['ligand_recommendations'] = ligands
            
            # Get top solvents for this reaction type
            solvents = recommend_solvents_for_reaction(
                reaction_type=reaction_type,
                top_n=5,
                min_compatibility=0.4
            )
            recommendations['solvent_recommendations'] = solvents
            
            # Try to get base recommendations if available
            try:
                from reagents.base import recommend_bases_for_reaction  # type: ignore
                bases = recommend_bases_for_reaction(
                    reaction_type=reaction_type,
                    top_n=5,
                    min_compatibility=0.4
                )
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
                    'reaction_types_supported': ['Cross-Coupling', 'Hydrogenation', 'Metathesis', 'C-H_Activation', 'Carbonylation']
                }
            except:
                pass
            
        except Exception as e:
            recommendations['error'] = f"Enhanced recommendation error: {str(e)}"
        
        return recommendations

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
            for fname in os.listdir(data_dir):
                if not fname.lower().endswith('.csv'):
                    continue
                path = os.path.join(data_dir, fname)
                try:
                    import csv
                    with open(path, 'r', encoding='utf-8') as f:
                        reader = csv.DictReader(f)
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
                                name = it
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
