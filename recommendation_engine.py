"""
Reaction Condition Recommendation Engine
=======================================

A modular recommendation system that can be extended for different reaction types.
Currently focused on Buchwald-Hartwig amination but designed for future expansion.
"""

import pandas as pd
import numpy as np
from abc import ABC, abstractmethod
from typing import Dict, List, Optional, Tuple
import json
import os
from dataset_registry import resolve_dataset_path

class ReactionRecommender(ABC):
    """Base class for reaction condition recommenders"""
    
    def __init__(self, name: str):
        self.name = name
        self.dataset = None
        self.is_loaded = False
    
    @abstractmethod
    def can_handle_reaction(self, reaction_smiles: str, reaction_type: str = None) -> bool:
        """Check if this recommender can handle the given reaction"""
        pass
    
    @abstractmethod
    def get_recommendations(self, reaction_smiles: str, top_k: int = 5) -> Dict:
        """Get condition recommendations for a reaction"""
        pass
    
    @abstractmethod
    def load_data(self, dataset_path: str) -> bool:
        """Load reaction dataset"""
        pass

class CatalystFamily:
    """Represents a family of related catalysts"""
    
    def __init__(self, family_id: str, name: str):
        self.family_id = family_id
        self.name = name
        self.catalysts = []
        self.ligands = []
        self.properties = {}
        self.performance_stats = {}
    
    def add_catalyst(self, catalyst: str):
        if catalyst not in self.catalysts:
            self.catalysts.append(catalyst)
    
    def add_ligand(self, ligand: str):
        if ligand not in self.ligands:
            self.ligands.append(ligand)
    
    def get_similarity_score(self, other_family) -> float:
        """Calculate similarity with another catalyst family"""
        if not isinstance(other_family, CatalystFamily):
            return 0.0
        
        # Simple similarity based on shared components
        shared_ligands = set(self.ligands) & set(other_family.ligands)
        shared_catalysts = set(self.catalysts) & set(other_family.catalysts)
        
        total_ligands = len(set(self.ligands) | set(other_family.ligands))
        total_catalysts = len(set(self.catalysts) | set(other_family.catalysts))
        
        if total_ligands == 0 and total_catalysts == 0:
            return 0.0
        
        ligand_sim = len(shared_ligands) / total_ligands if total_ligands > 0 else 0
        catalyst_sim = len(shared_catalysts) / total_catalysts if total_catalysts > 0 else 0
        
        # Weight ligands more heavily (60%) than catalysts (40%)
        return ligand_sim * 0.6 + catalyst_sim * 0.4

class BuchwaldHartwigRecommender(ReactionRecommender):
    """Recommender for Buchwald-Hartwig amination reactions"""
    
    def __init__(self):
        super().__init__("Buchwald-Hartwig Amination")
        self.catalyst_families = {}
        self.condition_performance = {}
        self.ligand_families_config = None
        
    def can_handle_reaction(self, reaction_smiles: str, reaction_type: str = None) -> bool:
        """Check if this is a Buchwald-Hartwig reaction"""
        
        # Check by reaction type first
        if reaction_type and "buchwald" in reaction_type.lower():
            return True
        
        # Simple pattern matching for C-N coupling
        try:
            if ">>" not in reaction_smiles:
                return False
            
            reactants, products = reaction_smiles.split(">>")
            
            # Look for patterns: aryl halide + amine -> aryl amine
            has_halide = any(x in reactants for x in ['Br', 'Cl', 'I'])
            has_nitrogen = 'N' in reactants and 'N' in products
            
            return has_halide and has_nitrogen
        except:
            return False
    
    def load_data(self, dataset_path: str) -> bool:
        """Load Buchwald dataset and ligand families configuration"""
        try:
            if not os.path.exists(dataset_path):
                return False
            
            # Load main dataset
            self.dataset = pd.read_csv(dataset_path)
            
            # Load ligand families configuration
            self._load_ligand_families_config()
            
            # Process data
            self._clean_data()
            self._create_catalyst_families()
            self._calculate_performance_stats()
            self.is_loaded = True
            return True
            
        except Exception as e:
            print(f"Error loading dataset: {e}")
            return False
    
    def _clean_data(self):
        """Clean and standardize the dataset"""
        if self.dataset is None:
            return
        
        # Clean string columns, handling JSON arrays
        string_cols = ['CoreGeneric', 'Ligand', 'Solvent']
        for col in string_cols:
            if col in self.dataset.columns:
                # Handle JSON array format like ["XPhos"]
                cleaned_values = []
                for value in self.dataset[col]:
                    try:
                        if pd.isna(value) or value == '':
                            cleaned_values.append('none')
                        elif isinstance(value, str):
                            # Try to parse JSON array
                            if value.startswith('[') and value.endswith(']'):
                                import json
                                parsed = json.loads(value)
                                if isinstance(parsed, list) and len(parsed) > 0:
                                    cleaned_values.append(str(parsed[0]).strip())
                                else:
                                    cleaned_values.append('none')
                            else:
                                # Plain string
                                cleaned_values.append(str(value).strip())
                        else:
                            cleaned_values.append(str(value).strip())
                    except (json.JSONDecodeError, ValueError):
                        # Fallback: strip brackets and quotes manually
                        cleaned_val = str(value).strip('[]"').strip()
                        cleaned_values.append(cleaned_val if cleaned_val else 'none')
                
                self.dataset[f'{col}_clean'] = cleaned_values
        
        # Create combination strings
        self.dataset['catalyst_combination'] = (
            self.dataset['CoreGeneric_clean'].astype(str) + '/' + 
            self.dataset['Ligand_clean'].astype(str)
        )
    
    def _load_ligand_families_config(self):
        """Load ligand families configuration from CSV"""
        try:
            # Look for ligand families config in data directory
            config_path = os.path.join(os.path.dirname(__file__), 'data', 'ligand_families.csv')
            
            if os.path.exists(config_path):
                self.ligand_families_config = pd.read_csv(config_path)
                print(f"Loaded {len(self.ligand_families_config)} ligand families from config")
            else:
                print(f"Ligand families config not found at {config_path}, using defaults")
                self._create_default_families_config()
                
        except Exception as e:
            print(f"Error loading ligand families config: {e}")
            self._create_default_families_config()
    
    def _create_default_families_config(self):
        """Create default ligand families if config file not found"""
        default_data = [
            ['xphos_family', 'XPhos Family', 'XPhos - Excellent for electron-deficient aryl halides', 'XPhos', 'XPhos,X-Phos', 1.1, 1],
            ['sphos_family', 'SPhos Family', 'SPhos - Good balance of activity and selectivity', 'SPhos', 'SPhos,S-Phos', 1.05, 2],
            ['brettphos_family', 'BrettPhos Family', 'BrettPhos - Highly active for challenging substrates', 'BrettPhos', 'BrettPhos,Brett-Phos', 1.15, 3],
            ['tbu_family', 'tBuXPhos Family', 'tBuXPhos - Sterically demanding for selective reactions', 'tBuXPhos', 'tBuXPhos,tBu-XPhos', 1.0, 4],
            ['other', 'Alternative Conditions', 'Alternative Conditions (No Specific Ligand)', 'none,nan,', 'none,nan,null', 0.7, 99]
        ]
        
        self.ligand_families_config = pd.DataFrame(default_data, columns=[
            'family_id', 'family_name', 'description', 'ligands', 'aliases', 'performance_modifier', 'priority_rank'
        ])
    
    def _create_catalyst_families(self):
        """Create catalyst families based on CSV configuration"""
        if self.dataset is None or self.ligand_families_config is None:
            return

        # Create family objects from configuration
        for _, config_row in self.ligand_families_config.iterrows():
            family_id = config_row['family_id']
            family_name = config_row['family_name']
            description = config_row['description']
            ligands_str = config_row['ligands']
            
            # Parse ligands (comma-separated)
            ligands = [ligand.strip() for ligand in ligands_str.split(',') if ligand.strip()]
            
            # Create family object
            family = CatalystFamily(family_id, description)
            
            # Add ligands to family
            for ligand in ligands:
                family.add_ligand(ligand)
            
            # Add aliases if available
            if 'aliases' in config_row and pd.notna(config_row['aliases']):
                aliases = [alias.strip() for alias in str(config_row['aliases']).split(',') if alias.strip()]
                for alias in aliases:
                    family.add_ligand(alias)
            
            # Store additional properties
            family.properties = {
                'ligand_type': config_row.get('ligand_type', 'unknown'),
                'steric_bulk': config_row.get('steric_bulk', 'unknown'),
                'electronic_property': config_row.get('electronic_property', 'unknown'),
                'coordination_mode': config_row.get('coordination_mode', 'unknown'),
                'typical_applications': config_row.get('typical_applications', 'general'),
                'cost_category': config_row.get('cost_category', 'medium'),
                'performance_modifier': config_row.get('performance_modifier', 1.0),
                'priority_rank': config_row.get('priority_rank', 99)
            }
            
            # Find all catalysts used with these ligands in the dataset
            all_family_ligands = ligands + (aliases if 'aliases' in config_row and pd.notna(config_row['aliases']) else [])
            family_data = self.dataset[
                self.dataset['Ligand_clean'].isin(all_family_ligands)
            ]
            
            for catalyst in family_data['CoreGeneric_clean'].unique():
                if pd.notna(catalyst):
                    family.add_catalyst(str(catalyst))
            
            self.catalyst_families[family_id] = family
    
    def _calculate_performance_stats(self):
        """Calculate performance statistics for each condition combination"""
        if self.dataset is None:
            return
        
        # Group by catalyst combination
        grouped = self.dataset.groupby('catalyst_combination')
        
        for combination, group in grouped:
            yields = group['Yield_%'].dropna()
            
            if len(yields) > 0:
                stats = {
                    'avg_yield': yields.mean(),
                    'max_yield': yields.max(),
                    'min_yield': yields.min(),
                    'success_rate': (yields >= 80).mean(),
                    'reaction_count': len(yields),
                    'confidence': min(len(yields) / 10.0, 1.0),  # Confidence based on data volume
                    'yield_std': yields.std()
                }
                
                self.condition_performance[combination] = stats
    
    def get_recommendations(self, reaction_smiles: str, top_k: int = 5) -> Dict:
        """Get recommendations for Buchwald-Hartwig reaction"""
        
        if not self.is_loaded:
            return {
                'error': 'Dataset not loaded',
                'recommendations': [],
                'analysis_type': 'error'
            }
        
        try:
            # Get structural similarity matches (simplified for now)
            similar_conditions = self._find_similar_conditions(reaction_smiles)
            
            # Group by catalyst families
            family_recommendations = self._group_by_families(similar_conditions)
            
            # Rank and format recommendations
            ranked_recommendations = self._rank_recommendations(family_recommendations)
            
            return {
                'reaction_smiles': reaction_smiles,
                'reaction_type': self.name,
                'analysis_type': 'buchwald_hartwig',
                'status': 'Recommendations generated successfully',
                'recommendations': ranked_recommendations[:top_k],
                'total_found': len(ranked_recommendations),
                'catalyst_families': len(self.catalyst_families),
                'dataset_size': len(self.dataset) if self.dataset is not None else 0
            }
            
        except Exception as e:
            return {
                'error': f'Recommendation failed: {str(e)}',
                'recommendations': [],
                'analysis_type': 'error'
            }
    
    def _find_similar_conditions(self, reaction_smiles: str) -> List[Dict]:
        """Find conditions for similar reactions (simplified approach)"""
        
        # For now, return top performing conditions
        # In future versions, this will use molecular similarity
        
        similar_conditions = []
        
        for combination, stats in self.condition_performance.items():
            if stats['reaction_count'] >= 2:  # Require at least 2 data points
                
                try:
                    catalyst, ligand = combination.split('/', 1)
                except ValueError:
                    continue
                
                condition = {
                    'catalyst': catalyst,
                    'ligand': ligand,
                    'combination': combination,
                    'performance': stats,
                    'structural_similarity': 0.5  # Placeholder for now
                }
                
                similar_conditions.append(condition)
        
        return similar_conditions
    
    def _group_by_families(self, conditions: List[Dict]) -> Dict[str, List[Dict]]:
        """Group conditions by catalyst families"""
        
        family_groups = {}
        
        for condition in conditions:
            ligand = condition['ligand']
            family_id = self._get_ligand_family(ligand)
            
            if family_id not in family_groups:
                family_groups[family_id] = []
            
            family_groups[family_id].append(condition)
        
        return family_groups
    
    def _get_ligand_family(self, ligand: str) -> str:
        """Get the family ID for a ligand using CSV configuration"""
        
        if not ligand or ligand.lower() in ['none', 'nan', '']:
            return 'other'

        # Direct match first - check all families from CSV config
        for family_id, family in self.catalyst_families.items():
            if ligand in family.ligands:
                return family_id

        # Check for partial matches using CSV data
        if self.ligand_families_config is not None:
            ligand_lower = ligand.lower()
            
            for _, config_row in self.ligand_families_config.iterrows():
                family_id = config_row['family_id']
                
                # Check main ligands
                ligands_str = config_row['ligands']
                ligands = [l.strip().lower() for l in ligands_str.split(',') if l.strip()]
                
                # Check aliases
                if 'aliases' in config_row and pd.notna(config_row['aliases']):
                    aliases_str = str(config_row['aliases'])
                    aliases = [a.strip().lower() for a in aliases_str.split(',') if a.strip()]
                    ligands.extend(aliases)
                
                # Check for partial matches
                for known_ligand in ligands:
                    if (ligand_lower in known_ligand) or (known_ligand in ligand_lower):
                        return family_id
        
        # Fallback to pattern matching for backwards compatibility
        return self._fallback_ligand_classification(ligand)
    
    def _fallback_ligand_classification(self, ligand: str) -> str:
        """Fallback ligand classification using pattern matching"""
        ligand_lower = ligand.lower()
        
        # Specific ligand families
        if 'xphos' in ligand_lower and 'tbux' not in ligand_lower:
            return 'xphos_family'
        elif 'sphos' in ligand_lower:
            return 'sphos_family'
        elif 'brettphos' in ligand_lower or 'brett' in ligand_lower:
            return 'brettphos_family'
        elif 'tbux' in ligand_lower or 'tbu' in ligand_lower:
            return 'tbu_family'
        elif 'ruphos' in ligand_lower:
            return 'ruphos_family'
        elif 'davephos' in ligand_lower:
            return 'davephos_family'
        elif 'johnphos' in ligand_lower:
            return 'johnphos_family'
        elif 'cyjohnphos' in ligand_lower:
            return 'cyjohnphos_family'
        elif 'dppf' in ligand_lower:
            return 'dppf_family'
        elif 'dppp' in ligand_lower:
            return 'dppp_family'
        elif 'dppe' in ligand_lower:
            return 'dppe_family'
        elif 'binap' in ligand_lower:
            return 'binap_family'
        elif 'xantphos' in ligand_lower:
            return 'xantphos_family'
        elif 'pph3' in ligand_lower or 'triphenylphosphine' in ligand_lower:
            return 'pph3_family'
        elif 'tolyl' in ligand_lower:
            return 'ptolyl3_family'
        elif 'pcy3' in ligand_lower or 'tricyclohexyl' in ligand_lower:
            return 'pcy3_family'
        elif 'ipr' in ligand_lower and 'nhc' in ligand_lower or 'carbene' in ligand_lower:
            return 'ipr_family'
        elif 'imes' in ligand_lower:
            return 'imes_family'
        elif 'sipr' in ligand_lower:
            return 'sipr_family'
        elif 'phen' in ligand_lower:
            return 'phen_family'
        elif 'bipy' in ligand_lower or 'bipyridine' in ligand_lower:
            return 'bipy_family'
        else:
            return 'other'

    def _rank_recommendations(self, family_groups: Dict[str, List[Dict]]) -> List[Dict]:
        """Rank and format final recommendations"""
        
        recommendations = []
        
        for family_id, conditions in family_groups.items():
            if not conditions:
                continue
            
            # Calculate family-level statistics
            family_stats = self._calculate_family_stats(conditions)
            
            # Get best representative condition
            best_condition = max(conditions, key=lambda x: x['performance']['avg_yield'])
            
            # Get alternative conditions in the same family
            alternatives = [c for c in conditions if c != best_condition][:3]
            
            family_rec = {
                'family_id': family_id,
                'family_name': self.catalyst_families.get(family_id, CatalystFamily(family_id, 'Unknown')).name,
                'recommended_catalyst': best_condition['catalyst'],
                'recommended_ligand': best_condition['ligand'],
                'performance': family_stats,
                'alternatives': [
                    {
                        'catalyst': alt['catalyst'],
                        'ligand': alt['ligand'],
                        'avg_yield': alt['performance']['avg_yield'],
                        'success_rate': alt['performance']['success_rate']
                    }
                    for alt in alternatives
                ],
                'confidence_score': family_stats['confidence'],
                'ranking_score': self._calculate_ranking_score(family_stats, family_id)
            }
            
            recommendations.append(family_rec)
        
        # Sort by ranking score
        recommendations.sort(key=lambda x: x['ranking_score'], reverse=True)
        
        return recommendations
    
    def _calculate_family_stats(self, conditions: List[Dict]) -> Dict:
        """Calculate aggregate statistics for a catalyst family"""
        
        all_yields = []
        all_success_rates = []
        total_reactions = 0
        
        for condition in conditions:
            perf = condition['performance']
            weight = perf['reaction_count']
            
            all_yields.extend([perf['avg_yield']] * weight)
            all_success_rates.append(perf['success_rate'])
            total_reactions += perf['reaction_count']
        
        return {
            'avg_yield': np.mean(all_yields) if all_yields else 0,
            'max_yield': max(c['performance']['max_yield'] for c in conditions),
            'success_rate': np.mean(all_success_rates) if all_success_rates else 0,
            'total_reactions': total_reactions,
            'confidence': min(total_reactions / 20.0, 1.0),
            'num_variations': len(conditions)
        }
    
    def _calculate_ranking_score(self, family_stats: Dict, family_id: str = None) -> float:
        """Calculate a ranking score for a catalyst family using CSV properties"""
        
        # Base performance score
        yield_score = family_stats['avg_yield'] * 0.4
        success_score = family_stats['success_rate'] * 100 * 0.3
        confidence_score = family_stats['confidence'] * 100 * 0.2
        variety_score = min(family_stats['num_variations'] / 5.0, 1.0) * 100 * 0.1
        
        base_score = yield_score + success_score + confidence_score + variety_score
        
        # Apply family-specific modifiers from CSV
        if family_id and family_id in self.catalyst_families:
            family = self.catalyst_families[family_id]
            performance_modifier = family.properties.get('performance_modifier', 1.0)
            base_score *= performance_modifier
            
            # Cost bonus for cost-effective ligands
            cost_category = family.properties.get('cost_category', 'medium')
            if cost_category == 'very_low':
                base_score *= 1.05
            elif cost_category == 'low':
                base_score *= 1.02
            elif cost_category == 'high':
                base_score *= 0.98  # Small penalty for expensive ligands
        
        # Apply penalty for "other" family (less informative ligands)
        if family_id == 'other':
            base_score *= 0.7  # 30% penalty for "other" category
        
        return base_score

class GeneralRecommender(ReactionRecommender):
    """Fallback recommender for general reactions"""
    
    def __init__(self):
        super().__init__("General Reaction")
    
    def can_handle_reaction(self, reaction_smiles: str, reaction_type: str = None) -> bool:
        """Can handle any reaction as fallback"""
        return True
    
    def load_data(self, dataset_path: str) -> bool:
        """No specific dataset needed for general recommender"""
        self.is_loaded = True
        return True
    
    def get_recommendations(self, reaction_smiles: str, top_k: int = 5) -> Dict:
        """Provide general recommendations"""
        
        return {
            'reaction_smiles': reaction_smiles,
            'reaction_type': self.name,
            'analysis_type': 'general',
            'recommendations': [],
            'message': 'No specific recommendations available. Consider literature search for similar reactions.',
            'suggestions': [
                'Check for similar reaction types in literature',
                'Consider standard conditions for the reaction class',
                'Start with mild conditions and optimize',
                'Use common catalysts for the reaction type'
            ]
        }

class RecommendationEngine:
    """Main recommendation engine that coordinates different recommenders"""
    
    def __init__(self):
        self.recommenders = []
        self.default_recommender = GeneralRecommender()
        self.default_recommender.load_data("")
        
    def add_recommender(self, recommender: ReactionRecommender, dataset_path: str = None):
        """Add a specific recommender to the engine"""
        
        if dataset_path and os.path.exists(dataset_path):
            if recommender.load_data(dataset_path):
                self.recommenders.append(recommender)
                print(f"Added recommender: {recommender.name}")
            else:
                print(f"Failed to load data for recommender: {recommender.name}")
        else:
            if recommender.load_data(""):
                self.recommenders.append(recommender)
                print(f"Added recommender: {recommender.name}")
    
    def get_recommendations(self, reaction_smiles: str, reaction_type: str = None, top_k: int = 5) -> Dict:
        """Get recommendations using the most appropriate recommender"""
        
        # Find the best recommender for this reaction
        for recommender in self.recommenders:
            if recommender.can_handle_reaction(reaction_smiles, reaction_type):
                return recommender.get_recommendations(reaction_smiles, top_k)
        
        # Fall back to general recommender
        return self.default_recommender.get_recommendations(reaction_smiles, top_k)
    
    def get_available_recommenders(self) -> List[str]:
        """Get list of available recommenders"""
        return [rec.name for rec in self.recommenders] + [self.default_recommender.name]

# Convenience function for easy setup
def create_recommendation_engine(data_dir: str = "data") -> RecommendationEngine:
    """Create and configure the recommendation engine"""
    
    engine = RecommendationEngine()
    
    # Add Buchwald-Hartwig recommender if data is available (resolve via registry)
    base_dir = os.path.dirname(__file__)
    buchwald_path = resolve_dataset_path("C-N Coupling - Buchwald-Hartwig", base_dir) or os.path.join(data_dir, "buchwald_reactions.csv")
    if os.path.exists(buchwald_path):
        buchwald_rec = BuchwaldHartwigRecommender()
        engine.add_recommender(buchwald_rec, buchwald_path)
    
    # Future: Add other reaction-specific recommenders here
    # suzuki_rec = SuzukiRecommender()
    # engine.add_recommender(suzuki_rec, "data/suzuki_reactions.csv")
    
    return engine
