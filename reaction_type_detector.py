"""
Automatic Reaction Type Detection using rxn-insight
===================================================

This module integrates the rxn-insight library to automatically detect reaction types
and maps them to our internal reaction type nomenclature.

rxn-insight paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC10980627/
GitHub: https://github.com/mrodobbe/Rxn-INSIGHT
"""

from typing import Optional, Dict, Any
import logging

# Set up logging
logger = logging.getLogger(__name__)

# Mapping from rxn-insight classification to our internal reaction types
RXN_INSIGHT_MAPPING = {
    # Heteroatom Alkylation and Arylation -> C-N Coupling
    ("Heteroatom Alkylation and Arylation", "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)"): "C-N Coupling - Buchwald-Hartwig",
    ("Heteroatom Alkylation and Arylation", "Goldberg coupling aryl amine-aryl chloride"): "C-N Coupling - Ullmann", 
    ("Heteroatom Alkylation and Arylation", "Ullmann condensation with aryl halides"): "C-N Coupling - Ullmann",
    ("Heteroatom Alkylation and Arylation", "Chan-Lam coupling"): "C-N Oxidative Coupling - Chan-Lam",
    
    # C-C Coupling reactions
    ("C-C Coupling", "Suzuki coupling with boronic acids"): "C-C Coupling - Suzuki",
    ("C-C Coupling", "Suzuki coupling with boronic esters"): "C-C Coupling - Suzuki", 
    ("C-C Coupling", "Stille coupling"): "C-C Coupling - Stille",
    ("C-C Coupling", "Heck coupling"): "C-C Coupling - Heck",
    ("C-C Coupling", "Sonogashira coupling"): "C-C Coupling - Sonogashira",
    ("C-C Coupling", "Negishi coupling"): "C-C Coupling - Negishi",
    
    # Acylation -> Amide Formation
    ("Acylation", "Carboxylic acid with primary amine to amide"): "Amide Formation - Acid + Amine",
    ("Acylation", "Carboxylic acid with secondary amine to amide"): "Amide Formation - Acid + Amine",
    ("Acylation", "Acid chloride with primary amine to amide"): "Amide Formation",
    ("Acylation", "Acid chloride with secondary amine to amide"): "Amide Formation",
    
    # Other common reaction types that might appear
    ("Oxidation", None): "Oxidation",
    ("Reduction", None): "Reduction", 
    ("Substitution", None): "Substitution",
    ("Elimination", None): "Elimination",
    ("Cycloaddition", None): "Cycloaddition",
}

# Fallback mapping based on CLASS only (when NAME is not in the mapping)
CLASS_FALLBACK_MAPPING = {
    "Heteroatom Alkylation and Arylation": "C-N Coupling - Buchwald-Hartwig",  # Default to most common
    "C-C Coupling": "Cross-Coupling",  # Generic cross-coupling
    "Acylation": "Amide Formation",
    "Oxidation": "Oxidation", 
    "Reduction": "Reduction",
    "Substitution": "Substitution",
    "Elimination": "Elimination",
    "Cycloaddition": "Cycloaddition",
}


def detect_reaction_type(reaction_smiles: str) -> Dict[str, Any]:
    """
    Detect reaction type using rxn-insight and map to our internal nomenclature.
    
    Args:
        reaction_smiles: Reaction SMILES string
        
    Returns:
        Dict containing:
        - 'detected_type': Our internal reaction type string
        - 'confidence': Confidence level (high/medium/low)
        - 'rxn_insight_class': Original rxn-insight CLASS
        - 'rxn_insight_name': Original rxn-insight NAME  
        - 'functional_groups': Reactant -> Product functional group changes
        - 'error': Error message if detection failed
    """
    result = {
        'detected_type': None,
        'confidence': 'low',
        'rxn_insight_class': None,
        'rxn_insight_name': None,
        'functional_groups': None,
        'error': None
    }
    
    try:
        # Try to import and use rxn-insight
        from rxn_insight.reaction import Reaction
        
        # Create reaction object and get info
        rxn = Reaction(reaction_smiles)
        info = rxn.get_reaction_info()
        
        # Extract key information
        rxn_class = info.get('CLASS', '')
        rxn_name = info.get('NAME', '')
        fg_reactants = info.get('FG_REACTANTS', [])
        fg_products = info.get('FG_PRODUCTS', [])
        
        # Store raw rxn-insight results
        result['rxn_insight_class'] = rxn_class
        result['rxn_insight_name'] = rxn_name
        result['functional_groups'] = f"{fg_reactants} -> {fg_products}"
        
        # Map to our internal types
        mapping_key = (rxn_class, rxn_name)
        
        if mapping_key in RXN_INSIGHT_MAPPING:
            # Exact match found
            result['detected_type'] = RXN_INSIGHT_MAPPING[mapping_key]
            result['confidence'] = 'high'
            logger.info(f"Exact mapping found: {mapping_key} -> {result['detected_type']}")
            
        elif rxn_class in CLASS_FALLBACK_MAPPING:
            # Fallback to class-based mapping
            result['detected_type'] = CLASS_FALLBACK_MAPPING[rxn_class]
            result['confidence'] = 'medium'
            logger.info(f"Fallback mapping used: {rxn_class} -> {result['detected_type']}")
            
        else:
            # No mapping found
            result['detected_type'] = "Unknown"
            result['confidence'] = 'low'
            result['error'] = f"No mapping found for rxn-insight classification: {rxn_class} / {rxn_name}"
            logger.warning(result['error'])
            
    except ImportError:
        result['error'] = "rxn-insight library not available. Install with: pip install rxn-insight"
        logger.error(result['error'])
        
    except Exception as e:
        result['error'] = f"Error in rxn-insight detection: {str(e)}"
        logger.error(result['error'])
    
    return result


def get_supported_rxn_insight_types() -> Dict[str, list]:
    """
    Get all rxn-insight reaction types that we can map to our internal types.
    
    Returns:
        Dict with 'mapped' and 'unmapped' lists
    """
    mapped_types = list(RXN_INSIGHT_MAPPING.keys())
    
    return {
        'mapped_exact': [f"{cls} / {name}" for cls, name in mapped_types if name is not None],
        'mapped_fallback': list(CLASS_FALLBACK_MAPPING.keys()),
        'total_mappings': len(RXN_INSIGHT_MAPPING) + len(CLASS_FALLBACK_MAPPING)
    }


if __name__ == "__main__":
    # Test the detection system
    test_reactions = [
        "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1",  # Buchwald-Hartwig
        "Clc1nc2ccccc2[nH]1.Nc1ccccc1>>c1ccc(Nc2nc3ccccc3[nH]2)cc1",  # User example
        "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1",  # Suzuki
        "O=C(O)c1ccccc1.NCc1ccccc1>>O=C(NCc1ccccc1)c1ccccc1",  # Amide
    ]
    
    print("Testing reaction type detection...")
    print("=" * 60)
    
    for smiles in test_reactions:
        print(f"\nSMILES: {smiles}")
        result = detect_reaction_type(smiles)
        
        if result['error']:
            print(f"ERROR: {result['error']}")
        else:
            print(f"Detected Type: {result['detected_type']} (confidence: {result['confidence']})")
            print(f"rxn-insight: {result['rxn_insight_class']} / {result['rxn_insight_name']}")
            print(f"Functional Groups: {result['functional_groups']}")
        print("-" * 40)
    
    # Show supported mappings
    print(f"\nSupported mappings:")
    mappings = get_supported_rxn_insight_types()
    print(f"Total mappings: {mappings['total_mappings']}")
