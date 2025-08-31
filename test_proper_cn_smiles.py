#!/usr/bin/env python3
"""Test with better C-N coupling SMILES for proper rxn-insight detection"""

from enhanced_recommendation_engine import create_recommendation_engine

def test_proper_cn_coupling_smiles():
    """Test with SMILES that should properly trigger C-N coupling detection"""
    
    # Test reactions with proper C-N coupling SMILES
    test_reactions = [
        {
            'name': 'Classic Buchwald-Hartwig (should detect as C-N)',
            'smiles': 'Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1'
        },
        {
            'name': 'Ullmann-type with simple aniline',
            'smiles': 'Ic1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1'
        },
        {
            'name': 'Buchwald with electron-poor aryl bromide',
            'smiles': 'Brc1ccc(C#N)cc1.Nc1ccccc1>>N#Cc1ccc(Nc2ccccc2)cc1'
        },
        {
            'name': 'Primary amine + aryl chloride',
            'smiles': 'Clc1ccccc1.CCN>>CCNc1ccccc1'
        }
    ]
    
    engine = create_recommendation_engine()
    
    print("=== PROPER C-N COUPLING SMILES TESTING ===\n")
    
    for i, reaction in enumerate(test_reactions, 1):
        print(f"{i}. {reaction['name']}")
        print(f"   SMILES: {reaction['smiles']}")
        
        # Test with auto-detection
        result = engine.get_recommendations(reaction['smiles'], "Auto detect reaction type")
        
        if 'rxn_insight_detection' in result:
            detection = result['rxn_insight_detection']
            print(f"   rxn-insight: {detection.get('rxn_insight_detected', 'N/A')}")
            print(f"   Raw detection: {detection.get('raw_detected_type', 'N/A')}")
            print(f"   Our mapping: {detection.get('mapped_to', 'N/A')}")
            print(f"   Confidence: {detection.get('confidence', 'N/A')}")
        else:
            print("   No rxn-insight detection info available")
        
        print(f"   Final type: {result.get('reaction_type', 'Unknown')}")
        print(f"   Status: {result.get('status', 'Unknown')}")
        print()

if __name__ == "__main__":
    test_proper_cn_coupling_smiles()
