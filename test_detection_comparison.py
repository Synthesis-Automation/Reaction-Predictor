#!/usr/bin/env python3
"""Test rxn-insight detection vs our mapping for C-N coupling reactions"""

from enhanced_recommendation_engine import create_recommendation_engine

def test_cn_coupling_detection():
    """Test C-N coupling detection for Ullmann vs Buchwald differentiation"""
    
    # Test reactions that could be either Ullmann or Buchwald
    test_reactions = [
        {
            'name': 'Simple aniline arylation',
            'smiles': 'Brc1ccccc1.Nc1ccccc1>>Nc1ccc(c2ccccc2)cc1'
        },
        {
            'name': 'Heteroaryl coupling',
            'smiles': 'Clc1ccncc1.Nc1ccccc1>>Nc1ccc(c2ccncc2)cc1'
        },
        {
            'name': 'Amide formation (should not be C-N coupling)',
            'smiles': 'O=C(O)c1ccccc1.Nc1ccccc1>>O=C(Nc1ccccc1)c1ccccc1'
        },
        {
            'name': 'Electron-deficient aryl halide (typical Buchwald)',
            'smiles': 'Brc1ccc(C(F)(F)F)cc1.Nc1ccccc1>>Nc1ccc(c2ccc(C(F)(F)F)cc2)cc1'
        }
    ]
    
    engine = create_recommendation_engine()
    
    print("=== C-N COUPLING DETECTION COMPARISON ===\n")
    
    for i, reaction in enumerate(test_reactions, 1):
        print(f"{i}. {reaction['name']}")
        print(f"   SMILES: {reaction['smiles']}")
        
        # Test with auto-detection
        result = engine.get_recommendations(reaction['smiles'], "Auto detect reaction type")
        
        if 'rxn_insight_detection' in result:
            detection = result['rxn_insight_detection']
            print(f"   rxn-insight: {detection.get('rxn_insight_detected', 'N/A')}")
            print(f"   Our mapping: {detection.get('mapped_to', 'N/A')}")
            print(f"   Confidence: {detection.get('confidence', 'N/A')}")
        else:
            print("   No rxn-insight detection info available")
        
        print(f"   Final type: {result.get('reaction_type', 'Unknown')}")
        print(f"   Status: {result.get('status', 'Unknown')}")
        print()

if __name__ == "__main__":
    test_cn_coupling_detection()
