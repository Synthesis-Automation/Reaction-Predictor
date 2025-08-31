#!/usr/bin/env python3
"""Demonstrate the rxn-insight detection comparison feature"""

from enhanced_recommendation_engine import create_recommendation_engine

def demonstrate_detection_comparison():
    """Show the rxn-insight vs our mapping comparison"""
    
    print("🔍 RXN-INSIGHT DETECTION COMPARISON DEMO")
    print("=" * 50)
    print()
    
    test_reactions = [
        {
            'name': 'C-N Coupling (detected as Buchwald-Hartwig by rxn-insight)',
            'smiles': 'Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1'
        },
        {
            'name': 'Amide Formation (should be detected correctly)',
            'smiles': 'O=C(O)c1ccccc1.Nc1ccccc1>>O=C(Nc1ccccc1)c1ccccc1'
        }
    ]
    
    engine = create_recommendation_engine()
    
    for i, reaction in enumerate(test_reactions, 1):
        print(f"{i}. {reaction['name']}")
        print(f"   SMILES: {reaction['smiles']}")
        print()
        
        # Get recommendations with auto-detection
        result = engine.get_recommendations(reaction['smiles'], "Auto detect reaction type")
        
        # Show the comparison
        if 'rxn_insight_detection' in result:
            detection = result['rxn_insight_detection']
            print("   🤖 rxn-insight Detection:")
            print(f"      Raw Detection: {detection.get('rxn_insight_detected', 'N/A')}")
            print(f"      Confidence: {detection.get('confidence', 'N/A')}")
            print()
            print("   🔄 Our System Mapping:")
            print(f"      Mapped to: {detection.get('mapped_to', 'N/A')}")
            print(f"      Final Type: {result.get('reaction_type', 'Unknown')}")
        else:
            print("   ❌ No rxn-insight detection info available")
        
        print(f"   📊 Status: {result.get('status', 'Unknown')}")
        print()
        print("-" * 40)
        print()

if __name__ == "__main__":
    demonstrate_detection_comparison()
