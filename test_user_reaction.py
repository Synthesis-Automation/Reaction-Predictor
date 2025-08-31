#!/usr/bin/env python3
"""Test the exact reaction from the user's issue"""

from enhanced_recommendation_engine import create_recommendation_engine

def test_user_reaction():
    """Test the exact reaction that showed 'Not specified' in GUI"""
    
    reaction_smiles = 'Brc1ccc(C#N)cc1.Nc1ccccc1>>N#Cc1ccc(Nc2ccccc2)cc1'
    
    print("Testing exact user reaction:")
    print(f"SMILES: {reaction_smiles}")
    print()
    
    engine = create_recommendation_engine()
    
    # Test with "Auto detect reaction type" (like GUI default)
    result = engine.get_recommendations(reaction_smiles, "Auto detect reaction type")
    
    print(f"Status: {result.get('status')}")
    print(f"Analysis Type: {result.get('analysis_type')}")
    
    if result.get('analysis_type') == 'catalyst_selection_needed':
        print("✅ Catalyst selection needed - this should show in GUI!")
        print(f"rxn-insight: {result.get('rxn_insight_class')} / {result.get('rxn_insight_name')}")
        print(f"Message: {result.get('message')}")
        print("Catalyst options:")
        for cat, rxn_type in result.get('catalyst_options', {}).items():
            print(f"  {cat}: {rxn_type}")
    else:
        print(f"❌ Got analysis_type: {result.get('analysis_type')}")
        print(f"Reaction type: {result.get('reaction_type')}")
        print(f"Error: {result.get('error')}")
        
        if 'rxn_insight_detection' in result:
            det = result['rxn_insight_detection']
            print(f"rxn-insight detected: {det.get('rxn_insight_detected')}")

if __name__ == "__main__":
    test_user_reaction()
