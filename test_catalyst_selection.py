#!/usr/bin/env python3
"""Test the catalyst selection workflow for C-N coupling reactions"""

from enhanced_recommendation_engine import create_recommendation_engine

def test_catalyst_selection_workflow():
    """Demonstrate the new catalyst selection workflow"""
    
    print("🧪 CATALYST SELECTION WORKFLOW DEMO")
    print("=" * 50)
    print()
    
    # Test C-N coupling reaction that should trigger catalyst selection
    reaction_smiles = 'Brc1ccc(C(F)(F)F)cc1.Nc1ccccc1>>FC(F)(F)c1ccc(Nc2ccccc2)cc1'
    
    print(f"Testing reaction: {reaction_smiles}")
    print()
    
    engine = create_recommendation_engine()
    
    # Step 1: Initial detection
    print("Step 1: Initial reaction type detection")
    print("-" * 40)
    
    result = engine.get_recommendations(reaction_smiles, "Auto detect reaction type")
    
    if result.get('status') == 'needs_catalyst_selection':
        print("✅ Catalyst selection needed detected!")
        print(f"   rxn-insight detected: {result.get('rxn_insight_class')} / {result.get('rxn_insight_name')}")
        print(f"   Message: {result.get('message')}")
        print(f"   Available catalysts: {result.get('available_catalysts')}")
        print()
        print("   Catalyst options:")
        for catalyst, reaction_type in result.get('catalyst_options', {}).items():
            print(f"      {catalyst}: {reaction_type}")
        print()
        
        # Step 2: Test different catalyst selections
        print("Step 2: Testing catalyst selections")
        print("-" * 40)
        
        catalysts_to_test = ['Pd', 'Cu', 'Ni']
        
        for catalyst in catalysts_to_test:
            print(f"Testing with {catalyst} catalyst:")
            
            catalyst_result = engine.get_recommendations_with_catalyst(
                reaction_smiles, 
                catalyst, 
                "Auto detect reaction type"
            )
            
            if 'rxn_insight_detection' in catalyst_result:
                detection = catalyst_result['rxn_insight_detection']
                print(f"   Final reaction type: {catalyst_result.get('reaction_type')}")
                print(f"   Catalyst selected: {detection.get('catalyst_selected')}")
                print(f"   Status: {catalyst_result.get('status')}")
            else:
                print(f"   Result: {catalyst_result.get('status')} - {catalyst_result.get('error', 'Unknown error')}")
            print()
            
    else:
        print("❌ Expected catalyst selection request, but got:")
        print(f"   Status: {result.get('status')}")
        print(f"   Reaction type: {result.get('reaction_type')}")
        if 'rxn_insight_detection' in result:
            detection = result['rxn_insight_detection']
            print(f"   rxn-insight: {detection.get('rxn_insight_detected')}")
        print()

if __name__ == "__main__":
    test_catalyst_selection_workflow()
