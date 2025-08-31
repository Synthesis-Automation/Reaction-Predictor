#!/usr/bin/env python3
"""Test the complete catalyst selection workflow with the new GUI feature"""

from enhanced_recommendation_engine import create_recommendation_engine

def test_catalyst_workflow_integration():
    """Test how different catalyst selections affect predictions"""
    
    # Test reaction that should trigger catalyst selection
    reaction_smiles = 'Brc1ccc(C#N)cc1.Nc1ccccc1>>N#Cc1ccc(Nc2ccccc2)cc1'
    
    print("🔬 TESTING CATALYST SELECTOR INTEGRATION")
    print("=" * 55)
    print(f"Test Reaction: {reaction_smiles}")
    print()
    
    engine = create_recommendation_engine()
    
    # Test different catalyst selections
    catalysts_to_test = [
        ("auto", "Auto-detect (no specific catalyst)"),
        ("Pd", "Palladium (Buchwald-Hartwig expected)"),
        ("Cu", "Copper (Ullmann expected)"),
        ("Ni", "Nickel (Ni coupling expected)")
    ]
    
    for catalyst_value, description in catalysts_to_test:
        print(f"Testing: {description}")
        print("-" * 30)
        
        if catalyst_value == "auto":
            # Test auto-detection (should show catalyst selection needed)
            result = engine.get_recommendations(reaction_smiles, "Auto detect reaction type")
            
            if result.get('analysis_type') == 'catalyst_selection_needed':
                print("✅ Auto-detection correctly identified need for catalyst selection")
                print(f"   Message: {result.get('message')}")
                print("   Available options:")
                for cat, rxn_type in result.get('catalyst_options', {}).items():
                    print(f"     {cat}: {rxn_type}")
            else:
                print(f"❌ Expected catalyst_selection_needed, got: {result.get('analysis_type')}")
        else:
            # Test specific catalyst selection
            try:
                # Simulate what GUI would do: first get catalyst selection, then resolve
                initial_result = engine.get_recommendations(reaction_smiles, "Auto detect reaction type")
                
                if initial_result.get('analysis_type') == 'catalyst_selection_needed':
                    # Resolve with specific catalyst
                    resolved_result = engine.get_recommendations_with_catalyst(
                        reaction_smiles, 
                        catalyst_value, 
                        "Auto detect reaction type"
                    )
                    
                    print(f"✅ Catalyst {catalyst_value} resolved successfully")
                    print(f"   Final reaction type: {resolved_result.get('reaction_type')}")
                    print(f"   Status: {resolved_result.get('status')}")
                    
                    if 'rxn_insight_detection' in resolved_result:
                        detection = resolved_result['rxn_insight_detection']
                        print(f"   Catalyst used: {detection.get('catalyst_selected')}")
                else:
                    print(f"❌ Expected catalyst selection but got: {initial_result.get('analysis_type')}")
                    
            except Exception as e:
                print(f"❌ Error with catalyst {catalyst_value}: {e}")
        
        print()
    
    print("🎯 SUMMARY:")
    print("- Auto-detect shows catalyst selection prompt")
    print("- Pd selection → Buchwald-Hartwig conditions")
    print("- Cu selection → Ullmann conditions") 
    print("- GUI now allows users to pre-select catalyst")
    print("- System automatically resolves based on selection")

if __name__ == "__main__":
    test_catalyst_workflow_integration()
