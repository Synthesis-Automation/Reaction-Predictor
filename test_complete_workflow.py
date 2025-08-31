#!/usr/bin/env python3
"""Test the complete catalyst selection workflow in the GUI format"""

from enhanced_recommendation_engine import create_recommendation_engine

def test_gui_workflow():
    """Test the complete workflow as it would happen in GUI"""
    
    reaction_smiles = 'Brc1ccc(C#N)cc1.Nc1ccccc1>>N#Cc1ccc(Nc2ccccc2)cc1'
    reaction_type = "Auto detect reaction type"  # This is what GUI sends
    
    print("🔬 TESTING COMPLETE GUI WORKFLOW")
    print("=" * 50)
    print(f"Reaction SMILES: {reaction_smiles}")
    print(f"Reaction Type: {reaction_type}")
    print()
    
    engine = create_recommendation_engine()
    
    # Step 1: Initial prediction (like GUI worker does)
    print("Step 1: Initial prediction request")
    print("-" * 30)
    
    result = engine.get_recommendations(reaction_smiles, reaction_type)
    
    print(f"Analysis Type: {result.get('analysis_type')}")
    print(f"Status: {result.get('status')}")
    
    if result.get('analysis_type') == 'catalyst_selection_needed':
        print("✅ Catalyst selection correctly detected!")
        print()
        print("GUI should show:")
        print(f"• Message: {result.get('message')}")
        print("• Catalyst options:")
        for cat, rxn_type in result.get('catalyst_options', {}).items():
            print(f"  - {cat}: {rxn_type}")
        print()
        
        # Step 2: User selects Pd (Buchwald-Hartwig)
        print("Step 2: User selects Pd catalyst")
        print("-" * 30)
        
        pd_result = engine.get_recommendations_with_catalyst(reaction_smiles, 'Pd', reaction_type)
        print(f"Final reaction type: {pd_result.get('reaction_type')}")
        print(f"Status: {pd_result.get('status')}")
        
        if 'rxn_insight_detection' in pd_result:
            detection = pd_result['rxn_insight_detection']
            print(f"Catalyst selected: {detection.get('catalyst_selected')}")
            print(f"rxn-insight detected: {detection.get('rxn_insight_detected')}")
        
        print()
        
        # Step 3: User selects Cu (Ullmann)
        print("Step 3: User selects Cu catalyst")
        print("-" * 30)
        
        cu_result = engine.get_recommendations_with_catalyst(reaction_smiles, 'Cu', reaction_type)
        print(f"Final reaction type: {cu_result.get('reaction_type')}")
        print(f"Status: {cu_result.get('status')}")
        
        if 'rxn_insight_detection' in cu_result:
            detection = cu_result['rxn_insight_detection']
            print(f"Catalyst selected: {detection.get('catalyst_selected')}")
            print(f"rxn-insight detected: {detection.get('rxn_insight_detected')}")
            
    else:
        print("❌ Expected catalyst_selection_needed but got:")
        print(f"  Analysis type: {result.get('analysis_type')}")
        print(f"  Reaction type: {result.get('reaction_type')}")
        print(f"  Error: {result.get('error')}")

if __name__ == "__main__":
    test_gui_workflow()
