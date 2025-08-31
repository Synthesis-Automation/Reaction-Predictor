#!/usr/bin/env python3
"""Simple demonstration of the catalyst selection workflow"""

from enhanced_recommendation_engine import create_recommendation_engine

def demo_catalyst_workflow():
    """Simple demo of how the catalyst selection works"""
    
    reaction = 'Brc1ccc(C(F)(F)F)cc1.Nc1ccccc1>>FC(F)(F)c1ccc(Nc2ccccc2)cc1'
    
    print("🧪 CATALYST SELECTION DEMO")
    print("Reaction:", reaction)
    print()
    
    engine = create_recommendation_engine()
    
    # Initial detection
    result = engine.get_recommendations(reaction, "Auto detect reaction type")
    
    if result.get('status') == 'needs_catalyst_selection':
        print("✅ rxn-insight detected: N-arylation")
        print("❓ No catalyst found in SMILES - need user selection")
        print()
        print("Available options:")
        for cat, rxn_type in result.get('catalyst_options', {}).items():
            print(f"  • {cat} → {rxn_type}")
        print()
        
        # Test Buchwald-Hartwig (Pd)
        print("User selects: Pd (Buchwald-Hartwig)")
        pd_result = engine.get_recommendations_with_catalyst(reaction, 'Pd')
        
        if pd_result.get('status') == 'success':
            print(f"✅ Final type: {pd_result.get('reaction_type')}")
            if 'rxn_insight_detection' in pd_result:
                det = pd_result['rxn_insight_detection']
                print(f"   Catalyst: {det.get('catalyst_selected')}")
                print(f"   rxn-insight: {det.get('rxn_insight_detected')}")
        else:
            print(f"❌ Error: {pd_result.get('error')}")
            
        print()
        
        # Test Ullmann (Cu)  
        print("User selects: Cu (Ullmann)")
        cu_result = engine.get_recommendations_with_catalyst(reaction, 'Cu')
        
        if cu_result.get('status') == 'success':
            print(f"✅ Final type: {cu_result.get('reaction_type')}")
            if 'rxn_insight_detection' in cu_result:
                det = cu_result['rxn_insight_detection']
                print(f"   Catalyst: {det.get('catalyst_selected')}")
                print(f"   rxn-insight: {det.get('rxn_insight_detected')}")
        else:
            print(f"❌ Error: {cu_result.get('error')}")
            
    else:
        print("❌ Expected catalyst selection, got:", result.get('status'))

if __name__ == "__main__":
    demo_catalyst_workflow()
