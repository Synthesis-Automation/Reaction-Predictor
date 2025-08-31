#!/usr/bin/env python3

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/..")

def test_gui_formatting():
    """Test the exact formatting logic that was causing the GUI error"""
    try:
        from enhanced_recommendation_engine import create_recommendation_engine
        
        print("Testing GUI formatting compatibility...")
        e = create_recommendation_engine()
        recommendations = e.get_recommendations('O=C(O)c1ccccc1.Nc1ccccc1>>O=C(Nc1ccccc1)c1ccccc1','Amide Formation - Acid + Amine')
        
        # Simulate the GUI formatting logic that was failing
        ligand_recs = recommendations.get('ligand_recommendations', [])
        
        if ligand_recs:
            print("🧬 TOP LIGAND OPTIONS:")
            print("-" * 40)
            
            for i, lig in enumerate(ligand_recs[:5], 1):
                # This is the exact line that was causing the KeyError
                formatted_text = f"""  {i}. {lig['ligand']} (Score: {lig['compatibility_score']})
     • Applications: {lig['applications']}
     • Suitability: {lig['reaction_suitability']}
"""
                print(formatted_text)
            print("✅ SUCCESS: Ligand formatting works!")
        else:
            print("❌ No ligands found")
            
        # Test solvent formatting too
        solvent_recs = recommendations.get('solvent_recommendations', [])
        if solvent_recs:
            print("\n🧪 TOP SOLVENT OPTIONS:")
            print("-" * 40)
            
            for i, sol in enumerate(solvent_recs[:3], 1):
                # Test solvent formatting
                formatted_text = f"""  {i}. {sol['solvent']} ({sol['abbreviation']}) (Score: {sol['compatibility_score']})
     • Suitability: {sol['reaction_suitability']}
"""
                print(formatted_text)
            print("✅ SUCCESS: Solvent formatting works!")
        else:
            print("❌ No solvents found")
            
        print("\n🎉 OVERALL SUCCESS: GUI formatting test passed!")
        return True
        
    except KeyError as e:
        print(f"❌ KeyError still exists: {e}")
        return False
    except Exception as ex:
        print(f"❌ OTHER ERROR: {ex}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = test_gui_formatting()
    sys.exit(0 if success else 1)
