#!/usr/bin/env python3

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/..")

def test_all_required_fields():
    """Test that all GUI-required fields are present"""
    try:
        from enhanced_recommendation_engine import create_recommendation_engine
        
        print("Testing all required fields for GUI compatibility...")
        e = create_recommendation_engine()
        recommendations = e.get_recommendations('O=C(O)c1ccccc1.Nc1ccccc1>>O=C(Nc1ccccc1)c1ccccc1','Amide Formation - Acid + Amine')
        
        # Define required fields based on GUI usage
        required_ligand_fields = ['ligand', 'compatibility_score', 'applications', 'reaction_suitability']
        required_solvent_fields = ['solvent', 'abbreviation', 'compatibility_score', 'applications', 'reaction_suitability']
        
        # Check ligands
        ligands = recommendations.get('ligand_recommendations', [])
        ligand_issues = []
        
        if ligands:
            for field in required_ligand_fields:
                if field not in ligands[0]:
                    ligand_issues.append(field)
            
            if not ligand_issues:
                print(f"✅ LIGANDS: All {len(required_ligand_fields)} required fields present")
                print(f"   Sample: {ligands[0]['ligand']}")
            else:
                print(f"❌ LIGANDS: Missing fields: {ligand_issues}")
        else:
            print("❌ No ligands found")
            
        # Check solvents
        solvents = recommendations.get('solvent_recommendations', [])
        solvent_issues = []
        
        if solvents:
            for field in required_solvent_fields:
                if field not in solvents[0]:
                    solvent_issues.append(field)
            
            if not solvent_issues:
                print(f"✅ SOLVENTS: All {len(required_solvent_fields)} required fields present")
                print(f"   Sample: {solvents[0]['solvent']} ({solvents[0]['abbreviation']})")
            else:
                print(f"❌ SOLVENTS: Missing fields: {solvent_issues}")
        else:
            print("❌ No solvents found")
            
        # Overall result
        total_issues = len(ligand_issues) + len(solvent_issues)
        if total_issues == 0:
            print(f"\n🎉 SUCCESS: All GUI compatibility issues resolved!")
            print(f"   • Ligands: {len(ligands)} items with all required fields")
            print(f"   • Solvents: {len(solvents)} items with all required fields")
            return True
        else:
            print(f"\n❌ FAILED: {total_issues} field issues remain")
            return False
            
    except Exception as ex:
        print(f"❌ ERROR: {ex}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = test_all_required_fields()
    sys.exit(0 if success else 1)
