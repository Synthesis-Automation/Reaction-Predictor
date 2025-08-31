#!/usr/bin/env python3

import sys
import os
sys.path.insert(0, os.path.abspath('.'))

from enhanced_recommendation_engine import create_recommendation_engine

def test_full_pipeline():
    """Test full recommendation pipeline for abbreviation issues"""
    print("Testing Full Amide Formation Recommendation Pipeline")
    print("=" * 55)
    
    engine = create_recommendation_engine()
    reaction = 'O=C(O)c1ccccc1.Nc1ccccc1>>O=C(Nc1ccccc1)c1ccccc1'
    reaction_type = 'Amide Formation - Acid + Amine'
    
    try:
        # Get full recommendations
        recs = engine.get_recommendations(reaction, reaction_type)
        
        print("✅ Basic recommendations generated successfully")
        
        # Check ligand recommendations
        ligands = recs.get('ligand_recommendations', [])
        print(f"✅ Ligand recommendations: {len(ligands)} items")
        
        # Check solvent recommendations for abbreviations
        solvents = recs.get('solvent_recommendations', [])
        print(f"✅ Solvent recommendations: {len(solvents)} items")
        
        abbreviation_issues = []
        for i, solvent in enumerate(solvents):
            if 'abbreviation' not in solvent:
                abbreviation_issues.append(f"Solvent {i+1}: {solvent.get('solvent', 'Unknown')} missing abbreviation")
            elif solvent['abbreviation'] is None:
                abbreviation_issues.append(f"Solvent {i+1}: {solvent.get('solvent', 'Unknown')} has None abbreviation")
        
        if abbreviation_issues:
            print("❌ Abbreviation issues found:")
            for issue in abbreviation_issues:
                print(f"  - {issue}")
        else:
            print("✅ All solvents have valid abbreviations")
        
        # Check base recommendations
        bases = recs.get('base_recommendations', [])
        print(f"✅ Base recommendations: {len(bases)} items")
        
        # Check combined conditions (this uses solvent abbreviations)
        combined = recs.get('combined_conditions', [])
        print(f"✅ Combined conditions: {len(combined)} items")
        
        # Test if any combined condition has abbreviation issues
        combined_issues = []
        for i, combo in enumerate(combined):
            if 'solvent_abbreviation' not in combo:
                combined_issues.append(f"Combo {i+1}: missing solvent_abbreviation")
            elif combo['solvent_abbreviation'] is None:
                combined_issues.append(f"Combo {i+1}: None solvent_abbreviation")
        
        if combined_issues:
            print("❌ Combined condition abbreviation issues:")
            for issue in combined_issues:
                print(f"  - {issue}")
        else:
            print("✅ All combined conditions have valid solvent abbreviations")
        
        print(f"\n🎯 SUMMARY:")
        print(f"  • Ligands: {len(ligands)}")
        print(f"  • Solvents: {len(solvents)}")  
        print(f"  • Bases: {len(bases)}")
        print(f"  • Combined: {len(combined)}")
        print(f"  • Abbreviation errors: {len(abbreviation_issues + combined_issues)}")
        
        if len(abbreviation_issues + combined_issues) == 0:
            print("\n🎉 SUCCESS: No abbreviation errors found!")
            return True
        else:
            print(f"\n❌ FAILED: {len(abbreviation_issues + combined_issues)} abbreviation errors found!")
            return False
            
    except Exception as e:
        print(f"❌ ERROR: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = test_full_pipeline()
    sys.exit(0 if success else 1)
