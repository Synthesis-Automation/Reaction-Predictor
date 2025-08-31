#!/usr/bin/env python3

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/..")

try:
    from enhanced_recommendation_engine import create_recommendation_engine
    
    print("Testing reaction_suitability field fix...")
    e = create_recommendation_engine()
    r = e.get_recommendations('O=C(O)c1ccccc1.Nc1ccccc1>>O=C(Nc1ccccc1)c1ccccc1','Amide Formation - Acid + Amine')
    
    # Check ligands
    ligands = r.get('ligand_recommendations', [])
    if ligands:
        first_ligand = ligands[0]
        has_suitability = 'reaction_suitability' in first_ligand
        if has_suitability:
            print(f"✅ LIGANDS: reaction_suitability field present")
            print(f"  Sample: {first_ligand['ligand']}")
            print(f"  Suitability: {first_ligand['reaction_suitability']}")
        else:
            print("❌ LIGANDS: Missing reaction_suitability field")
            print(f"  Available fields: {list(first_ligand.keys())}")
    
    # Check solvents  
    solvents = r.get('solvent_recommendations', [])
    if solvents:
        first_solvent = solvents[0]
        has_suitability = 'reaction_suitability' in first_solvent
        if has_suitability:
            print(f"✅ SOLVENTS: reaction_suitability field present")
            print(f"  Sample: {first_solvent['solvent']}")
            print(f"  Suitability: {first_solvent['reaction_suitability']}")
        else:
            print("❌ SOLVENTS: Missing reaction_suitability field")
            print(f"  Available fields: {list(first_solvent.keys())}")
            
    # Test all required fields for GUI compatibility
    required_ligand_fields = ['ligand', 'compatibility_score', 'applications', 'reaction_suitability']
    required_solvent_fields = ['solvent', 'abbreviation', 'compatibility_score', 'reaction_suitability']
    
    ligand_missing = []
    solvent_missing = []
    
    if ligands:
        for field in required_ligand_fields:
            if field not in ligands[0]:
                ligand_missing.append(field)
                
    if solvents:
        for field in required_solvent_fields:
            if field not in solvents[0]:
                solvent_missing.append(field)
    
    if not ligand_missing and not solvent_missing:
        print("\n🎉 SUCCESS: All required fields present for GUI compatibility!")
    else:
        print(f"\n❌ MISSING FIELDS:")
        if ligand_missing:
            print(f"  Ligands missing: {ligand_missing}")
        if solvent_missing:
            print(f"  Solvents missing: {solvent_missing}")
        
except Exception as ex:
    print(f"❌ ERROR: {ex}")
    import traceback
    traceback.print_exc()
