#!/usr/bin/env python3

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/..")

try:
    from enhanced_recommendation_engine import create_recommendation_engine
    
    print("Testing abbreviation fix...")
    e = create_recommendation_engine()
    r = e.get_recommendations('O=C(O)c1ccccc1.Nc1ccccc1>>O=C(Nc1ccccc1)c1ccccc1','Amide Formation - Acid + Amine')
    
    # Check if we have solvents with abbreviations
    solvents = r.get('solvent_recommendations', [])
    if solvents:
        first_solvent = solvents[0]
        has_abbrev = 'abbreviation' in first_solvent
        if has_abbrev:
            print(f"✅ SUCCESS: Abbreviation error FIXED!")
            print(f"Sample: {first_solvent['solvent']} ({first_solvent['abbreviation']})")
        else:
            print("❌ FAILED: Still missing abbreviation field")
    else:
        print("❌ FAILED: No solvents found")
        
except Exception as ex:
    print(f"❌ ERROR: {ex}")
    import traceback
    traceback.print_exc()
