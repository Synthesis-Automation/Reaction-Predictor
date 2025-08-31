#!/usr/bin/env python3

import sys
import os
sys.path.insert(0, os.path.abspath('.'))

from enhanced_recommendation_engine import create_recommendation_engine

def test_abbreviations():
    """Test if abbreviations are working correctly"""
    print("Testing Abbreviations")
    print("=" * 30)
    
    engine = create_recommendation_engine()
    recs = engine.get_recommendations('O=C(O)c1ccccc1.Nc1ccccc1>>O=C(Nc1ccccc1)c1ccccc1','Amide Formation - Acid + Amine')
    
    solvents = recs.get('solvent_recommendations', [])
    print("Solvents with abbreviations:")
    for i, s in enumerate(solvents[:3], 1):
        solvent_name = s.get('solvent', 'Unknown')
        abbrev = s.get('abbreviation', 'NO_ABBREV') 
        score = s.get('compatibility_score', 0)
        print(f"  {i}. {solvent_name} ({abbrev}) - Score: {score:.3f}")
    
    # Test if any solvent is missing abbreviation
    missing_abbrev = [s for s in solvents if 'abbreviation' not in s]
    if missing_abbrev:
        print(f"\n❌ WARNING: {len(missing_abbrev)} solvents missing abbreviation field!")
        for s in missing_abbrev:
            print(f"  - {s}")
    else:
        print(f"\n✅ All {len(solvents)} solvents have abbreviation field")

if __name__ == "__main__":
    test_abbreviations()
