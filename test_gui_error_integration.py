#!/usr/bin/env python3

"""Test GUI error handling integration"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from enhanced_recommendation_engine import create_recommendation_engine

def test_gui_error_integration():
    """Test that error responses contain all the fields the GUI expects"""
    print("Testing GUI error integration...")
    
    engine = create_recommendation_engine()
    
    # Test unsupported reaction type
    result = engine.get_recommendations("CC(=O)C>>CC(O)C", "Hydrogenation - Catalytic")
    
    print("Error response structure:")
    for key, value in result.items():
        print(f"  {key}: {value}")
    
    # Check that all GUI-expected fields are present
    expected_fields = [
        'analysis_type', 'error', 'status', 'reaction_type', 'detected_from',
        'ligand_recommendations', 'solvent_recommendations', 'base_recommendations',
        'combined_conditions', 'property_based_alternatives', 'reaction_specific_notes'
    ]
    
    missing_fields = [field for field in expected_fields if field not in result]
    
    if missing_fields:
        print(f"❌ Missing fields that GUI might expect: {missing_fields}")
        return False
    else:
        print("✅ All expected GUI fields present in error response")
    
    # Verify empty recommendations are lists
    list_fields = ['ligand_recommendations', 'solvent_recommendations', 'base_recommendations']
    for field in list_fields:
        if not isinstance(result.get(field), list):
            print(f"❌ {field} should be a list but is {type(result.get(field))}")
            return False
    
    print("✅ All recommendation fields are properly formatted as empty lists")
    return True

if __name__ == "__main__":
    if test_gui_error_integration():
        print("\n🎉 GUI error integration test passed")
    else:
        print("\n❌ GUI error integration test failed")
