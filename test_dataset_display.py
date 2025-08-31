#!/usr/bin/env python3
"""
Test dataset name display in prediction results
"""
import sys
import os
sys.path.insert(0, os.path.dirname(__file__))

from enhanced_recommendation_engine import create_recommendation_engine

def test_dataset_name_display():
    """Test that dataset names are shown in predictions"""
    print("🧪 TESTING DATASET NAME DISPLAY")
    print("=" * 50)
    
    engine = create_recommendation_engine()
    
    # Test Ullmann reaction
    print("\n1. Testing Ullmann reaction:")
    ullmann_smiles = "Brc1ccc(C#N)cc1.Nc1ccccc1>>N#Cc1ccc(Nc2ccccc2)cc1"
    result = engine.get_recommendations_with_catalyst(
        ullmann_smiles, 
        "Cu",  # catalyst_choice
        "Auto detect reaction type"  # reaction_type
    )
    
    # Check if dataset info exists
    dataset_info = result.get('dataset_info', {})
    dataset_name = dataset_info.get('dataset_name', 'Not found')
    
    print(f"   - Reaction Type: {result.get('reaction_type', 'Unknown')}")
    print(f"   - Dataset Name: {dataset_name}")
    print(f"   - Status: {'✅ Found' if dataset_name != 'Not found' else '❌ Missing'}")
    
    # Test Buchwald reaction  
    print("\n2. Testing Buchwald reaction:")
    buchwald_smiles = "Brc1ccccc1.NC1CCCCC1>>c1ccc(N2CCCCC2)cc1"
    result = engine.get_recommendations_with_catalyst(
        buchwald_smiles,
        "Pd",  # catalyst_choice
        "Auto detect reaction type"  # reaction_type
    )
    
    dataset_info = result.get('dataset_info', {})
    dataset_name = dataset_info.get('dataset_name', 'Not found')
    
    print(f"   - Reaction Type: {result.get('reaction_type', 'Unknown')}")
    print(f"   - Dataset Name: {dataset_name}")
    print(f"   - Status: {'✅ Found' if dataset_name != 'Not found' else '❌ Missing'}")
    
    # Test manual reaction type selection
    print("\n3. Testing manual Buchwald selection:")
    result = engine.get_recommendations(
        buchwald_smiles,
        "Cross-Coupling"
    )
    
    dataset_info = result.get('dataset_info', {})
    dataset_name = dataset_info.get('dataset_name', 'Not found')
    
    print(f"   - Reaction Type: {result.get('reaction_type', 'Unknown')}")
    print(f"   - Dataset Name: {dataset_name}")
    print(f"   - Status: {'✅ Found' if dataset_name != 'Not found' else '❌ Missing'}")
    
    print("\n📊 Test complete!")

if __name__ == "__main__":
    test_dataset_name_display()
