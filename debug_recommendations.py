#!/usr/bin/env python3
"""
Debug recommendation engine results
"""
import sys
import os
sys.path.insert(0, os.path.dirname(__file__))

from enhanced_recommendation_engine import create_recommendation_engine

def debug_recommendations():
    """Debug what the recommendation engine returns"""
    print("🐛 DEBUGGING RECOMMENDATION ENGINE RESULTS")
    print("=" * 50)
    
    engine = create_recommendation_engine()
    
    # Test simple manual selection first
    print("\n1. Testing manual 'Cross-Coupling' selection:")
    smiles = "Brc1ccccc1.NC1CCCCC1>>c1ccc(N2CCCCC2)cc1"
    result = engine.get_recommendations(smiles, "Cross-Coupling")
    
    print(f"   Full result keys: {list(result.keys())}")
    print(f"   Reaction Type: {result.get('reaction_type', 'NOT_FOUND')}")
    print(f"   Status: {result.get('status', 'NOT_FOUND')}")
    
    if 'dataset_info' in result:
        dataset_info = result['dataset_info']
        print(f"   Dataset Info: {dataset_info}")
        print(f"   Dataset Name: {dataset_info.get('dataset_name', 'NOT_FOUND')}")
    else:
        print("   ❌ No dataset_info in result")
    
    # Test Ullmann manual selection
    print("\n2. Testing manual 'Ullmann' selection:")
    ullmann_smiles = "Brc1ccc(C#N)cc1.Nc1ccccc1>>N#Cc1ccc(Nc2ccccc2)cc1"
    result = engine.get_recommendations(ullmann_smiles, "Ullmann")
    
    print(f"   Full result keys: {list(result.keys())}")
    print(f"   Reaction Type: {result.get('reaction_type', 'NOT_FOUND')}")
    print(f"   Status: {result.get('status', 'NOT_FOUND')}")
    
    if 'dataset_info' in result:
        dataset_info = result['dataset_info']
        print(f"   Dataset Info: {dataset_info}")
        print(f"   Dataset Name: {dataset_info.get('dataset_name', 'NOT_FOUND')}")
    else:
        print("   ❌ No dataset_info in result")

    # Test C-N Coupling manual selection
    print("\n3. Testing manual 'C-N Coupling - Ullmann' selection:")
    result = engine.get_recommendations(ullmann_smiles, "C-N Coupling - Ullmann")
    
    print(f"   Full result keys: {list(result.keys())}")
    print(f"   Reaction Type: {result.get('reaction_type', 'NOT_FOUND')}")
    print(f"   Status: {result.get('status', 'NOT_FOUND')}")
    
    if 'dataset_info' in result:
        dataset_info = result['dataset_info']
        print(f"   Dataset Info: {dataset_info}")
        print(f"   Dataset Name: {dataset_info.get('dataset_name', 'NOT_FOUND')}")
    else:
        print("   ❌ No dataset_info in result")

if __name__ == "__main__":
    debug_recommendations()
