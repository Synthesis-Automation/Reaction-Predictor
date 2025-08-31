#!/usr/bin/env python3
"""
Test dataset-specific ligand and solvent counts
"""
import sys
import os
sys.path.insert(0, os.path.dirname(__file__))

from enhanced_recommendation_engine import create_recommendation_engine

def test_dataset_specific_counts():
    """Test that different datasets show different ligand/solvent counts"""
    print("🧪 TESTING DATASET-SPECIFIC COUNTS")
    print("=" * 60)
    
    engine = create_recommendation_engine()
    
    test_cases = [
        ('Ullmann', 'Brc1ccc(C#N)cc1.Nc1ccccc1>>N#Cc1ccc(Nc2ccccc2)cc1'),
        ('Cross-Coupling', 'Brc1ccccc1.NC1CCCCC1>>c1ccc(N2CCCCC2)cc1'),
        ('C-N Coupling - Buchwald-Hartwig', 'Brc1ccccc1.NC1CCCCC1>>c1ccc(N2CCCCC2)cc1')
    ]
    
    for reaction_type, smiles in test_cases:
        print(f"\n🧪 Testing {reaction_type}:")
        result = engine.get_recommendations(smiles, reaction_type)
        
        info = result.get('dataset_info', {})
        dataset_name = info.get('dataset_name', 'Unknown')
        ligands = info.get('ligands_available', 'Unknown')
        solvents = info.get('solvents_available', 'Unknown')
        
        print(f"   Dataset: {dataset_name}")
        print(f"   Ligands: {ligands}")
        print(f"   Solvents: {solvents}")
        print(f"   Status: {result.get('status', 'Unknown')}")
        
        # Check if counts are reaction-specific (should be 10 each) or global (123/72)
        if ligands == 10 and solvents == 10:
            print("   ✅ Using reaction-specific counts!")
        elif ligands == 123 and solvents == 72:
            print("   ⚠️  Still using global counts")
        else:
            print(f"   ❓ Unexpected counts: {ligands}/{solvents}")
    
    print(f"\n📊 Expected Result:")
    print(f"   Each dataset should show ~10 ligands and ~10 solvents")
    print(f"   (reaction-specific counts, not global database size)")

if __name__ == "__main__":
    test_dataset_specific_counts()
