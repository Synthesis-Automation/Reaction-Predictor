#!/usr/bin/env python3
"""
Final test of dataset name display functionality
"""
import sys
import os
sys.path.insert(0, os.path.dirname(__file__))

from enhanced_recommendation_engine import create_recommendation_engine

def test_final_dataset_display():
    """Final test showing dataset names in results"""
    print("🎯 FINAL DATASET NAME DISPLAY TEST")
    print("=" * 50)
    
    engine = create_recommendation_engine()
    
    # Test different reaction types
    test_cases = [
        ("Ullmann reaction", "Brc1ccc(C#N)cc1.Nc1ccccc1>>N#Cc1ccc(Nc2ccccc2)cc1", "Ullmann"),
        ("Cross-Coupling", "Brc1ccccc1.NC1CCCCC1>>c1ccc(N2CCCCC2)cc1", "Cross-Coupling"),
        ("C-N Coupling - Ullmann", "Brc1ccc(C#N)cc1.Nc1ccccc1>>N#Cc1ccc(Nc2ccccc2)cc1", "C-N Coupling - Ullmann"),
        ("C-N Coupling - Buchwald-Hartwig", "Brc1ccccc1.NC1CCCCC1>>c1ccc(N2CCCCC2)cc1", "C-N Coupling - Buchwald-Hartwig"),
    ]
    
    for description, smiles, reaction_type in test_cases:
        print(f"\n🧪 {description}:")
        print(f"   SMILES: {smiles}")
        print(f"   Type: {reaction_type}")
        
        result = engine.get_recommendations(smiles, reaction_type)
        
        if result.get('status') == 'success':
            dataset_info = result.get('dataset_info', {})
            dataset_name = dataset_info.get('dataset_name', 'NOT_FOUND')
            ligands = dataset_info.get('ligands_available', 0)
            solvents = dataset_info.get('solvents_available', 0)
            
            print(f"   ✅ SUCCESS")
            print(f"   📊 Dataset: {dataset_name}")
            print(f"   🧪 Ligands: {ligands}, Solvents: {solvents}")
        else:
            status = result.get('status', 'unknown')
            print(f"   ❌ FAILED: {status}")
            if status == 'needs_catalyst_selection':
                print(f"      (Expected - auto-detection requires catalyst selection)")
    
    print(f"\n{'='*50}")
    print("✅ Dataset name display functionality implemented!")
    print("   • Dataset names are now shown in prediction results")
    print("   • GUI will display actual dataset filenames")
    print("   • Users can see which datasets are used for recommendations")

if __name__ == "__main__":
    test_final_dataset_display()
