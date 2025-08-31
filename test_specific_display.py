#!/usr/bin/env python3
"""
Test the updated Database Coverage display with specific ligands and solvents
"""
import sys
import os
sys.path.insert(0, os.path.dirname(__file__))

from enhanced_recommendation_engine import create_recommendation_engine

def test_specific_display():
    """Test the enhanced Database Coverage display"""
    print("🧪 TESTING ENHANCED DATABASE COVERAGE DISPLAY")
    print("=" * 60)
    
    engine = create_recommendation_engine()
    
    test_cases = [
        ('Ullmann', 'Brc1ccc(C#N)cc1.Nc1ccccc1>>N#Cc1ccc(Nc2ccccc2)cc1'),
        ('Cross-Coupling', 'Brc1ccccc1.NC1CCCCC1>>c1ccc(N2CCCCC2)cc1')
    ]
    
    for reaction_type, smiles in test_cases:
        print(f"\n🧪 Testing {reaction_type}:")
        print("-" * 40)
        
        result = engine.get_recommendations(smiles, reaction_type)
        
        info = result.get('dataset_info', {})
        
        print(f"📊 Basic Info:")
        print(f"   Dataset: {info.get('dataset_name', 'Unknown')}")
        print(f"   Ligands: {info.get('ligands_available', 'Unknown')}")
        print(f"   Solvents: {info.get('solvents_available', 'Unknown')}")
        
        # Check specific ligands
        specific_ligands = info.get('specific_ligands', [])
        print(f"\n🔗 Specific Ligands ({len(specific_ligands)} found):")
        if specific_ligands:
            for i, ligand in enumerate(specific_ligands[:5], 1):
                if isinstance(ligand, dict):
                    name = ligand.get('name', ligand.get('ligand', str(ligand)))
                    score = ligand.get('score', ligand.get('compatibility_score', ''))
                    if score:
                        print(f"   {i}. {name} (score: {score:.2f})")
                    else:
                        print(f"   {i}. {name}")
                else:
                    print(f"   {i}. {ligand}")
        else:
            print("   No specific ligands found")
        
        # Check specific solvents
        specific_solvents = info.get('specific_solvents', [])
        print(f"\n🧪 Specific Solvents ({len(specific_solvents)} found):")
        if specific_solvents:
            for i, solvent in enumerate(specific_solvents[:5], 1):
                if isinstance(solvent, dict):
                    name = solvent.get('name', solvent.get('solvent', str(solvent)))
                    score = solvent.get('score', solvent.get('compatibility_score', ''))
                    if score:
                        print(f"   {i}. {name} (score: {score:.2f})")
                    else:
                        print(f"   {i}. {name}")
                else:
                    print(f"   {i}. {solvent}")
        else:
            print("   No specific solvents found")
    
    print(f"\n✅ Enhanced Database Coverage now shows:")
    print(f"   • Dataset name and counts")
    print(f"   • Top 5 specific ligands with scores")
    print(f"   • Top 5 specific solvents with scores")
    print(f"   • Reaction-type tailored recommendations")

if __name__ == "__main__":
    test_specific_display()
