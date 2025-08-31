#!/usr/bin/env python3

import sys
import os
sys.path.insert(0, os.path.abspath('.'))

from enhanced_recommendation_engine import create_recommendation_engine

def test_full_recommendations():
    """Test full recommendation system"""
    print("Testing Full Recommendation System for Amide Formation")
    print("=" * 55)
    
    engine = create_recommendation_engine()
    
    # Test with a simple amide formation reaction
    test_reaction = "O=C(O)c1ccccc1.Nc1ccccc1>>O=C(Nc1ccccc1)c1ccccc1"
    reaction_type = "Amide Formation - Acid + Amine"
    
    print(f"\nTest reaction: {test_reaction}")
    print(f"Reaction type: {reaction_type}")
    
    # Get recommendations
    recommendations = engine.get_recommendations(test_reaction, reaction_type)
    
    print(f"\nDetected reaction type: {recommendations.get('reaction_type', 'Unknown')}")
    
    # Check what's in the recommendations
    print(f"\nKeys in recommendations: {list(recommendations.keys())}")
    
    # Show ligand recommendations
    ligands = recommendations.get('ligand_recommendations', [])
    print(f"\nLigand recommendations count: {len(ligands)}")
    if ligands:
        print("Ligand details:")
        for i, lig in enumerate(ligands, 1):
            print(f"  {i}. {lig}")
    
    # Show solvent recommendations  
    solvents = recommendations.get('solvent_recommendations', [])
    print(f"\nSolvent recommendations count: {len(solvents)}")
    if solvents:
        print("Solvent details:")
        for i, sol in enumerate(solvents, 1):
            print(f"  {i}. {sol}")
    
    # Show base recommendations
    bases = recommendations.get('base_recommendations', [])
    print(f"\nBase recommendations count: {len(bases)}")
    if bases:
        print("Base details:")
        for i, base in enumerate(bases, 1):
            print(f"  {i}. {base}")

if __name__ == "__main__":
    test_full_recommendations()
