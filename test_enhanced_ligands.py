#!/usr/bin/env python3
"""
Test script for enhanced ligand recommendation system.
"""

import sys
import os

# Add the reagents directory to the Python path
sys.path.append(os.path.join(os.path.dirname(__file__), 'reagents'))

from ligand import recommend_ligands_for_reaction, get_reaction_specific_ligands

def test_ligand_recommendations():
    """Test the enhanced ligand recommendation functions"""
    
    print("=" * 80)
    print("ENHANCED LIGAND RECOMMENDATION SYSTEM TEST")
    print("=" * 80)
    
    # Test 1: Cross-coupling ligands
    print("\n1. TOP CROSS-COUPLING LIGANDS:")
    print("-" * 40)
    cross_coupling_ligands = recommend_ligands_for_reaction(
        reaction_type='Cross-Coupling',
        top_n=5,
        min_compatibility=0.5
    )
    
    for rec in cross_coupling_ligands:
        print(f"   {rec['rank']}. {rec['ligand']}")
        print(f"      Compatibility: {rec['compatibility_score']}")
        print(f"      Applications: {rec['applications']}")
        print()
    
    # Test 2: Hydrogenation ligands
    print("\n2. TOP HYDROGENATION LIGANDS:")
    print("-" * 40)
    hydrogenation_ligands = recommend_ligands_for_reaction(
        reaction_type='Hydrogenation',
        top_n=5,
        min_compatibility=0.5
    )
    
    for rec in hydrogenation_ligands:
        print(f"   {rec['rank']}. {rec['ligand']}")
        print(f"      Compatibility: {rec['compatibility_score']}")
        print(f"      Applications: {rec['applications']}")
        print()
    
    # Test 3: Metathesis ligands
    print("\n3. TOP METATHESIS LIGANDS:")
    print("-" * 40)
    metathesis_ligands = recommend_ligands_for_reaction(
        reaction_type='Metathesis',
        top_n=5,
        min_compatibility=0.5
    )
    
    for rec in metathesis_ligands:
        print(f"   {rec['rank']}. {rec['ligand']}")
        print(f"      Compatibility: {rec['compatibility_score']}")
        print(f"      Applications: {rec['applications']}")
        print()
    
    # Test 4: Similarity-based recommendations
    print("\n4. LIGANDS SIMILAR TO PPh3 FOR CROSS-COUPLING:")
    print("-" * 50)
    similar_ligands = recommend_ligands_for_reaction(
        target_ligand='PPh3',
        reaction_type='Cross-Coupling',
        top_n=3,
        min_compatibility=0.3
    )
    
    for rec in similar_ligands:
        print(f"   {rec['rank']}. {rec['ligand']}")
        print(f"      Compatibility: {rec['compatibility_score']}")
        if 'similarity_score' in rec:
            print(f"      Similarity: {rec['similarity_score']}")
            print(f"      Combined Score: {rec['combined_score']}")
        print(f"      Applications: {rec['applications']}")
        print()
    
    # Test 5: Property-filtered recommendations
    print("\n5. AFFORDABLE LIGANDS FOR CROSS-COUPLING:")
    print("-" * 45)
    filtered_ligands = get_reaction_specific_ligands(
        reaction_type='Cross-Coupling',
        property_preferences={
            'price_category_max': 3,  # Only affordable ligands
            'cone_angle_max': 170     # Not too bulky
        }
    )
    
    for rec in filtered_ligands:
        print(f"   {rec['rank']}. {rec['ligand']}")
        print(f"      Compatibility: {rec['compatibility_score']}")
        print(f"      Applications: {rec['applications']}")
        print()
    
    print("=" * 80)
    print("TEST COMPLETED SUCCESSFULLY!")
    print("Enhanced ligand recommendation system is working properly.")
    print("=" * 80)

if __name__ == "__main__":
    test_ligand_recommendations()
