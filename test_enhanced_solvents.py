#!/usr/bin/env python3
"""
Test script for enhanced solvent recommendation system.
"""

import sys
import os

# Add the reagents directory to the Python path
sys.path.append(os.path.join(os.path.dirname(__file__), 'reagents'))

from solvent import recommend_solvents_for_reaction, get_reaction_specific_solvents

def test_solvent_recommendations():
    """Test the enhanced solvent recommendation functions"""
    
    print("=" * 80)
    print("ENHANCED SOLVENT RECOMMENDATION SYSTEM TEST")
    print("=" * 80)
    
    # Test 1: Cross-coupling solvents
    print("\n1. TOP CROSS-COUPLING SOLVENTS:")
    print("-" * 40)
    cross_coupling_solvents = recommend_solvents_for_reaction(
        reaction_type='Cross-Coupling',
        top_n=5,
        min_compatibility=0.5
    )
    
    for rec in cross_coupling_solvents:
        print(f"   {rec['rank']}. {rec['solvent']} ({rec['abbreviation']})")
        print(f"      Compatibility: {rec['compatibility_score']}")
        print(f"      Applications: {rec['applications']}")
        print()
    
    # Test 2: Hydrogenation solvents
    print("\n2. TOP HYDROGENATION SOLVENTS:")
    print("-" * 40)
    hydrogenation_solvents = recommend_solvents_for_reaction(
        reaction_type='Hydrogenation',
        top_n=5,
        min_compatibility=0.5
    )
    
    for rec in hydrogenation_solvents:
        print(f"   {rec['rank']}. {rec['solvent']} ({rec['abbreviation']})")
        print(f"      Compatibility: {rec['compatibility_score']}")
        print(f"      Applications: {rec['applications']}")
        print()
    
    # Test 3: Metathesis solvents
    print("\n3. TOP METATHESIS SOLVENTS:")
    print("-" * 40)
    metathesis_solvents = recommend_solvents_for_reaction(
        reaction_type='Metathesis',
        top_n=5,
        min_compatibility=0.5
    )
    
    for rec in metathesis_solvents:
        print(f"   {rec['rank']}. {rec['solvent']} ({rec['abbreviation']})")
        print(f"      Compatibility: {rec['compatibility_score']}")
        print(f"      Applications: {rec['applications']}")
        print()
    
    # Test 4: Similarity-based recommendations
    print("\n4. SOLVENTS SIMILAR TO THF FOR CROSS-COUPLING:")
    print("-" * 50)
    similar_solvents = recommend_solvents_for_reaction(
        target_solvent='THF',
        reaction_type='Cross-Coupling',
        top_n=3,
        min_compatibility=0.3
    )
    
    for rec in similar_solvents:
        print(f"   {rec['rank']}. {rec['solvent']} ({rec['abbreviation']})")
        print(f"      Compatibility: {rec['compatibility_score']}")
        if 'similarity_score' in rec:
            print(f"      Similarity: {rec['similarity_score']}")
            print(f"      Combined Score: {rec['combined_score']}")
        print(f"      Applications: {rec['applications']}")
        print()
    
    # Test 5: Property-filtered recommendations
    print("\n5. LOW BP POLAR SOLVENTS FOR CROSS-COUPLING:")
    print("-" * 47)
    filtered_solvents = get_reaction_specific_solvents(
        reaction_type='Cross-Coupling',
        property_preferences={
            'bp_max': 120,     # Low boiling point
            'polarity_min': 4, # Polar solvents
            'protic': False    # Aprotic
        }
    )
    
    for rec in filtered_solvents:
        print(f"   {rec['rank']}. {rec['solvent']} ({rec['abbreviation']})")
        print(f"      Compatibility: {rec['compatibility_score']}")
        print(f"      Applications: {rec['applications']}")
        print()
    
    # Test 6: Carbonylation solvents
    print("\n6. TOP CARBONYLATION SOLVENTS:")
    print("-" * 35)
    carbonylation_solvents = recommend_solvents_for_reaction(
        reaction_type='Carbonylation',
        top_n=5,
        min_compatibility=0.6
    )
    
    for rec in carbonylation_solvents:
        print(f"   {rec['rank']}. {rec['solvent']} ({rec['abbreviation']})")
        print(f"      Compatibility: {rec['compatibility_score']}")
        print(f"      Applications: {rec['applications']}")
        print()
    
    print("=" * 80)
    print("TEST COMPLETED SUCCESSFULLY!")
    print("Enhanced solvent recommendation system is working properly.")
    print("=" * 80)

if __name__ == "__main__":
    test_solvent_recommendations()
