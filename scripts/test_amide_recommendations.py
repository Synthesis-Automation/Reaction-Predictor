#!/usr/bin/env python3
"""
Test the amide formation recommendation system.
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from enhanced_recommendation_engine import create_recommendation_engine

def test_amide_recommendations():
    """Test amide formation recommendations."""
    
    print("Testing Amide Formation Recommendation System")
    print("=" * 50)
    
    # Create engine
    engine = create_recommendation_engine()
    
    # Test reaction: benzoic acid + aniline -> N-phenylbenzamide
    test_reaction = "O=C(O)c1ccccc1.Nc1ccccc1>>O=C(Nc1ccccc1)c1ccccc1"
    reaction_type = "Amide Formation - Acid + Amine"
    
    print(f"Test reaction: {test_reaction}")
    print(f"Reaction type: {reaction_type}")
    print()
    
    # Get recommendations
    result = engine.get_recommendations(test_reaction, reaction_type)
    
    print(f"Detected reaction type: {result.get('reaction_type', 'Unknown')}")
    print()
    
    # Display coupling reagents (stored as ligand_recommendations)
    print("=== COUPLING REAGENTS ===")
    ligands = result.get('ligand_recommendations', [])
    if ligands:
        for i, reagent in enumerate(ligands, 1):
            # Try both 'ligand' and 'name' fields for compatibility
            name = reagent.get('ligand', reagent.get('name', 'Unknown'))
            score = reagent.get('compatibility_score', 0)
            print(f"{i}. {name} (score: {score:.2f})")
    else:
        print("No coupling reagent recommendations found")
    print()
    
    # Display solvents
    print("=== SOLVENTS ===")
    solvents = result.get('solvent_recommendations', [])
    if solvents:
        for i, solvent in enumerate(solvents, 1):
            # Try both 'solvent' and 'name' fields for compatibility
            name = solvent.get('solvent', solvent.get('name', 'Unknown'))
            score = solvent.get('compatibility_score', 0)
            print(f"{i}. {name} (score: {score:.2f})")
    else:
        print("No solvent recommendations found")
    print()
    
    # Display bases
    print("=== BASES ===")
    bases = result.get('base_recommendations', [])
    if bases:
        for i, base in enumerate(bases, 1):
            # Try both 'base' and 'name' fields for compatibility  
            name = base.get('base', base.get('name', 'Unknown'))
            score = base.get('compatibility_score', 0)
            print(f"{i}. {name} (score: {score:.2f})")
    else:
        print("No base recommendations found")
    print()
    
    # Display conditions if available
    conditions = result.get('reaction_conditions', {})
    if conditions:
        print("=== RECOMMENDED CONDITIONS ===")
        if 'temperature' in conditions:
            temp = conditions['temperature']
            print(f"Temperature: {temp.get('recommended', 'N/A')} °C")
        if 'time' in conditions:
            time = conditions['time']
            print(f"Time: {time.get('recommended', 'N/A')} hours")
        print()
    
    print("Test completed successfully!")
    return result

if __name__ == "__main__":
    test_amide_recommendations()
