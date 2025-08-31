#!/usr/bin/env python3

import sys
import os
sys.path.insert(0, os.path.abspath('.'))

from enhanced_recommendation_engine import create_recommendation_engine

def test_amide_logic_debugging():
    """Test amide formation logic with detailed debugging"""
    print("Testing Amide Formation Logic with Debugging")
    print("=" * 50)
    
    engine = create_recommendation_engine()
    reaction_type = "Amide Formation - Acid + Amine"
    
    # Check what happens in the logic
    print(f"Reaction type: {reaction_type}")
    print(f"Starts with 'amide formation': {(reaction_type or '').lower().startswith('amide formation')}")
    
    # Test evidence harvesting
    evidence_ligands = engine._harvest_evidence_ligands(reaction_type)
    evidence_solvents = engine._harvest_evidence_solvents(reaction_type)
    evidence_bases = engine._harvest_evidence_bases(reaction_type)
    
    print(f"\nEvidence harvested:")
    print(f"  Ligands: {len(evidence_ligands)} items")
    print(f"  Solvents: {len(evidence_solvents)} items")
    print(f"  Bases: {len(evidence_bases)} items")
    
    if evidence_ligands:
        print(f"  Top ligands: {list(evidence_ligands.keys())[:3]}")
    if evidence_solvents:
        print(f"  Top solvents: {list(evidence_solvents.keys())[:3]}")
    if evidence_bases:
        print(f"  Top bases: {list(evidence_bases.keys())[:3]}")
    
    # Test the recommendation logic manually
    print(f"\n=== Manual Ligand Conversion ===")
    if evidence_ligands:
        total_count = sum(evidence_ligands.values())
        print(f"Total count: {total_count}")
        
        ligands = []
        for reagent_name, count in sorted(evidence_ligands.items(), key=lambda x: x[1], reverse=True)[:5]:
            score = min(0.95, 0.5 + (count / total_count) * 0.45)
            ligand_rec = {
                'ligand': reagent_name,
                'compatibility_score': score,
                'applications': 'Coupling reagent for amide formation',
                'type': 'coupling_reagent'
            }
            ligands.append(ligand_rec)
            print(f"  {reagent_name}: count={count}, score={score:.3f}")
        
        print(f"Generated {len(ligands)} ligand recommendations")
    else:
        print("No evidence_ligands found!")

if __name__ == "__main__":
    test_amide_logic_debugging()
