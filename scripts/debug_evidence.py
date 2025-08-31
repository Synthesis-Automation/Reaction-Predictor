#!/usr/bin/env python3

import sys
import os
sys.path.insert(0, os.path.abspath('.'))

from enhanced_recommendation_engine import create_recommendation_engine

def test_evidence_harvesting():
    """Test evidence harvesting directly"""
    print("Testing Evidence Harvesting for Amide Formation")
    print("=" * 50)
    
    engine = create_recommendation_engine()
    
    # Test each evidence method directly
    print("\n=== Testing Ligand Evidence ===")
    ligand_evidence = engine._harvest_evidence_ligands("Amide Formation - Acid + Amine")
    print(f"Ligand evidence: {ligand_evidence}")
    
    print("\n=== Testing Solvent Evidence ===")
    solvent_evidence = engine._harvest_evidence_solvents("Amide Formation - Acid + Amine")
    print(f"Solvent evidence: {solvent_evidence}")
    
    print("\n=== Testing Base Evidence ===")
    base_evidence = engine._harvest_evidence_bases("Amide Formation - Acid + Amine")
    print(f"Base evidence: {base_evidence}")
    
    # Check if dataset file exists
    data_dir = os.path.join('.', 'data', 'reaction_dataset')
    print(f"\n=== Dataset Files ===")
    print(f"Data directory: {data_dir}")
    if os.path.exists(data_dir):
        files = [f for f in os.listdir(data_dir) if f.endswith('.jsonl')]
        print(f"JSONL files: {files}")
        
        # Check if amide dataset exists
        amide_file = os.path.join(data_dir, 'amide-formation-2021-2024.jsonl')
        if os.path.exists(amide_file):
            print(f"Amide dataset exists: {amide_file}")
            # Count lines
            with open(amide_file, 'r', encoding='utf-8') as f:
                line_count = sum(1 for line in f if line.strip())
            print(f"Lines in dataset: {line_count}")
        else:
            print(f"Amide dataset NOT found: {amide_file}")
    else:
        print(f"Data directory does not exist: {data_dir}")

if __name__ == "__main__":
    test_evidence_harvesting()
