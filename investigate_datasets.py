#!/usr/bin/env python3
"""
Investigate dataset-specific ligand and solvent counts
"""
import sys
import os
sys.path.insert(0, os.path.dirname(__file__))

def investigate_dataset_specificity():
    """Check if datasets have different ligand/solvent counts"""
    print("🔍 INVESTIGATING DATASET-SPECIFIC COUNTS")
    print("=" * 60)
    
    # Check global dataframes
    try:
        from reagents.ligand import create_ligand_dataframe
        from reagents.solvent import create_solvent_dataframe
        
        print("📊 Global Database Sizes:")
        ligand_df = create_ligand_dataframe()
        solvent_df = create_solvent_dataframe()
        print(f"   Total Ligands: {len(ligand_df)}")
        print(f"   Total Solvents: {len(solvent_df)}")
        
        # Check if there are reaction-specific functions
        try:
            from reagents.ligand import get_reaction_specific_ligands
            from reagents.solvent import get_reaction_specific_solvents
            print("\n✅ Found reaction-specific functions!")
            
            # Test for different reaction types
            test_reactions = ['Ullmann', 'Cross-Coupling', 'Buchwald-Hartwig']
            
            print("\n📋 Reaction-Specific Counts:")
            for reaction in test_reactions:
                print(f"\n   {reaction}:")
                try:
                    specific_ligands = get_reaction_specific_ligands(reaction)
                    specific_solvents = get_reaction_specific_solvents(reaction)
                    print(f"      Ligands: {len(specific_ligands) if specific_ligands else 'N/A'}")
                    print(f"      Solvents: {len(specific_solvents) if specific_solvents else 'N/A'}")
                except Exception as e:
                    print(f"      Error: {e}")
                    
        except ImportError:
            print("\n❌ No reaction-specific functions found")
            
    except Exception as e:
        print(f"❌ Error: {e}")
    
    # Check actual dataset files
    print("\n📁 Checking Dataset Files:")
    data_dir = os.path.join(os.path.dirname(__file__), 'data', 'reaction_dataset')
    if os.path.exists(data_dir):
        files = [f for f in os.listdir(data_dir) if f.endswith(('.jsonl', '.csv'))]
        print(f"   Dataset files: {files}")
        
        # Try to load and check a dataset file
        for file in files[:2]:  # Check first 2 files
            filepath = os.path.join(data_dir, file)
            print(f"\n   📄 {file}:")
            try:
                if file.endswith('.csv'):
                    import pandas as pd
                    df = pd.read_csv(filepath)
                    print(f"      Rows: {len(df)}")
                    if 'ligand' in df.columns:
                        unique_ligands = df['ligand'].nunique()
                        print(f"      Unique ligands: {unique_ligands}")
                    if 'solvent' in df.columns:
                        unique_solvents = df['solvent'].nunique()
                        print(f"      Unique solvents: {unique_solvents}")
                elif file.endswith('.jsonl'):
                    import json
                    ligands = set()
                    solvents = set()
                    count = 0
                    with open(filepath, 'r') as f:
                        for line in f:
                            if line.strip():
                                data = json.loads(line)
                                count += 1
                                if 'ligand' in data:
                                    ligands.add(str(data['ligand']))
                                if 'solvent' in data:
                                    solvents.add(str(data['solvent']))
                    print(f"      Reactions: {count}")
                    print(f"      Unique ligands: {len(ligands)}")
                    print(f"      Unique solvents: {len(solvents)}")
                    
            except Exception as e:
                print(f"      Error reading file: {e}")
    else:
        print(f"   ❌ Dataset directory not found: {data_dir}")

if __name__ == "__main__":
    investigate_dataset_specificity()
