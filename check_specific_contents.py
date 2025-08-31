#!/usr/bin/env python3
"""
Check what the reaction-specific ligands and solvents actually contain
"""
import sys
import os
sys.path.insert(0, os.path.dirname(__file__))

def check_reaction_specific_contents():
    """Check the actual contents of reaction-specific ligands and solvents"""
    print("🔍 CHECKING REACTION-SPECIFIC CONTENTS")
    print("=" * 60)
    
    try:
        from reagents.ligand import get_reaction_specific_ligands
        from reagents.solvent import get_reaction_specific_solvents
        
        test_reactions = ['Ullmann', 'Cross-Coupling', 'C-N Coupling - Buchwald-Hartwig']
        
        for reaction in test_reactions:
            print(f"\n🧪 {reaction}:")
            print("-" * 40)
            
            # Check ligands
            try:
                ligands = get_reaction_specific_ligands(reaction)
                print(f"📎 Ligands ({len(ligands)} total):")
                if ligands:
                    for i, ligand in enumerate(ligands[:10], 1):  # Show first 10
                        # Check if ligand is a dict or string
                        if isinstance(ligand, dict):
                            name = ligand.get('name', ligand.get('ligand', str(ligand)))
                            score = ligand.get('score', ligand.get('compatibility_score', ''))
                            if score:
                                print(f"   {i:2d}. {name} (score: {score})")
                            else:
                                print(f"   {i:2d}. {name}")
                        else:
                            print(f"   {i:2d}. {ligand}")
                else:
                    print("   No ligands found")
                    
            except Exception as e:
                print(f"   ❌ Error getting ligands: {e}")
            
            # Check solvents
            try:
                solvents = get_reaction_specific_solvents(reaction)
                print(f"\n🧪 Solvents ({len(solvents)} total):")
                if solvents:
                    for i, solvent in enumerate(solvents[:10], 1):  # Show first 10
                        # Check if solvent is a dict or string
                        if isinstance(solvent, dict):
                            name = solvent.get('name', solvent.get('solvent', str(solvent)))
                            score = solvent.get('score', solvent.get('compatibility_score', ''))
                            if score:
                                print(f"   {i:2d}. {name} (score: {score})")
                            else:
                                print(f"   {i:2d}. {name}")
                        else:
                            print(f"   {i:2d}. {solvent}")
                else:
                    print("   No solvents found")
                    
            except Exception as e:
                print(f"   ❌ Error getting solvents: {e}")
                
    except ImportError as e:
        print(f"❌ Import error: {e}")

if __name__ == "__main__":
    check_reaction_specific_contents()
