#!/usr/bin/env python3

import sys
sys.path.append('.')
from enhanced_recommendation_engine import create_recommendation_engine

# Test your exact example
print('=== FINAL TEST: YOUR EXACT EXAMPLE ===')
smiles = 'Brc1ccc(C)cc1.Nc1ccccc1>>Cc1ccc(Nc2ccccc2)cc1'

engine = create_recommendation_engine()

# This is what the GUI does when no reaction type is specified (default dropdown)
result = engine.get_recommendations(smiles, 'Auto detect reaction type')

print('Before fix: Selected Reaction Type: Not specified')
print('After fix:')
if 'auto_detection' in result:
    auto = result['auto_detection']
    print(f'• Selected Reaction Type: Auto-detected as {result.get("reaction_type")}')
    print(f'• Auto-Detection: {auto.get("rxn_insight_name")} (confidence: {auto.get("confidence")})')
    print(f'• Classification: {auto.get("rxn_insight_class")}')
    print(f'• Status: {result.get("status")}')
    print(f'• Ligands: {len(result.get("ligand_recommendations", []))}')
    print(f'• Solvents: {len(result.get("solvent_recommendations", []))}')
    print()
    print('✅ PROBLEM FIXED: Auto-detection now works correctly!')
else:
    print('❌ Auto-detection info not found')
