#!/usr/bin/env python3

import sys
sys.path.append('.')
from enhanced_recommendation_engine import create_recommendation_engine

# Test with a separator line
engine = create_recommendation_engine()
separator_line = '───────── Cross Coupling Reactions ─────────'
smiles = 'Brc1ccc(C)cc1.Nc1ccccc1>>Cc1ccc(Nc2ccccc2)cc1'

print('=== TESTING SEPARATOR LINE HANDLING ===')
print(f'Separator: {separator_line}')
print(f'SMILES: {smiles}')

result = engine.get_recommendations(smiles, separator_line)

print(f'Status: {result.get("status")}')
print(f'Detected Type: {result.get("reaction_type")}')
print(f'Error: {result.get("error", "None")}')
print(f'Analysis Type: {result.get("analysis_type")}')

# Check if it gets treated as auto-detect
if 'auto_detection' in result:
    print('Treated as auto-detection: YES')
    auto = result['auto_detection']
    print(f'Auto-detection method: {auto.get("method")}')
else:
    print('Treated as auto-detection: NO')

# Check if _is_auto_detect recognizes it
print(f'Is separator considered auto-detect: {engine._is_auto_detect(separator_line)}')
