#!/usr/bin/env python3
"""
Debug dataset name functionality
"""
import sys
import os
sys.path.insert(0, os.path.dirname(__file__))

from dataset_registry import DATASET_MAP, resolve_dataset_path

def debug_dataset_names():
    """Debug dataset name resolution"""
    print("🐛 DEBUGGING DATASET NAME RESOLUTION")
    print("=" * 50)
    
    print("\n1. DATASET_MAP contents:")
    for key, value in DATASET_MAP.items():
        print(f"   '{key}' -> '{value}'")
    
    print("\n2. Testing resolve_dataset_path function:")
    test_types = ['Cross-Coupling', 'Ullmann', 'C-N Coupling - Ullmann (Cu)', 'C-N Coupling - Buchwald-Hartwig (Pd)']
    
    for reaction_type in test_types:
        try:
            path = resolve_dataset_path(reaction_type)
            print(f"   '{reaction_type}' -> '{path}'")
            if path:
                basename = os.path.basename(path)
                print(f"      basename: '{basename}'")
        except Exception as e:
            print(f"   '{reaction_type}' -> ERROR: {e}")
    
    print("\n3. Testing direct DATASET_MAP lookups:")
    for reaction_type in test_types:
        if reaction_type in DATASET_MAP:
            print(f"   '{reaction_type}' found: '{DATASET_MAP[reaction_type]}'")
        else:
            print(f"   '{reaction_type}' NOT found")
            
            # Try normalized lookup
            normalized = reaction_type.replace(' (Pd)', '').replace(' (Cu)', '').replace(' (Ni)', '')
            if normalized in DATASET_MAP:
                print(f"      normalized '{normalized}' found: '{DATASET_MAP[normalized]}'")
            else:
                print(f"      normalized '{normalized}' NOT found")

if __name__ == "__main__":
    debug_dataset_names()
