#!/usr/bin/env python3

import sys
import os
sys.path.insert(0, os.path.abspath('.'))

from enhanced_recommendation_engine import create_recommendation_engine

def test_analytics_vs_evidence():
    """Test analytics vs evidence harvesting"""
    print("Testing Analytics vs Evidence Harvesting")
    print("=" * 45)
    
    engine = create_recommendation_engine()
    reaction_type = "Amide Formation - Acid + Amine"
    
    # Check analytics loading
    print("\n=== Testing Analytics Loading ===")
    priors = engine._load_analytics_summary(reaction_type)
    print(f"Priors loaded: {priors is not None}")
    if priors:
        print(f"Priors type: {type(priors)}")
        print(f"Priors keys: {list(priors.keys()) if isinstance(priors, dict) else 'Not a dict'}")
        
        # Try extracting from analytics
        try:
            analytics_ligands = engine._extract_priors(priors, 'ligands')
            print(f"Analytics ligands: {analytics_ligands}")
        except Exception as e:
            print(f"Analytics ligands error: {e}")
            
        try:
            analytics_solvents = engine._extract_priors(priors, 'solvents')
            print(f"Analytics solvents: {analytics_solvents}")
        except Exception as e:
            print(f"Analytics solvents error: {e}")
            
        try:
            analytics_bases = engine._extract_priors(priors, 'bases')
            print(f"Analytics bases: {analytics_bases}")
        except Exception as e:
            print(f"Analytics bases error: {e}")
    
    # Compare with direct evidence harvesting
    print("\n=== Direct Evidence Harvesting ===")
    evidence_ligands = engine._harvest_evidence_ligands(reaction_type)
    evidence_solvents = engine._harvest_evidence_solvents(reaction_type)
    evidence_bases = engine._harvest_evidence_bases(reaction_type)
    
    print(f"Evidence ligands count: {len(evidence_ligands)}")
    print(f"Evidence solvents count: {len(evidence_solvents)}")
    print(f"Evidence bases count: {len(evidence_bases)}")
    
    # Check analytics configuration
    print(f"\n=== Analytics Configuration ===")
    print(f"Analytics enabled: {engine._analytics_cfg['enabled']}")
    print(f"Apply to ligands: {engine._analytics_cfg['apply_to']['ligands']}")
    print(f"Apply to solvents: {engine._analytics_cfg['apply_to']['solvents']}")
    print(f"Apply to bases: {engine._analytics_cfg['apply_to']['bases']}")

if __name__ == "__main__":
    test_analytics_vs_evidence()
