#!/usr/bin/env python3

"""Test Buchwald reaction handling consistency and error handling for unsupported reactions"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from enhanced_recommendation_engine import create_recommendation_engine

def test_buchwald_consistency():
    """Test if Buchwald reactions are handled consistently with other reaction types"""
    print("Testing Buchwald reaction consistency...")
    
    engine = create_recommendation_engine()
    
    # Test with explicit Buchwald reaction type
    buchwald_smiles = "Brc1ccccc1.Nc1ccccc1>>Nc1ccccc1c1ccccc1"  # Simple aryl halide + amine -> arylamine
    
    print("\n1. Testing explicit Buchwald reaction type:")
    try:
        result = engine.get_recommendations(buchwald_smiles, "C-N Coupling - Buchwald-Hartwig")
        print(f"   Status: {result.get('status', 'unknown')}")
        print(f"   Reaction type: {result.get('reaction_type', 'none')}")
        print(f"   Ligand recommendations count: {len(result.get('ligand_recommendations', []))}")
        print(f"   Solvent recommendations count: {len(result.get('solvent_recommendations', []))}")
        
        # Check if ligands have all required fields
        ligands = result.get('ligand_recommendations', [])
        if ligands:
            sample_ligand = ligands[0]
            required_fields = ['ligand', 'compatibility_score', 'applications', 'reaction_suitability']
            missing_fields = [field for field in required_fields if field not in sample_ligand]
            if missing_fields:
                print(f"   ❌ Missing ligand fields: {missing_fields}")
            else:
                print(f"   ✅ All required ligand fields present")
                
        # Check if solvents have all required fields
        solvents = result.get('solvent_recommendations', [])
        if solvents:
            sample_solvent = solvents[0]
            required_fields = ['solvent', 'abbreviation', 'compatibility_score', 'applications', 'reaction_suitability']
            missing_fields = [field for field in required_fields if field not in sample_solvent]
            if missing_fields:
                print(f"   ❌ Missing solvent fields: {missing_fields}")
            else:
                print(f"   ✅ All required solvent fields present")
        
        if 'error' in result:
            print(f"   ❌ Error: {result['error']}")
            
    except Exception as e:
        print(f"   ❌ Exception: {e}")

def test_unsupported_reaction():
    """Test error handling for unsupported reaction types"""
    print("\n2. Testing unsupported reaction type:")
    
    engine = create_recommendation_engine()
    
    # Test with a completely unsupported reaction type
    test_smiles = "CC(=O)C>>CC(O)C"  # Simple ketone reduction
    
    try:
        result = engine.get_recommendations(test_smiles, "Hydrogenation - Catalytic")
        print(f"   Status: {result.get('status', 'unknown')}")
        print(f"   Reaction type: {result.get('reaction_type', 'none')}")
        print(f"   Ligand recommendations count: {len(result.get('ligand_recommendations', []))}")
        print(f"   Solvent recommendations count: {len(result.get('solvent_recommendations', []))}")
        
        if 'error' in result:
            print(f"   ❌ Error: {result['error']}")
        elif len(result.get('ligand_recommendations', [])) == 0 and len(result.get('solvent_recommendations', [])) == 0:
            print(f"   ⚠️ No recommendations returned (but no error either)")
        else:
            print(f"   ⚠️ Recommendations returned for unsupported reaction type")
            
    except Exception as e:
        print(f"   ❌ Exception: {e}")

def test_reaction_type_detection():
    """Test what happens with auto-detection"""
    print("\n3. Testing auto-detection with different reaction patterns:")
    
    engine = create_recommendation_engine()
    
    test_cases = [
        ("Brc1ccccc1.Nc1ccccc1>>Nc1ccccc1c1ccccc1", "Aryl halide + amine (should detect Buchwald/Ullmann)"),
        ("CC(=O)O.NCc1ccccc1>>CC(=O)NCc1ccccc1", "Acid + amine (should detect amide formation)"),
        ("CCO.CC(=O)Cl>>CC(=O)OCC", "Alcohol + acid chloride (esterification - unsupported?)"),
    ]
    
    for smiles, description in test_cases:
        print(f"\n   Testing: {description}")
        try:
            result = engine.get_recommendations(smiles, "Auto-detect")
            print(f"   Detected type: {result.get('reaction_type', 'none')}")
            print(f"   Status: {result.get('status', 'unknown')}")
            print(f"   Recommendations: L={len(result.get('ligand_recommendations', []))}, S={len(result.get('solvent_recommendations', []))}")
            
            if 'error' in result:
                print(f"   ❌ Error: {result['error']}")
                
        except Exception as e:
            print(f"   ❌ Exception: {e}")

if __name__ == "__main__":
    test_buchwald_consistency()
    test_unsupported_reaction()
    test_reaction_type_detection()
    print("\nTest completed.")
