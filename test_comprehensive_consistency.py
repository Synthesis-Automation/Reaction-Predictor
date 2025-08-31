#!/usr/bin/env python3

"""Comprehensive test for consistency across all reaction types and error handling"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from enhanced_recommendation_engine import create_recommendation_engine

def test_reaction_type(engine, smiles, reaction_type, description):
    """Test a specific reaction type and return results"""
    print(f"\n{description}:")
    try:
        result = engine.get_recommendations(smiles, reaction_type)
        status = result.get('status', 'unknown')
        rt = result.get('reaction_type', 'none')
        ligand_count = len(result.get('ligand_recommendations', []))
        solvent_count = len(result.get('solvent_recommendations', []))
        
        print(f"  Status: {status}")
        print(f"  Reaction type: {rt}")
        print(f"  Recommendations: L={ligand_count}, S={solvent_count}")
        
        if 'error' in result:
            print(f"  Error: {result['error']}")
            print(f"  Message: {result.get('message', 'N/A')}")
            return 'error'
        elif ligand_count > 0 or solvent_count > 0:
            # Check field consistency
            ligands = result.get('ligand_recommendations', [])
            solvents = result.get('solvent_recommendations', [])
            
            fields_ok = True
            if ligands:
                required_ligand_fields = ['ligand', 'compatibility_score', 'applications', 'reaction_suitability']
                missing = [f for f in required_ligand_fields if f not in ligands[0]]
                if missing:
                    print(f"  ❌ Missing ligand fields: {missing}")
                    fields_ok = False
                    
            if solvents:
                required_solvent_fields = ['solvent', 'abbreviation', 'compatibility_score', 'applications', 'reaction_suitability']
                missing = [f for f in required_solvent_fields if f not in solvents[0]]
                if missing:
                    print(f"  ❌ Missing solvent fields: {missing}")
                    fields_ok = False
                    
            if fields_ok:
                print("  ✅ All required fields present")
                return 'success'
            else:
                return 'field_error'
        else:
            print("  ⚠️ No recommendations returned")
            return 'empty'
            
    except Exception as e:
        print(f"  ❌ Exception: {e}")
        return 'exception'

def main():
    """Test all supported and unsupported reaction types"""
    print("Testing reaction type consistency and error handling...")
    
    engine = create_recommendation_engine()
    
    test_cases = [
        # Supported reaction types
        ("Brc1ccccc1.Nc1ccccc1>>Nc1ccccc1c1ccccc1", "C-N Coupling - Buchwald-Hartwig", "Buchwald-Hartwig (explicit)"),
        ("Brc1ccccc1.Nc1ccccc1>>Nc1ccccc1c1ccccc1", "C-N Coupling - Ullmann", "Ullmann (explicit)"),
        ("CC(=O)O.NCc1ccccc1>>CC(=O)NCc1ccccc1", "Amide Formation - Acid + Amine", "Amide Formation (explicit)"),
        
        # Auto-detection
        ("Brc1ccccc1.Nc1ccccc1>>Nc1ccccc1c1ccccc1", "Auto-detect", "Cross-coupling (auto-detect)"),
        ("CC(=O)O.NCc1ccccc1>>CC(=O)NCc1ccccc1", "Auto-detect", "Amide formation (auto-detect)"),
        
        # Unsupported reaction types
        ("CC(=O)C>>CC(O)C", "Hydrogenation - Catalytic", "Hydrogenation (unsupported)"),
        ("CC(=O)Cl.CCO>>CC(=O)OCC", "Esterification", "Esterification (unsupported)"),
        ("CC=CC>>CCC", "Alkene Hydrogenation", "Alkene Hydrogenation (unsupported)"),
        ("CC(C)=O>>CC(C)O", "Ketone Reduction", "Ketone Reduction (unsupported)"),
    ]
    
    results = {}
    for smiles, reaction_type, description in test_cases:
        result = test_reaction_type(engine, smiles, reaction_type, description)
        results[description] = result
    
    # Summary
    print("\n" + "="*60)
    print("SUMMARY:")
    
    supported_success = 0
    supported_total = 0
    unsupported_errors = 0
    unsupported_total = 0
    
    for desc, result in results.items():
        if 'unsupported' in desc.lower():
            unsupported_total += 1
            if result == 'error':
                unsupported_errors += 1
                print(f"✅ {desc}: Properly rejected")
            else:
                print(f"❌ {desc}: Should have been rejected but got {result}")
        else:
            supported_total += 1
            if result == 'success':
                supported_success += 1
                print(f"✅ {desc}: Working correctly")
            else:
                print(f"❌ {desc}: Failed with {result}")
    
    print(f"\nSupported reactions: {supported_success}/{supported_total} working")
    print(f"Unsupported reactions: {unsupported_errors}/{unsupported_total} properly rejected")
    
    if supported_success == supported_total and unsupported_errors == unsupported_total:
        print("\n🎉 ALL TESTS PASSED - System is working correctly!")
        return True
    else:
        print(f"\n❌ Some tests failed")
        return False

if __name__ == "__main__":
    main()
