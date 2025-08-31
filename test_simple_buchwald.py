#!/usr/bin/env python3

"""Simple test to check Buchwald reaction consistency"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from enhanced_recommendation_engine import create_recommendation_engine

def test_buchwald_fields():
    """Test if Buchwald reactions have consistent field structures"""
    print("Testing Buchwald reaction field consistency...")
    
    engine = create_recommendation_engine()
    
    # Test with explicit Buchwald reaction type
    buchwald_smiles = "Brc1ccccc1.Nc1ccccc1>>Nc1ccccc1c1ccccc1"
    
    try:
        result = engine.get_recommendations(buchwald_smiles, "C-N Coupling - Buchwald-Hartwig")
        print(f"Status: {result.get('status', 'unknown')}")
        print(f"Reaction type: {result.get('reaction_type', 'none')}")
        
        # Check ligands
        ligands = result.get('ligand_recommendations', [])
        print(f"Ligand count: {len(ligands)}")
        if ligands:
            sample_ligand = ligands[0]
            print(f"Sample ligand fields: {list(sample_ligand.keys())}")
            
            # Check required fields
            required_fields = ['ligand', 'compatibility_score', 'applications', 'reaction_suitability']
            missing_fields = [field for field in required_fields if field not in sample_ligand]
            if missing_fields:
                print(f"❌ Missing ligand fields: {missing_fields}")
                return False
            else:
                print(f"✅ All required ligand fields present")
                
        # Check solvents
        solvents = result.get('solvent_recommendations', [])
        print(f"Solvent count: {len(solvents)}")
        if solvents:
            sample_solvent = solvents[0]
            print(f"Sample solvent fields: {list(sample_solvent.keys())}")
            
            # Check required fields
            required_fields = ['solvent', 'abbreviation', 'compatibility_score', 'applications', 'reaction_suitability']
            missing_fields = [field for field in required_fields if field not in sample_solvent]
            if missing_fields:
                print(f"❌ Missing solvent fields: {missing_fields}")
                return False
            else:
                print(f"✅ All required solvent fields present")
        
        if 'error' in result:
            print(f"❌ Error: {result['error']}")
            return False
            
        return True
            
    except Exception as e:
        print(f"❌ Exception: {e}")
        return False

def test_unsupported_reaction():
    """Test handling of unsupported reaction"""
    print("\nTesting unsupported reaction handling...")
    
    engine = create_recommendation_engine()
    
    try:
        result = engine.get_recommendations("CC(=O)C>>CC(O)C", "Hydrogenation - Catalytic")
        print(f"Status: {result.get('status', 'unknown')}")
        print(f"Ligand count: {len(result.get('ligand_recommendations', []))}")
        print(f"Solvent count: {len(result.get('solvent_recommendations', []))}")
        
        if 'error' in result:
            print(f"Error message: {result['error']}")
        
        # This should either have an error or empty recommendations
        has_error = 'error' in result
        has_empty_recommendations = (len(result.get('ligand_recommendations', [])) == 0 and 
                                   len(result.get('solvent_recommendations', [])) == 0)
        
        if has_error or has_empty_recommendations:
            print("✅ Unsupported reaction handled appropriately")
            return True
        else:
            print("⚠️ Unsupported reaction returned recommendations")
            return False
            
    except Exception as e:
        print(f"❌ Exception: {e}")
        return False

if __name__ == "__main__":
    success1 = test_buchwald_fields()
    success2 = test_unsupported_reaction()
    
    if success1 and success2:
        print("\n✅ All tests passed")
    else:
        print("\n❌ Some tests failed")
