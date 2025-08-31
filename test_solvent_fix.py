#!/usr/bin/env python3
"""
Test the fixed solvent count
"""
import sys
import os
sys.path.insert(0, os.path.dirname(__file__))

from enhanced_recommendation_engine import create_recommendation_engine

def test_fixed_solvent_count():
    """Test that solvent count now shows correctly"""
    print("✅ SOLVENT COUNT FIX VERIFICATION")
    print("=" * 50)
    
    engine = create_recommendation_engine()
    
    # Test Ullmann reaction
    result = engine.get_recommendations(
        'Brc1ccc(C#N)cc1.Nc1ccccc1>>N#Cc1ccc(Nc2ccccc2)cc1', 
        'Ullmann'
    )
    
    info = result.get('dataset_info', {})
    solvents_count = info.get('solvents_available', 'ERROR')
    
    print(f"🧪 Ullmann Reaction Test:")
    print(f"   Status: {result.get('status', 'Unknown')}")
    print(f"   Dataset: {info.get('dataset_name', 'Unknown')}")
    print(f"   Ligands: {info.get('ligands_available', 'Unknown')}")
    print(f"   Solvents: {solvents_count}")
    
    if solvents_count == 72:
        print("   ✅ SUCCESS: Solvent count is now correct!")
    elif solvents_count == 0:
        print("   ❌ STILL BROKEN: Solvent count is still 0")
    else:
        print(f"   ⚠️  UNEXPECTED: Solvent count is {solvents_count}")
    
    print(f"\n📊 Expected Database Coverage Display:")
    print(f"• Dataset: {info.get('dataset_name', 'Unknown')}")
    print(f"• Available Ligands: {info.get('ligands_available', 'Unknown')}")
    print(f"• Available Solvents: {solvents_count}")
    
    return solvents_count == 72

if __name__ == "__main__":
    success = test_fixed_solvent_count()
    if success:
        print("\n🎉 ISSUE RESOLVED!")
        print("The 'Available Solvents: 0' problem has been fixed.")
        print("GUI will now show the correct solvent count: 72")
    else:
        print("\n❌ Issue still exists - needs further investigation")
