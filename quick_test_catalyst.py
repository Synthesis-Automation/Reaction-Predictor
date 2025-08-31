#!/usr/bin/env python3
"""Quick test of catalyst selection for the user's reaction"""

import sys
import os
sys.path.append(os.path.dirname(__file__))

# Direct import without GUI
from enhanced_recommendation_engine import EnhancedRecommendationEngine

def quick_test():
    """Quick test without GUI overhead"""
    
    reaction_smiles = 'Brc1ccc(C#N)cc1.Nc1ccccc1>>N#Cc1ccc(Nc2ccccc2)cc1'
    
    print("Quick test of catalyst selection:")
    print(f"SMILES: {reaction_smiles}")
    print()
    
    # Create engine directly
    engine = EnhancedRecommendationEngine()
    
    # Test detection
    result = engine.get_recommendations(reaction_smiles, "Auto detect reaction type")
    
    print(f"Status: {result.get('status')}")
    print(f"Analysis Type: {result.get('analysis_type')}")
    
    if result.get('analysis_type') == 'catalyst_selection_needed':
        print("✅ SUCCESS - Catalyst selection working!")
        print(f"Message: {result.get('message')}")
        print("Options:")
        for cat, rxn in result.get('catalyst_options', {}).items():
            print(f"  {cat}: {rxn}")
    else:
        print(f"❌ Unexpected result:")
        print(f"  Analysis type: {result.get('analysis_type')}")
        print(f"  Reaction type: {result.get('reaction_type')}")
        if result.get('error'):
            print(f"  Error: {result.get('error')}")

if __name__ == "__main__":
    quick_test()
