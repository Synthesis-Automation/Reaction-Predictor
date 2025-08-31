#!/usr/bin/env python3
"""Test the improved reaction type display"""

from enhanced_recommendation_engine import create_recommendation_engine

def test_improved_display():
    """Test the improved reaction type display with the problematic reaction"""
    
    # The reaction that was showing just "Pd"
    test_reaction = 'Brc1ccc(C(F)(F)F)cc1.Nc1ccccc1>>FC(F)(F)c1ccc(Nc2ccccc2)cc1'
    
    print("🧪 TESTING IMPROVED REACTION TYPE DISPLAY")
    print("=" * 50)
    print(f"Reaction: {test_reaction}")
    print()
    
    engine = create_recommendation_engine()
    result = engine.get_recommendations(test_reaction, "Auto detect reaction type")
    
    print("Raw result data:")
    print(f"  reaction_type: {result.get('reaction_type')}")
    print(f"  detected_from: {result.get('detected_from')}")
    if 'rxn_insight_detection' in result:
        detection = result['rxn_insight_detection']
        print(f"  rxn_insight_detected: {detection.get('rxn_insight_detected')}")
        print(f"  mapped_to: {detection.get('mapped_to')}")
    print()
    
    # Test what the display function would show
    recommendations = result.get('recommendations', result)  # Fallback for structure
    
    # Simulate the display function logic
    rt = recommendations.get('reaction_type', result.get('reaction_type', 'Unknown')) or 'Unknown'
    rxn_insight_info = result.get('rxn_insight_detection', {})
    rxn_insight_detected = rxn_insight_info.get('rxn_insight_detected', '')
    
    print(f"Display would show:")
    print(f"  Detected Type: Should now show better info instead of just 'Pd'")
    print(f"  rxn-insight detected: {rxn_insight_detected}")
    print(f"  Our final type: {rt}")

if __name__ == "__main__":
    test_improved_display()
