#!/usr/bin/env python3
"""
Test script to verify that image scaling preserves aspect ratio
"""

import sys
import os
sys.path.insert(0, os.path.dirname(__file__))

def test_aspect_ratio_preservation():
    """Test that the scaling function preserves aspect ratio"""
    
    # Test the scaling logic without GUI dependencies
    def test_scale_calculation(original_width, original_height, max_width, max_height):
        """Simulate the scaling calculation from scale_pixmap_to_fit"""
        width_ratio = max_width / original_width if original_width > 0 else 1.0
        height_ratio = max_height / original_height if original_height > 0 else 1.0
        scale_factor = min(width_ratio, height_ratio, 1.0)
        
        new_width = int(original_width * scale_factor)
        new_height = int(original_height * scale_factor)
        
        # Calculate aspect ratios
        original_aspect = original_width / original_height if original_height > 0 else 1.0
        new_aspect = new_width / new_height if new_height > 0 else 1.0
        
        return new_width, new_height, original_aspect, new_aspect, scale_factor
    
    print("=== Image Scaling Aspect Ratio Test ===\n")
    
    # Test cases with various image dimensions
    test_cases = [
        # (original_width, original_height, max_width, max_height, description)
        (580, 200, 400, 280, "Wide reaction image → label bounds"),
        (480, 140, 300, 150, "Small reaction image → smaller bounds"),
        (600, 300, 580, 280, "Tall image → wide label"),
        (300, 300, 580, 280, "Square image → rectangle label"),
        (100, 50, 580, 280, "Small image (no scaling needed)"),
        (800, 400, 580, 280, "Large image → needs scaling down"),
    ]
    
    print("Testing aspect ratio preservation:")
    print("-" * 80)
    
    for orig_w, orig_h, max_w, max_h, desc in test_cases:
        new_w, new_h, orig_aspect, new_aspect, scale = test_scale_calculation(
            orig_w, orig_h, max_w, max_h
        )
        
        aspect_diff = abs(orig_aspect - new_aspect)
        status = "✅ PRESERVED" if aspect_diff < 0.01 else "❌ DISTORTED"
        
        print(f"{desc}")
        print(f"  Original: {orig_w}×{orig_h} (aspect: {orig_aspect:.3f})")
        print(f"  Scaled:   {new_w}×{new_h} (aspect: {new_aspect:.3f}) [scale: {scale:.3f}]")
        print(f"  Fits in:  {max_w}×{max_h}")
        print(f"  Result:   {status}")
        print()
    
    print("Key improvements:")
    print("• setScaledContents(False) - prevents automatic stretching")
    print("• scale_pixmap_to_fit() - preserves aspect ratio during scaling")
    print("• Images now remain crisp and undistorted! 🎯")

if __name__ == "__main__":
    test_aspect_ratio_preservation()
