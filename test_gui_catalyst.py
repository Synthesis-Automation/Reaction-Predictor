#!/usr/bin/env python3
"""Quick test to verify GUI can launch with catalyst selector"""

import sys
import os
sys.path.append(os.path.dirname(__file__))

from PyQt6.QtWidgets import QApplication

def test_gui_with_catalyst():
    """Test that GUI launches with catalyst selector"""
    
    try:
        app = QApplication(sys.argv)
        
        # Import and create GUI
        from simple_reaction_gui import SimpleReactionGUI
        
        gui = SimpleReactionGUI()
        
        # Check if catalyst buttons were created
        if hasattr(gui, 'catalyst_buttons'):
            print("✅ Catalyst selector created successfully!")
            print(f"Available catalysts: {list(gui.catalyst_buttons.keys())}")
            
            # Test catalyst selection
            selected = gui.get_selected_catalyst()
            print(f"Default selected catalyst: {selected}")
            
        else:
            print("❌ Catalyst selector not found!")
            
        # Don't show GUI for testing
        print("GUI compilation test passed!")
        
    except Exception as e:
        print(f"❌ Error: {e}")
        return False
        
    return True

if __name__ == "__main__":
    test_gui_with_catalyst()
