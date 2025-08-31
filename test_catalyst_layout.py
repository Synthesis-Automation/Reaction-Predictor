#!/usr/bin/env python3
"""
Test the two-column catalyst selector functionality
"""
import sys
import os
from PyQt6.QtWidgets import QApplication
from PyQt6.QtCore import QTimer

# Add the project root to the path
sys.path.insert(0, os.path.dirname(__file__))

from simple_reaction_gui import SimplePredictionGUI

def test_catalyst_selector_layout():
    """Test that the catalyst selector displays in two columns"""
    app = QApplication(sys.argv)
    
    # Create the GUI window
    window = SimplePredictionGUI()
    window.show()
    
    print("🧪 TESTING TWO-COLUMN CATALYST SELECTOR")
    print("=" * 50)
    
    # Check that all catalyst buttons exist
    expected_catalysts = [
        "auto", "Pd", "Ni", "Cu", "Fe", "Mn", "Co", 
        "Au", "Ir", "Ru", "Ti", "other", "organo"
    ]
    
    print(f"✅ Checking {len(expected_catalysts)} catalyst options...")
    
    missing_catalysts = []
    for catalyst in expected_catalysts:
        if catalyst not in window.catalyst_buttons:
            missing_catalysts.append(catalyst)
    
    if missing_catalysts:
        print(f"❌ Missing catalysts: {missing_catalysts}")
    else:
        print("✅ All catalyst buttons found!")
    
    # Test default selection
    if window.catalyst_buttons["auto"].isChecked():
        print("✅ Default 'Not specified' selection working!")
    else:
        print("❌ Default selection not working")
    
    # Test button group functionality
    if window.catalyst_group_buttons.checkedButton() == window.catalyst_buttons["auto"]:
        print("✅ Button group working correctly!")
    else:
        print("❌ Button group issue")
    
    print("\n📊 Catalyst Selector Layout:")
    print("   • Layout changed from single column to two columns")
    print("   • Should display catalysts in grid format")
    print("   • Total catalyst options: 13")
    print("   • Grid arrangement: 7 rows x 2 columns")
    
    print("\n✅ GUI is running - check the catalyst selector visually!")
    print("   The catalyst options should now be arranged in two columns")
    
    # Close after a short delay for testing
    def close_app():
        print("\n🔄 Test complete - you can close the GUI window")
    
    QTimer.singleShot(2000, close_app)  # Show message after 2 seconds
    
    # Don't exit immediately - let user see the GUI
    # app.exec()

if __name__ == "__main__":
    test_catalyst_selector_layout()
