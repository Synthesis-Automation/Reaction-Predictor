#!/usr/bin/env python3
"""
Final Integration Test - Complete Structure Drawing Workflow
"""

import sys
import os
from PyQt6.QtWidgets import (
    QApplication,
    QMainWindow,
    QPushButton,
    QVBoxLayout,
    QWidget,
    QTextEdit,
    QLabel,
    QHBoxLayout,
)
from PyQt6.QtCore import Qt, QTimer

# Add the parent directory to sys.path to import our modules
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from simple_reaction_gui import SmilesDrawingDialog


class IntegrationTestWindow(QMainWindow):
    """Integration test window for the complete workflow"""

    def __init__(self):
        super().__init__()
        self.setWindowTitle("Structure Drawing Integration - Final Test")
        self.setGeometry(200, 200, 800, 600)

        # Create central widget and layout
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        layout = QVBoxLayout(central_widget)

        # Title
        title = QLabel("🧪 Structure Drawing Integration Test")
        title.setStyleSheet(
            "font-size: 24px; font-weight: bold; color: #4CAF50; padding: 20px;"
        )
        title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(title)

        # Instructions
        instructions = QLabel(
            """
<h3>Complete Workflow Test</h3>
<p>This test demonstrates the complete structure drawing integration workflow:</p>
<ol>
<li><b>Click "Test Structure Drawing"</b> to open the drawing dialog</li>
<li><b>Click "Launch Structure Editor"</b> in the dialog</li>
<li><b>Draw a molecule</b> in the web interface (benzene ring recommended)</li>
<li><b>Click "Save to Python"</b> in the web interface</li>
<li><b>Click "Use This Structure"</b> to complete the test</li>
<li><b>View the result</b> in the text area below</li>
</ol>
<p><b>Expected Result:</b> The SMILES string should appear in the result area</p>
"""
        )
        instructions.setWordWrap(True)
        instructions.setStyleSheet(
            """
            QLabel {
                background-color: #f0f8ff;
                border: 2px solid #4CAF50;
                border-radius: 8px;
                padding: 15px;
                margin: 10px;
                color: #333;
            }
        """
        )
        layout.addWidget(instructions)

        # Test button
        button_layout = QHBoxLayout()
        self.test_button = QPushButton("🎯 Test Structure Drawing")
        self.test_button.clicked.connect(self.test_structure_drawing)
        self.test_button.setMinimumHeight(50)
        self.test_button.setStyleSheet(
            """
            QPushButton {
                background-color: #4CAF50;
                border: none;
                color: white;
                padding: 15px 30px;
                border-radius: 8px;
                font-size: 16px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #45a049;
            }
            QPushButton:pressed {
                background-color: #3d8b40;
            }
        """
        )
        button_layout.addStretch()
        button_layout.addWidget(self.test_button)
        button_layout.addStretch()
        layout.addLayout(button_layout)

        # Results area
        results_label = QLabel("📊 Test Results:")
        results_label.setStyleSheet(
            "font-size: 16px; font-weight: bold; color: #333; margin-top: 20px;"
        )
        layout.addWidget(results_label)

        self.results_area = QTextEdit()
        self.results_area.setPlaceholderText(
            "Test results will appear here after completing the workflow..."
        )
        self.results_area.setMinimumHeight(200)
        self.results_area.setStyleSheet(
            """
            QTextEdit {
                background-color: #f8f9fa;
                border: 2px solid #ddd;
                border-radius: 8px;
                padding: 10px;
                font-family: monospace;
                font-size: 12px;
            }
        """
        )
        layout.addWidget(self.results_area)

        # Status area
        self.status_label = QLabel("Ready for testing")
        self.status_label.setStyleSheet(
            """
            QLabel {
                background-color: #e8f5e8;
                border: 1px solid #4CAF50;
                color: #2e7d32;
                padding: 8px;
                border-radius: 4px;
                font-weight: bold;
            }
        """
        )
        layout.addWidget(self.status_label)

        # Apply overall dark theme
        self.setStyleSheet(
            """
            QMainWindow {
                background-color: #ffffff;
            }
        """
        )

        # Log initial state
        self.log_message("🟢 Integration test window initialized")
        self.log_message("📁 Checking required files...")
        self.check_files()

    def check_files(self):
        """Check if all required files exist"""
        base_path = os.path.dirname(__file__)
        required_files = [
            "SMILES-Drawer/index_enhanced.html",
            "SMILES-Drawer/js/enhanced_script.js",
            "SMILES-Drawer/jsme/jsme.nocache.js",
            "SMILES-Drawer/css/styles.css",
        ]

        all_files_exist = True
        for file_path in required_files:
            full_path = os.path.join(base_path, file_path)
            if os.path.exists(full_path):
                self.log_message(f"✅ Found: {file_path}")
            else:
                self.log_message(f"❌ Missing: {file_path}")
                all_files_exist = False

        if all_files_exist:
            self.log_message("🎉 All required files found - ready for testing!")
            self.status_label.setText("All files present - ready for testing")
        else:
            self.log_message("⚠️ Some files are missing - test may fail")
            self.status_label.setText("Warning: Some files missing")
            self.status_label.setStyleSheet(
                """
                QLabel {
                    background-color: #fff3cd;
                    border: 1px solid #ffc107;
                    color: #856404;
                    padding: 8px;
                    border-radius: 4px;
                    font-weight: bold;
                }
            """
            )

    def log_message(self, message):
        """Add a message to the results area"""
        current_text = self.results_area.toPlainText()
        if current_text:
            new_text = current_text + "\n" + message
        else:
            new_text = message
        self.results_area.setPlainText(new_text)

        # Auto-scroll to bottom
        cursor = self.results_area.textCursor()
        cursor.movePosition(cursor.MoveOperation.End)
        self.results_area.setTextCursor(cursor)

    def test_structure_drawing(self):
        """Test the complete structure drawing workflow"""
        self.log_message("\n" + "=" * 60)
        self.log_message("🚀 Starting structure drawing test...")
        self.status_label.setText("Testing structure drawing workflow...")

        try:
            # Create and show the drawing dialog
            self.log_message("📋 Opening structure drawing dialog...")
            dialog = SmilesDrawingDialog(self)

            # Log dialog creation
            self.log_message("✅ Dialog created successfully")

            # Execute dialog
            self.log_message("⏳ Waiting for user to complete drawing workflow...")
            self.log_message("   1. Click 'Launch Structure Editor'")
            self.log_message("   2. Draw a molecule in the web interface")
            self.log_message("   3. Click 'Save to Python'")
            self.log_message("   4. Close browser and click 'Use This Structure'")

            result = dialog.exec()

            if result == dialog.DialogCode.Accepted:
                # Test completed successfully
                smiles = dialog.get_smiles()
                self.log_message("\n🎉 TEST COMPLETED SUCCESSFULLY!")
                self.log_message(f"📄 Received SMILES: {smiles}")

                if smiles:
                    self.log_message("✅ SMILES string is not empty")
                    self.log_message(
                        "✅ Structure drawing integration working correctly"
                    )
                    self.status_label.setText("✅ Test PASSED - Integration working!")
                    self.status_label.setStyleSheet(
                        """
                        QLabel {
                            background-color: #d4edda;
                            border: 1px solid #c3e6cb;
                            color: #155724;
                            padding: 8px;
                            border-radius: 4px;
                            font-weight: bold;
                        }
                    """
                    )
                else:
                    self.log_message("⚠️ SMILES string is empty")
                    self.log_message(
                        "   This might indicate the user didn't draw anything"
                    )
                    self.status_label.setText("⚠️ Test completed but no structure drawn")
            else:
                # Test was cancelled
                self.log_message("\n⏹️ Test cancelled by user")
                self.log_message("   Dialog was closed without completing the workflow")
                self.status_label.setText("Test cancelled by user")

        except Exception as e:
            # Test failed
            self.log_message(f"\n❌ TEST FAILED!")
            self.log_message(f"💥 Error: {str(e)}")
            self.status_label.setText("❌ Test FAILED - Check results")
            self.status_label.setStyleSheet(
                """
                QLabel {
                    background-color: #f8d7da;
                    border: 1px solid #f5c6cb;
                    color: #721c24;
                    padding: 8px;
                    border-radius: 4px;
                    font-weight: bold;
                }
            """
            )

        self.log_message("=" * 60)


def main():
    """Main function"""
    app = QApplication(sys.argv)

    print("🧪 Starting Structure Drawing Integration Test")
    print("=" * 50)
    print("This will test the complete workflow of:")
    print("1. Opening the structure drawing dialog")
    print("2. Launching the web-based editor")
    print("3. Drawing a structure")
    print("4. Transferring SMILES back to the GUI")
    print("5. Completing the integration")
    print("=" * 50)

    window = IntegrationTestWindow()
    window.show()

    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
