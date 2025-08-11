#!/usr/bin/env python3
"""
Simple test to debug the HTTP server
"""

import os
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from simple_reaction_gui import SmilesDrawingServer
import time
import webbrowser

def test_server():
    """Test the HTTP server directly"""
    print("Starting HTTP server test...")
    
    # Create and start server
    server = SmilesDrawingServer(port=8002)  # Use different port to avoid conflicts
    
    try:
        server.start_server()
        print("Server started on port 8002")
        
        # Give server time to start
        time.sleep(1)
        
        # List files that should be available
        smiles_drawer_path = os.path.join(os.path.dirname(__file__), 'SMILES-Drawer')
        print(f"SMILES-Drawer path: {smiles_drawer_path}")
        print(f"Directory exists: {os.path.exists(smiles_drawer_path)}")
        
        if os.path.exists(smiles_drawer_path):
            print("Files in SMILES-Drawer:")
            for root, dirs, files in os.walk(smiles_drawer_path):
                level = root.replace(smiles_drawer_path, '').count(os.sep)
                indent = ' ' * 2 * level
                print(f"{indent}{os.path.basename(root)}/")
                subindent = ' ' * 2 * (level + 1)
                for file in files:
                    print(f"{subindent}{file}")
        
        print("\nTesting URLs:")
        print("http://localhost:8002/")
        print("http://localhost:8002/index.html")
        print("http://localhost:8002/index_enhanced.html")
        print("http://localhost:8002/js/enhanced_script.js")
        print("http://localhost:8002/jsme/jsme.nocache.js")
        
        print("\nTry opening one of these URLs in your browser to test.")
        print("Press Enter to stop the server...")
        input()
        
    except Exception as e:
        print(f"Error: {e}")
    finally:
        server.stop_server()
        print("Server stopped")

if __name__ == '__main__':
    test_server()
