#!/usr/bin/env python3
"""
Debug solvent dataframe creation
"""
import sys
import os
sys.path.insert(0, os.path.dirname(__file__))

def debug_solvent_dataframe():
    """Debug why solvent count shows 0"""
    print("🧪 DEBUGGING SOLVENT DATAFRAME")
    print("=" * 50)
    
    # Test import
    try:
        from reagents.solvent import create_solvent_dataframe
        print("✅ Successfully imported create_solvent_dataframe")
    except ImportError as e:
        print(f"❌ Import error: {e}")
        return
    
    # Test dataframe creation
    try:
        print("\n📊 Creating solvent dataframe...")
        solvent_df = create_solvent_dataframe()
        print(f"✅ Solvent dataframe created successfully")
        print(f"   Shape: {solvent_df.shape}")
        print(f"   Columns: {list(solvent_df.columns)}")
        
        if len(solvent_df) > 0:
            print(f"   First few solvents:")
            print(solvent_df.head())
        else:
            print("❌ Dataframe is empty!")
            
    except Exception as e:
        print(f"❌ Error creating solvent dataframe: {e}")
        import traceback
        traceback.print_exc()
    
    # Test ligand dataframe for comparison
    try:
        from reagents.ligand import create_ligand_dataframe
        print("\n📊 Creating ligand dataframe for comparison...")
        ligand_df = create_ligand_dataframe()
        print(f"✅ Ligand dataframe shape: {ligand_df.shape}")
        
    except Exception as e:
        print(f"❌ Error creating ligand dataframe: {e}")
    
    # Check if solvent data file exists
    print("\n📁 Checking solvent data files...")
    data_dir = os.path.join(os.path.dirname(__file__), 'data')
    if os.path.exists(data_dir):
        files = os.listdir(data_dir)
        solvent_files = [f for f in files if 'solvent' in f.lower()]
        print(f"   Solvent-related files: {solvent_files}")
        
        # Check for solvents.json specifically
        solvents_json = os.path.join(data_dir, 'solvents.json')
        if os.path.exists(solvents_json):
            print(f"   ✅ solvents.json exists")
            # Check file size
            size = os.path.getsize(solvents_json)
            print(f"   File size: {size} bytes")
            
            if size < 100:
                print("   ⚠️  File seems very small - might be empty or corrupted")
                # Try to read first few lines
                try:
                    with open(solvents_json, 'r') as f:
                        content = f.read(200)  # First 200 chars
                        print(f"   First 200 chars: {repr(content)}")
                except Exception as e:
                    print(f"   ❌ Error reading file: {e}")
        else:
            print(f"   ❌ solvents.json not found")
    else:
        print(f"   ❌ Data directory not found: {data_dir}")

if __name__ == "__main__":
    debug_solvent_dataframe()
