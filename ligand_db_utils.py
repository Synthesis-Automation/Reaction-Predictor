"""
Ligand Database Management Utilities
===================================

Utilities for managing and updating the ligand families database.
"""

import pandas as pd
import os

def load_ligand_database(csv_path=None):
    """Load the ligand families database"""
    if csv_path is None:
        csv_path = os.path.join(os.path.dirname(__file__), 'data', 'ligand_families.csv')
    
    if os.path.exists(csv_path):
        return pd.read_csv(csv_path)
    else:
        print(f"Ligand database not found at {csv_path}")
        return None

def add_ligand_family(family_id, family_name, description, ligands, **kwargs):
    """Add a new ligand family to the database"""
    csv_path = os.path.join(os.path.dirname(__file__), 'data', 'ligand_families.csv')
    
    # Load existing database
    df = load_ligand_database(csv_path)
    if df is None:
        print("Could not load database")
        return False
    
    # Check if family already exists
    if family_id in df['family_id'].values:
        print(f"Family {family_id} already exists")
        return False
    
    # Create new row
    new_row = {
        'family_id': family_id,
        'family_name': family_name,
        'description': description,
        'ligands': ligands,
        'aliases': kwargs.get('aliases', ''),
        'ligand_type': kwargs.get('ligand_type', 'unknown'),
        'steric_bulk': kwargs.get('steric_bulk', 'unknown'),
        'electronic_property': kwargs.get('electronic_property', 'unknown'),
        'coordination_mode': kwargs.get('coordination_mode', 'unknown'),
        'bite_angle_range': kwargs.get('bite_angle_range', ''),
        'typical_applications': kwargs.get('typical_applications', 'general'),
        'cost_category': kwargs.get('cost_category', 'medium'),
        'performance_modifier': kwargs.get('performance_modifier', 1.0),
        'priority_rank': kwargs.get('priority_rank', 50),
        'notes': kwargs.get('notes', '')
    }
    
    # Add to dataframe
    df = pd.concat([df, pd.DataFrame([new_row])], ignore_index=True)
    
    # Save back to CSV
    df.to_csv(csv_path, index=False)
    print(f"Added ligand family: {family_id}")
    return True

def update_ligand_family(family_id, **updates):
    """Update an existing ligand family"""
    csv_path = os.path.join(os.path.dirname(__file__), 'data', 'ligand_families.csv')
    
    # Load existing database
    df = load_ligand_database(csv_path)
    if df is None:
        return False
    
    # Check if family exists
    if family_id not in df['family_id'].values:
        print(f"Family {family_id} not found")
        return False
    
    # Update the row
    for column, value in updates.items():
        if column in df.columns:
            df.loc[df['family_id'] == family_id, column] = value
    
    # Save back to CSV
    df.to_csv(csv_path, index=False)
    print(f"Updated ligand family: {family_id}")
    return True

def list_ligand_families():
    """List all ligand families in the database"""
    df = load_ligand_database()
    if df is None:
        return
    
    print("Ligand Families Database:")
    print("=" * 80)
    
    for _, row in df.iterrows():
        print(f"ID: {row['family_id']}")
        print(f"Name: {row['family_name']}")
        print(f"Description: {row['description']}")
        print(f"Ligands: {row['ligands']}")
        print(f"Type: {row.get('ligand_type', 'unknown')}")
        print(f"Priority: {row.get('priority_rank', 99)}")
        print("-" * 40)

def export_to_json(output_path=None):
    """Export ligand database to JSON format"""
    if output_path is None:
        output_path = os.path.join(os.path.dirname(__file__), 'data', 'ligand_families.json')
    
    df = load_ligand_database()
    if df is None:
        return False
    
    # Convert to JSON
    df.to_json(output_path, orient='records', indent=2)
    print(f"Exported ligand database to {output_path}")
    return True

def validate_database():
    """Validate the ligand database for consistency"""
    df = load_ligand_database()
    if df is None:
        return False
    
    issues = []
    
    # Check for required columns
    required_columns = ['family_id', 'family_name', 'description', 'ligands']
    for col in required_columns:
        if col not in df.columns:
            issues.append(f"Missing required column: {col}")
    
    # Check for duplicate family IDs
    if df['family_id'].duplicated().any():
        issues.append("Duplicate family IDs found")
    
    # Check for empty essential fields
    for idx, row in df.iterrows():
        if pd.isna(row['family_id']) or row['family_id'].strip() == '':
            issues.append(f"Row {idx}: Empty family_id")
        if pd.isna(row['ligands']) or row['ligands'].strip() == '':
            issues.append(f"Row {idx}: Empty ligands field")
    
    if issues:
        print("Database validation issues:")
        for issue in issues:
            print(f"  - {issue}")
        return False
    else:
        print("Database validation passed!")
        return True

if __name__ == "__main__":
    # Example usage
    list_ligand_families()
    validate_database()
