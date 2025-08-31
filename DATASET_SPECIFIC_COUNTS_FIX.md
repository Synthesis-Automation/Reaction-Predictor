# Dataset-Specific Ligand and Solvent Count Fix

## ✅ ISSUE RESOLVED

### Original Problem
Both datasets were showing identical counts:
```
📊 Buchwald Dataset:
• Available Ligands: 123
• Available Solvents: 72

📊 Ullmann Dataset:  
• Available Ligands: 123
• Available Solvents: 72
```

### Root Cause
The enhanced recommendation engine was using **global database counts** instead of **dataset-specific counts**:

- **Global database**: 123 total ligands, 72 total solvents
- **Each reaction dataset**: ~10 specific ligands, ~10 specific solvents

### The Fix
Modified `enhanced_recommendation_engine.py` to use reaction-specific functions:

**BEFORE** (Global counts):
```python
ligand_df = create_ligand_dataframe()          # All 123 ligands
solvent_df = create_solvent_dataframe()        # All 72 solvents
ligands_available = len(ligand_df)             # Always 123
solvents_available = len(solvent_df)           # Always 72
```

**AFTER** (Dataset-specific counts):
```python
specific_ligands = get_reaction_specific_ligands(reaction_type)    # ~10 for each dataset
specific_solvents = get_reaction_specific_solvents(reaction_type)  # ~10 for each dataset
ligands_available = len(specific_ligands) if specific_ligands else 0
solvents_available = len(specific_solvents) if specific_solvents else 0
```

### Current Result
Now each dataset correctly shows its specific counts:
```
📊 Ullmann Dataset:
• Dataset: Ullman-2020-2024.jsonl
• Available Ligands: 10
• Available Solvents: 10

📊 Cross-Coupling Dataset:
• Dataset: Buchwald-2021-2024.jsonl  
• Available Ligands: 10
• Available Solvents: 10
```

### Why This Makes Sense

1. **Different Reaction Types**: Each reaction type (Ullmann vs Buchwald-Hartwig) has different optimal ligands and solvents
2. **Dataset Specificity**: The counts now reflect the actual reagents available for that specific reaction type
3. **Accurate Information**: Users see the true number of options for their selected reaction type, not the entire database

### Benefits

- **Accurate Representation**: Shows actual available options for each reaction type
- **Better User Understanding**: Users know exactly how many ligands/solvents are optimized for their reaction
- **Dataset Transparency**: Different datasets now show their true sizes
- **Realistic Expectations**: Users see the focused, curated options rather than overwhelming global counts

## ✅ IMPLEMENTATION COMPLETE

The system now correctly displays dataset-specific ligand and solvent counts, providing users with accurate information about the available options for each reaction type.
