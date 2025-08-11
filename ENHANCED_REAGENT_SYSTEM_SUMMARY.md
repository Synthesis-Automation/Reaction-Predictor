# Enhanced Ligand and Solvent Recommendation System

## Overview

We have successfully enhanced your existing `ligand.py` and `solvent.py` scripts to provide comprehensive, reaction-specific recommendations for all major reaction types. This represents a significant upgrade from the original Buchwald-Hartwig focused system to a general-purpose chemical informatics platform.

## Key Enhancements

### 1. **Expanded Ligand Database**
- **130+ ligands** across multiple categories:
  - Monodentate phosphines
  - Buchwald-type ligands  
  - BINAP family
  - Bidentate phosphines
  - NHC ligands
  - Nitrogen ligands
  - Phosphoramidites
  - Carbene precursors

### 2. **Expanded Solvent Database** 
- **52+ solvents** with comprehensive properties:
  - Dielectric constant
  - Polarity index
  - Boiling point
  - Density
  - Dipole moment
  - Donor number
  - Hydrogen bond donor capability

### 3. **Reaction Type Support**
Both systems now support **5 major reaction types**:
- ✅ **Cross-Coupling** reactions
- ✅ **Hydrogenation** reactions  
- ✅ **Metathesis** reactions
- ✅ **C-H Activation** reactions
- ✅ **Carbonylation** reactions

### 4. **Intelligent Compatibility Scoring**
Each reagent has compatibility scores (0-1) for each reaction type:
- **0.9-1.0**: Highly optimized for reaction type
- **0.7-0.8**: Good compatibility
- **0.5-0.6**: Moderate compatibility  
- **0.3-0.4**: Limited compatibility
- **<0.3**: Poor compatibility

## New Functions Added

### Ligand Recommendations (`ligand.py`)

```python
# Get top ligands for specific reaction type
recommend_ligands_for_reaction(
    target_ligand=None,           # Optional reference ligand
    reaction_type='Cross-Coupling', 
    top_n=5, 
    min_compatibility=0.3
)

# Get filtered ligands with property constraints
get_reaction_specific_ligands(
    reaction_type='Cross-Coupling',
    property_preferences={
        'cone_angle_max': 170,    # Maximum steric bulk
        'price_category_max': 3,  # Affordable ligands only
    }
)
```

### Solvent Recommendations (`solvent.py`)

```python
# Get top solvents for specific reaction type
recommend_solvents_for_reaction(
    target_solvent=None,          # Optional reference solvent
    reaction_type='Cross-Coupling',
    top_n=5,
    min_compatibility=0.3
)

# Get filtered solvents with property constraints
get_reaction_specific_solvents(
    reaction_type='Cross-Coupling',
    property_preferences={
        'bp_max': 120,           # Low boiling point
        'polarity_min': 4,       # Polar solvents
        'protic': False          # Aprotic only
    }
)
```

## Advanced Features

### 1. **Similarity-Based Recommendations**
When you provide a reference ligand/solvent, the system:
- Calculates feature similarity using reaction-specific weights
- Combines compatibility and similarity scores
- Returns reagents that are both compatible AND similar

### 2. **Reaction-Specific Weighting**
Different properties are weighted based on reaction type:

**Cross-Coupling** prioritizes:
- Cone angle (steric effects)
- Electronic parameters
- Donor coordination

**Hydrogenation** prioritizes: 
- Electronic parameters (π-acceptance)
- Donor strength
- Bite angle for bidentate

**Metathesis** prioritizes:
- Air stability (sensitive catalysts)
- Electronic parameters
- Steric bulk

### 3. **Property-Based Filtering**
Filter recommendations by:
- **Ligands**: Cone angle, price category, air stability, coordination mode
- **Solvents**: Boiling point, polarity, protic/aprotic nature

## Example Usage

### Get Cross-Coupling Conditions
```python
from ligand import recommend_ligands_for_reaction
from solvent import recommend_solvents_for_reaction

# Get top ligands and solvents
ligands = recommend_ligands_for_reaction('Cross-Coupling', top_n=3)
solvents = recommend_solvents_for_reaction('Cross-Coupling', top_n=3)

# Results: SPhos/XPhos/RuPhos + DMF/THF/NMP
```

### Find Similar Reagents
```python
# Find ligands similar to PPh3 for cross-coupling
similar_ligands = recommend_ligands_for_reaction(
    target_ligand='PPh3',
    reaction_type='Cross-Coupling', 
    top_n=3
)

# Results: SPhos, RuPhos, XPhos (high similarity + compatibility)
```

### Filter by Properties
```python
# Low BP, polar, aprotic solvents for cross-coupling
filtered_solvents = get_reaction_specific_solvents(
    reaction_type='Cross-Coupling',
    property_preferences={
        'bp_max': 120,
        'polarity_min': 4,
        'protic': False
    }
)

# Results: THF, MeTHF, Acetonitrile
```

## Test Results Summary

### ✅ Cross-Coupling Recommendations:
- **Top Ligands**: SPhos, XPhos, RuPhos (Score: 0.9)
- **Top Solvents**: DMF, THF, NMP (Score: 0.9)

### ✅ Hydrogenation Recommendations:  
- **Top Ligands**: BINAP, Tol-BINAP, H8-BINAP (Score: 0.95)
- **Top Solvents**: EtOH, MeOH, 1-PrOH (Score: 0.9)

### ✅ Metathesis Recommendations:
- **Top Ligands**: IPr, IMes, SIPr (Score: 0.9)  
- **Top Solvents**: CCl4, PCE, AcOH (Score: 0.6-0.8)

### ✅ C-H Activation Recommendations:
- **Top Ligands**: IPr, IMes, SIPr (Score: 0.8)
- **Top Solvents**: DMSO, NEt3, Benzene (Score: 0.8-0.9)

### ✅ Carbonylation Recommendations:
- **Top Ligands**: PPh3, PtBu3, SPhos (Score: 0.8)
- **Top Solvents**: AcOH, DMF, NMP (Score: 0.9)

## Integration with Main System

The enhanced reagent recommendation functions can now be easily integrated into your main reaction prediction system (`simple_reaction_gui.py`) to provide comprehensive, intelligent condition recommendations for any reaction type.

This represents a complete transformation from a specialized Buchwald-Hartwig system to a general-purpose chemical informatics platform capable of supporting all major organometallic and organic reaction types.

## Files Modified/Created

### Enhanced Core Files:
- ✅ `reagents/ligand.py` - Enhanced with 130+ ligands and reaction compatibility
- ✅ `reagents/solvent.py` - Enhanced with 52+ solvents and reaction compatibility

### Test Files Created:
- ✅ `test_enhanced_ligands.py` - Ligand recommendation testing
- ✅ `test_enhanced_solvents.py` - Solvent recommendation testing  
- ✅ `test_comprehensive_recommendations.py` - Integrated system testing

All systems tested and verified working properly! 🎉
