# Reaction Condition Recommendation System

A modular and extensible system for recommending reaction conditions based on literature data and catalyst family analysis.

## Features

### Current Implementation (v1.0)
- **Buchwald-Hartwig Amination Specialist**: Provides specific catalyst recommendations based on literature data
- **Catalyst Family Analysis**: Groups related catalysts and recommends alternatives
- **Performance-Based Ranking**: Ranks recommendations by yield statistics and literature support
- **Extensible Architecture**: Designed to easily add new reaction types
- **🆕 Reaction Scheme Visualization**: Real-time reaction drawing using RDKit

### GUI Interface
- Dark-themed PyQt6 interface
- Real-time condition recommendations
- **🆕 Interactive reaction scheme display**
- Sample reaction browser
- Example reactions for testing

## System Architecture

### Core Components

#### 1. Recommendation Engine (`recommendation_engine.py`)
```python
RecommendationEngine
├── BuchwaldHartwigRecommender  # Specific for B-H aminations
├── GeneralRecommender          # Fallback for other reactions
└── Future: SuzukiRecommender, HeckRecommender, etc.
```

#### 2. Catalyst Family System
- **Ligand-Centric Organization**: Ligands have 60% weight in similarity calculations
- **Metal Precursor Flexibility**: Different Pd sources (Pd(OAc)2, Pd2(dba)3, PdCl2) treated as similar
- **Family Relationships**: Related catalyst families provide backup options

#### 3. Performance Analysis
- **Yield Statistics**: Average, maximum, success rates
- **Confidence Scoring**: Based on amount of literature data
- **Literature Support**: Number of reactions and references

### Catalyst Families (Buchwald-Hartwig)

1. **Bulky Phosphines** (`bulky_phosphines`)
   - Ligands: XPhos, SPhos, RuPhos, BrettPhos, tBuXPhos
   - Best for: Challenging, sterically hindered substrates

2. **Moderate Phosphines** (`moderate_phosphines`)
   - Ligands: DavePhos, JohnPhos, CyJohnPhos
   - Best for: General applications

3. **Bidentate Ligands** (`bidentate_ligands`)
   - Ligands: dppf, dppp, dppe, BINAP, Xantphos
   - Best for: Specific electronic requirements

4. **Simple Phosphines** (`simple_phosphines`)
   - Ligands: PPh3, P(o-tolyl)3, PCy3
   - Best for: Simple, electron-rich substrates

## Usage

### Basic Usage
```python
from recommendation_engine import create_recommendation_engine

# Create engine (auto-loads available datasets)
engine = create_recommendation_engine()

# Get recommendations
result = engine.get_recommendations(
    "c1ccc(Br)cc1.CN>>c1ccc(N(C)C)cc1",  # Reaction SMILES
    "Buchwald-Hartwig Amination"          # Optional reaction type
)

# Result contains ranked catalyst family recommendations
for rec in result['recommendations']:
    print(f"Family: {rec['family_name']}")
    print(f"Catalyst: {rec['recommended_catalyst']}/{rec['recommended_ligand']}")
    print(f"Performance: {rec['performance']['avg_yield']:.1f}% average yield")
```

### GUI Usage
```bash
python simple_reaction_gui.py
```

1. Enter reaction SMILES or use "Load Example"
2. **🆕 Watch the reaction scheme appear automatically**
3. Select reaction type (or use "Auto-detect")
4. Click "Predict Conditions"
5. View ranked catalyst recommendations

### Testing
```bash
python test_recommendations.py
```

## Data Requirements

### Current Dataset
- **Buchwald Reactions**: `data/buchwald_reactions.csv`
  - ~500 literature reactions
  - Columns: ReactantSMILES, ProductSMILES, CoreGeneric, Ligand, Solvent, Yield_%

### Expected Format for New Datasets
```csv
ReactantSMILES,ProductSMILES,CoreGeneric,Ligand,Solvent,Base,Yield_%,Reference
[substrate_smiles],[product_smiles],[catalyst],[ligand],[solvent],[base],[yield],[reference]
```

## Extending the System

### Adding New Reaction Types

1. **Create New Recommender Class**:
```python
class SuzukiRecommender(ReactionRecommender):
    def can_handle_reaction(self, reaction_smiles, reaction_type=None):
        # Detection logic for Suzuki reactions
        return "suzuki" in reaction_type.lower() if reaction_type else False
    
    def get_recommendations(self, reaction_smiles, top_k=5):
        # Suzuki-specific recommendation logic
        pass
```

2. **Register in Engine**:
```python
def create_recommendation_engine(data_dir="data"):
    engine = RecommendationEngine()
    
    # Add Suzuki recommender
    suzuki_rec = SuzukiRecommender()
    engine.add_recommender(suzuki_rec, "data/suzuki_reactions.csv")
    
    return engine
```

3. **Update GUI** (if needed):
```python
def _format_suzuki_recommendations(self, result):
    # Format Suzuki-specific results
    pass
```

### Key Design Principles

1. **Modular**: Each reaction type has its own recommender
2. **Extensible**: Easy to add new reaction types without changing existing code
3. **Data-Driven**: Recommendations based on literature performance data
4. **Family-Aware**: Understands relationships between catalyst components
5. **Confidence-Scored**: Provides reliability indicators for recommendations

### Performance Metrics

- **Average Yield**: Mean yield across literature examples
- **Success Rate**: Percentage of reactions with ≥80% yield
- **Confidence**: Based on number of supporting literature examples
- **Ranking Score**: Weighted combination of all metrics

### Future Enhancements

1. **Structural Similarity**: Use RDKit molecular fingerprints for substrate matching
2. **Machine Learning**: Train models to predict optimal conditions
3. **More Reaction Types**: Suzuki, Heck, Sonogashira, etc.
4. **Condition Optimization**: Temperature, solvent, base recommendations
5. **Literature Integration**: Direct links to original papers

## Requirements

```
PyQt6
pandas
numpy
rdkit-pypi (optional, for reaction visualization)
```

### Installation Notes
- **RDKit** is optional but recommended for chemical structure visualization
- Without RDKit, the system will show placeholder reaction schemes
- All recommendation features work regardless of RDKit installation

## Installation

```bash
# Clone or download the repository
cd Reaction-Predictor

# Install dependencies
pip install PyQt6 pandas numpy

# Run the application
python simple_reaction_gui.py
```

---

**Note**: This is version 1.0 focused on Buchwald-Hartwig aminations. The system is designed to be easily extended for other catalytic and non-catalytic reactions.
