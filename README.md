# Reaction Condition Predictor — Enhanced Edition

A GUI and CLI toolkit for reaction condition prediction with an enhanced recommendation engine, dataset-aware boosts, and a comprehensive sample browser (including Ullmann and Buchwald-Hartwig scopes).

## 🆕 What's New

### Enhanced Visual Reaction Browser

- Integrated image + text details panel for comprehensive viewing
- Instant visual feedback with high-quality chemical structure rendering
- Seamless browsing experience with visual confirmation

### Comprehensive Buchwald-Hartwig Collection

- 65+ specialized reactions covering diverse substrate scope
- Dedicated category filter for easy browsing
- Literature-based examples including challenging cases
- Perfect for testing the recommendation engine

### Evidence-aware recommendations (New)

- The enhanced engine mines small, local datasets under `data/reaction_dataset/` to nudge recommendations:
  - Ligands: evidence-aware boost (existing)
  - Solvents: new gentle boost (+0.05 to +0.15) for solvents frequently used in related reactions
  - Bases: new gentle boost (+0.05 to +0.15) applied and re-ranked when evidence is present
- Non-invasive: boosts are bounded and keep scores within [0, 1]

### Ullmann support and mapping

- GUI reaction labels normalized and mapped to internal types (e.g., "C-N Coupling - Ullmann" → `Ullmann`)
- Ullmann-specific solvent/base heuristics plus evidence-aware nudges

## 📁 Project Structure

```text
Reaction-Predictor/
├── enhanced_recommendation_engine.py   # Enhanced engine with evidence-aware boosts
├── recommendation_engine.py            # Legacy/base engine (still available)
├── simple_reaction_gui.py              # Main GUI application with RDKit visualization
├── predict_cli.py                      # CLI wrapper for predictions
├── reagents/
│   ├── ligand.py
│   ├── solvent.py
│   └── base.py
├── data/
│   ├── ligands.json
│   ├── solvents.json
│   └── reaction_dataset/               # Small CSV datasets used for evidence-aware nudges
├── scripts/
│   ├── smoke_ullmann.py                # Quick Ullmann smoke
│   └── export_solvents_full.py
├── requirements.txt
└── README.md
```

## 🔧 Dependencies

### Required Files (critical)

- ✅ enhanced_recommendation_engine.py — Enhanced condition recommendations (default entry point)
- ✅ simple_reaction_gui.py — Main application with reaction visualization
- ✅ reaction_types.py — Reaction type dropdown and taxonomy
- ✅ sample_reactions.py — Sample library
- ✅ recommendation_engine.py — Base system (kept for compatibility)

### Python packages

- pandas, numpy
- PyQt6 (GUI)
- RDKit (optional; for molecule/reaction rendering in GUI)
- scikit-learn, scipy (optional; used for similarity/scoring when available)
- matplotlib, networkx, openpyxl (optional; exports and visuals)

Install everything from the pinned requirements:

```powershell
pip install -r requirements.txt
```

## 🚀 Quick Start

### 1) Install dependencies

```powershell
pip install -r requirements.txt
```

Optional: RDKit for high-quality drawings in the GUI

```powershell
pip install rdkit
```

### 2) Launch the GUI

```powershell
python simple_reaction_gui.py
```

If imports fail in your shell, set PYTHONPATH to the repo root when launching:

```powershell
$env:PYTHONPATH='.'; python simple_reaction_gui.py
```

## 🎯 Buchwald-Hartwig Collection Features

### Accessing the Collection

1. Click "Browse Samples" in the main interface
2. Check the "Buchwald-Hartwig" checkbox
3. Browse 65+ specialized amination reactions
4. Select reactions for testing with the recommendation engine

## 🎯 Enhanced Visual Browser Features

### Real-Time Structure Preview

1. Click "Browse Samples" in the main interface
2. Select any category (e.g., "Buchwald-Hartwig")
3. Click any reaction in the list
4. See the reaction structure instantly in the details panel
5. Review both visual structure and text details
6. Double-click to load into the main interface

### Visual Browsing Experience

- Instant structure rendering using RDKit or placeholder graphics
- High-quality chemical drawings with proper reaction arrows
- Integrated preview panel showing structure + reaction details
- Real-time visual feedback for better reaction selection

### Collection Contents

- Primary Amines: 32 reactions with aniline derivatives
- Secondary Amines: 9 reactions with dialkyl/cyclic amines
- Chloro Substrates: 8 challenging aryl chloride examples
- Iodo Substrates: 3 highly reactive aryl iodide cases
- Benzyl Amines: 6 benzylic amine substrates
- Heteroaryl: 12 pyridine/thiazole/furan examples
- Sterically Hindered: 9 ortho-substituted difficult cases
- Pharmaceutical: 3 drug-like intermediate examples

### Testing the Collection

```powershell
# Test the collection
python test_buchwald_collection.py

# Demo the browser features
python demo_buchwald_browser.py
```

### Example Reactions

```text
# Simple diphenylamine formation
Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1

# Challenging ortho-substituted case
Brc1ccccc1C.Nc1ccccc1>>Cc1ccccc1Nc1ccccc1

# Heteroaryl substrate
Brc1ccncc1.Nc1ccccc1>>c1ccc(Nc2ccncc2)cc1

# Pharmaceutical intermediate
Brc1ccc(S(=O)(=O)N)cc1.Nc1ccccc1>>NS(=O)(=O)c1ccc(Nc2ccccc2)cc1
```

## 📋 System Requirements

- Python: 3.7 or higher
- Operating System: Windows, macOS, or Linux
- Memory: 50MB+ available RAM
- Display: 900x700 minimum resolution

## 🔍 Troubleshooting

### Import errors

If you get import errors:

```text
ImportError: No module named 'reaction_types'
ImportError: No module named 'sample_reactions'
```

Solution: Ensure these files are in the same directory as `simple_reaction_gui.py`.

### PyQt6 not found

```text
ModuleNotFoundError: No module named 'PyQt6'
```

Solution: Install PyQt6

```powershell
pip install PyQt6
```

### Missing functions

If you get errors about missing functions like `get_reaction_types()`:

Solution: Ensure you have the complete `reaction_types.py` and `sample_reactions.py` files.

### ModuleNotFoundError: enhanced_recommendation_engine

Python can't see the repo root on sys.path.

```powershell
$env:PYTHONPATH='.'; python scripts/smoke_ullmann.py
```

If you use a full Python path, keep the `PYTHONPATH='.'` prefix.

## 🎯 Features

- SMILES Input: Enter reaction SMILES notation
- Reaction Type Selection: 60+ reaction types
- Sample Library: 77+ example reactions
- Advanced Filtering: Search and category filters
- Dark Theme: Professional interface
- Cross-Platform: Works on Windows, macOS, Linux

## 📊 Sample Library Statistics

- Total Examples: 77 reactions
- Coupling Reactions: 21 examples
- Reduction Reactions: 9 examples
- Oxidation Reactions: 5 examples
- Other Reactions: 42 examples (Substitution, Elimination, etc.)

## 🔧 Development

### File dependencies

```text
simple_reaction_gui.py
    ├── imports PyQt6 (external)
    ├── imports reaction_types.get_reaction_types()
    └── imports sample_reactions.get_*_reactions()

reaction_types.py
    └── REACTION_TYPES list + helper functions

sample_reactions.py
    └── SAMPLE_REACTIONS list + filtering functions
```

### Adding new features

1. New Reaction Types: Edit `REACTION_TYPES` in `reaction_types.py`
2. New Sample Reactions: Edit `SAMPLE_REACTIONS` in `sample_reactions.py`
3. UI Changes: Modify `simple_reaction_gui.py`

## 📦 Creating Executable

To create a standalone executable:

```powershell
pip install pyinstaller
pyinstaller --onefile --windowed simple_reaction_gui.py
```

## 🖥️ Command-line predictor

In addition to the GUI, a CLI app is available and shares the same enhanced recommendation engine and export format.

```powershell
# Inline JSON
python predict_cli.py '{"reaction_smiles": "Brc1ccc2c(c1)cccc2C(C)(C)C.Nc1ccccc1>>CC(C)(C)c1cccc2cc(Nc3ccccc3)ccc12", "selected_reaction_type": "C-N Coupling - Ullmann"}'

# From a file
Get-Content .\input.json | python predict_cli.py
```

The CLI prints the JSON export to stdout.

## 🚀 Next Steps

This is a basic interface ready for:

- ✅ Prediction System Integration
- ✅ Database Connectivity
- ✅ Advanced Analytics
- ✅ Export/Import Features

## 📞 Support

For issues or questions:

1. Check that required files are present
2. Verify PyQt6 is installed
3. Ensure Python 3.7+ is being used
4. Check file permissions for the directory

## 🎉 Success Verification

When working correctly, you should see:

1. GUI window opens with "Reaction Condition Predictor" title
2. Dropdown shows 60+ reaction types
3. "Browse Samples" button opens library with 77+ reactions
4. Dark theme applied throughout interface

If any of these fail, check the dependencies above.

---

## 📈 Evidence-aware behavior (details)

- The enhanced engine scans small CSVs under `data/reaction_dataset/` for the selected reaction type.
- It builds lightweight frequency maps and applies a gentle, bounded boost:
  - Solvents: passed directly into `reagents/solvent.py` for +0.05..+0.15 nudges.
  - Bases: boosted in the engine post-processing by +0.05..+0.15 then re-ranked.
- Scores remain within [0,1]. If no evidence is found, behavior falls back to intrinsic compatibility scoring.

### Quick Ullmann smoke

```powershell
# Ensure repo root is on sys.path
$env:PYTHONPATH='.'; python scripts/smoke_ullmann.py

# Or inline
$env:PYTHONPATH='.'; python scripts/quick_smoke_inline.py
```
