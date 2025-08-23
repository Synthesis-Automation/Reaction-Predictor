# Reaction Condition Predictor — Enhanced Edition

A GUI and CLI toolkit for reaction condition prediction with an enhanced recommendation engine, dataset-aware boosts, and a comprehensive sample browser (including Ullmann and Buchwald-Hartwig scopes).

## What changed (2025-08-24)

- Auto‑detect: added cross‑dataset similarity fallback using RDKit (Morgan FPs) to suggest ligands/solvents/bases and surface top related hits across all CSV datasets.
  - GUI shows a compact "Cross‑dataset similarity suggestions" section.
  - Related Reactions panel now falls back to these hits if dataset‑specific lookup finds none.
  - Export JSON includes a `general` block with suggestions and `top_hits` (each with `reaction_smiles`).
- Ullmann path refinements: demoted phosphines/NHCs (e.g., BenzP*, FerroTANE), added QUINAP penalty, and boosted N‑ligands (phen/bipy/L‑proline/diamines) with evidence support.
- Evidence‑aware boosting expanded to ligands, solvents, and bases, mined from local CSVs or analytics when available.
- GUI browser focus on coupling reactions: family/subtype filters and metal‑tagged labels (e.g., "Suzuki (Pd)", "Ullmann (Cu)"); detected type display includes metal.
- Solvent module cleanup: removed import‑time Excel side effects; Excel export is opt‑in via scripts.
- Export polish: simplified JSON with top_conditions; optional `dataset.analytics` snippet; new `general` block for Auto‑detect helper.

## 🆕 What's New

### General cross-dataset similarity (Auto‑detect)

- When the reaction type is left as Auto‑detect, the engine scans all CSVs under `data/reaction_dataset/` using RDKit similarity to surface:
  - Suggested ligands/solvents/bases (normalized 0..1 scores)
  - Top related reaction hits (now including `reaction_smiles`)
- The GUI shows a compact "Cross‑dataset similarity suggestions" block, and will fall back to these hits to populate the Related Reactions panel if dataset‑specific search yields none.

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
- Metal‑tagged labels throughout the browser (e.g., "Suzuki (Pd)", "Ullmann (Cu)") for clarity

### Solvent module cleanup

- Removed import‑time Excel writes; solvent Excel export is opt‑in via scripts.

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

Tip: Leave Reaction Type on "Auto‑detect" to see the cross‑dataset similarity helper alongside enhanced recommendations.
```

### 3) CLI: predict_cli.py usage

Use the CLI to generate a structured JSON export from the terminal.

Set PYTHONPATH (PowerShell):

```powershell
$env:PYTHONPATH='.'
```

Run with inline JSON (PowerShell):

```powershell
python predict_cli.py '{"reaction_smiles":"Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1","selected_reaction_type":"C-N Coupling - Ullmann"}'

Auto‑detect mode (shows `general` block in export):

```powershell
python predict_cli.py '{"reaction_smiles":"Brc1ccccc1.NH2c1ccccc1>>Nc1ccccc1","selected_reaction_type":"Auto-detect"}'
```
```

Read from a file and write output to exports (PowerShell):

```powershell
Get-Content .\input.json | python predict_cli.py | Set-Content -Path .\exports\prediction.json -Encoding utf8
```

Input payload (minimum):

```json
{
  "reaction_smiles": "<reactant(s)>><product>",
  "selected_reaction_type": "C-N Coupling - Ullmann"
}
```

Output (simplified):

- meta, input, detection
- dataset (includes `analytics` when available)
- top_conditions (up to 3 flattened sets with CAS where possible)
- related_reactions
- general (Auto‑detect helper): top ligands/solvents/bases and up to 10 similarity hits (with `reaction_smiles`)

Notes:

- `dataset.analytics` is loaded from `data/analytics/<Type>/latest.json` (currently `Ullmann`).
- Generate analytics via the VS Code tasks or `scripts/analyze_dataset.py`.

#### Capture output in other apps

`predict_cli.py` writes only JSON to stdout. Feed input via argv or stdin and capture stdout as a UTF‑8 string.

- PowerShell

```powershell
$env:PYTHONPATH='.'
$payload = '{"reaction_smiles":"R>>P","selected_reaction_type":"C-N Coupling - Ullmann"}'
$json = python predict_cli.py $payload; $json = ($json -join "`n"); $json

# Or from file
$json = Get-Content .\input.json | python predict_cli.py; $json = ($json -join "`n"); $json
```

- Python

```python
import subprocess
payload = '{"reaction_smiles":"R>>P","selected_reaction_type":"C-N Coupling - Ullmann"}'
out = subprocess.run([
  'python','predict_cli.py', payload
], capture_output=True, text=True, encoding='utf-8').stdout
```

- Node.js

```js
const { execFileSync } = require("child_process");
const payload =
  '{"reaction_smiles":"R>>P","selected_reaction_type":"C-N Coupling - Ullmann"}';
const out = execFileSync("python", ["predict_cli.py", payload], {
  encoding: "utf8",
});
```

- C# (.NET)

```csharp
using System.Diagnostics;
var payload = "{\"reaction_smiles\":\"R>>P\",\"selected_reaction_type\":\"C-N Coupling - Ullmann\"}";
var psi = new ProcessStartInfo("python", $"predict_cli.py \"{payload}\"")
{ RedirectStandardOutput = true, UseShellExecute = false };
using var p = Process.Start(psi)!;
string json = p.StandardOutput.ReadToEnd();
p.WaitForExit();
```

Tip: Check exit code (0 = success) and parse JSON on your side.

### Export JSON (Simplified)

Both the GUI and CLI emit a simplified JSON without a `recommendations` section. It focuses on:

- `meta`, `input`, `detection`, `dataset`
- `top_conditions`: up to 3 flattened condition sets (with CAS when available)
- `dataset.analytics`: dataset-driven summary (when available), including top ligands/solvents/bases, best pairs, numeric_stats (temperature_c, time_h, yield_pct), and typical_catalyst_loading if observed in suggestions
- `related_reactions`
- `general` (when Auto‑detect is used):
  - `ligands[]`, `solvents[]`, `bases[]` with normalized scores
  - `top_hits[]` with `ReactionID`, `ReactionType`, `CondKey`, `similarity`, and `reaction_smiles`

The GUI will also use `general.top_hits` to populate the Related Reactions panel when dataset lookup finds no matches.

The `dataset.analytics` block is sourced from `data/analytics/<Type>/latest.json` (currently `Ullmann`). Generate it via the analyzer (see Dataset analytics section below).

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
- Heteroaryl: 12 pyridine/thiazole/furan examples
- Sterically Hindered: 9 ortho-substituted difficult cases
- Pharmaceutical: 3 drug-like intermediate examples

# Test the collection

python test_buchwald_collection.py

python demo_buchwald_browser.py

````


```text
# Simple diphenylamine formation
# Challenging ortho-substituted case

## Dataset analytics (Milestone 1)

This repo includes a first-pass analytics pipeline (Ullmann only for now) to compute frequency summaries and numeric stats from CSV datasets.

- Run analyzer (PowerShell):

```powershell
$env:PYTHONPATH='.'
python scripts/analyze_dataset.py --reaction-type "Ullmann"
````

Artifacts will be written under `data/analytics/Ullmann/<timestamp>/summary.json` and `data/analytics/Ullmann/latest.json`.

Brc1ccccc1C.Nc1ccccc1>>Cc1ccccc1Nc1ccccc1

# Heteroaryl substrate

Brc1ccncc1.Nc1ccccc1>>c1ccc(Nc2ccncc2)cc1

# Pharmaceutical intermediate

Brc1ccc(S(=O)(=O)N)cc1.Nc1ccccc1>>NS(=O)(=O)c1ccc(Nc2ccccc2)cc1

````

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
````

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
