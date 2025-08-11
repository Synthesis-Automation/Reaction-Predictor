# Simple Reaction Condition Predictor - Enhanced Edition

A standalone GUI application for reaction condition prediction with comprehensive sample reactions library and specialized Buchwald-Hartwig amination collection.

## 🆕 NEW FEATURES

### **Enhanced Visual Reaction Browser**
- **Real-time structure preview** when selecting reactions
- **Integrated image + text details** panel for comprehensive viewing
- **Instant visual feedback** with high-quality chemical structure rendering
- **Seamless browsing experience** with visual confirmation

### **Comprehensive Buchwald-Hartwig Collection**
- **65+ specialized reactions** covering diverse substrate scope
- **Dedicated category filter** for easy browsing
- **Literature-based examples** including challenging cases
- **Perfect for testing** the recommendation engine

## 📁 Project Structure

```
Reaction-Predictor/
├── simple_reaction_gui.py          # Main GUI application with RDKit visualization
├── reaction_types.py               # Reaction type definitions (REQUIRED)
├── sample_reactions.py             # Enhanced sample library: 77+ general + 65+ B-H specific
├── recommendation_engine.py        # Advanced recommendation engine with catalyst families
├── run_simple_gui.bat             # Windows launcher script
├── test_buchwald_collection.py    # Test script for new B-H collection
├── demo_buchwald_browser.py       # Demo for new browser features
├── requirements.txt               # Python dependencies
└── README.md                      # This file
```

## 🔧 Dependencies

### **Required Files (CRITICAL):**
- ✅ **`simple_reaction_gui.py`** - Main application with reaction visualization
- ✅ **`reaction_types.py`** - Provides reaction type dropdown options
- ✅ **`sample_reactions.py`** - Provides 140+ sample reactions (77 general + 65 Buchwald-Hartwig)
- ✅ **`recommendation_engine.py`** - Advanced condition recommendation system

### **Python Packages:**
- **PyQt6** - GUI framework
- **pandas** - Dataset handling
- **RDKit** (optional) - Chemical structure visualization
- **Standard Library** - threading, time, re (built-in)

## 🚀 Quick Start

### **1. Install Python Dependencies**
```bash
pip install PyQt6
```

### **2. Run the Application**

**Option A: Python Command**
```bash
python simple_reaction_gui.py
```

**Option B: Windows Batch File**
```bash
run_simple_gui.bat
```

**Option C: Double-click** `simple_reaction_gui.py` (if Python is associated)

## 🎯 Buchwald-Hartwig Collection Features

### **Accessing the Collection**
1. Click **"Browse Samples"** button in the main interface
2. Check the **"Buchwald-Hartwig"** checkbox (new!)
3. Browse 65+ specialized amination reactions
4. Select reactions for testing with the recommendation engine

## 🎯 Enhanced Visual Browser Features

### **Real-Time Structure Preview**
1. Click **"Browse Samples"** button in the main interface
2. Select any category (e.g., **"Buchwald-Hartwig"**)
3. **Click on any reaction** in the list
4. **See the reaction structure instantly** in the details panel
5. Review both **visual structure and text details**
6. **Double-click to load** into main interface

### **Visual Browsing Experience**
- **Instant structure rendering** using RDKit or placeholder graphics
- **High-quality chemical drawings** with proper reaction arrows
- **Integrated preview panel** showing structure + reaction details
- **Real-time visual feedback** for better reaction selection

### **Collection Contents**
- **Primary Amines**: 32 reactions with aniline derivatives
- **Secondary Amines**: 9 reactions with dialkyl/cyclic amines  
- **Chloro Substrates**: 8 challenging aryl chloride examples
- **Iodo Substrates**: 3 highly reactive aryl iodide cases
- **Benzyl Amines**: 6 benzylic amine substrates
- **Heteroaryl**: 12 pyridine/thiazole/furan examples
- **Sterically Hindered**: 9 ortho-substituted difficult cases
- **Pharmaceutical**: 3 drug-like intermediate examples

### **Testing the Collection**
```bash
# Test the new collection
python test_buchwald_collection.py

# Demo the browser features  
python demo_buchwald_browser.py
```

### **Example Reactions**
```
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

- **Python**: 3.7 or higher
- **Operating System**: Windows, macOS, or Linux
- **Memory**: 50MB+ available RAM
- **Display**: 900x700 minimum resolution

## 🔍 Troubleshooting

### **Import Errors**
If you get import errors:

```python
ImportError: No module named 'reaction_types'
ImportError: No module named 'sample_reactions'
```

**Solution**: Ensure these files are in the same directory as `simple_reaction_gui.py`

### **PyQt6 Not Found**
```bash
ModuleNotFoundError: No module named 'PyQt6'
```

**Solution**: Install PyQt6
```bash
pip install PyQt6
```

### **Missing Functions**
If you get errors about missing functions like `get_reaction_types()`:

**Solution**: Ensure you have the complete `reaction_types.py` and `sample_reactions.py` files

## 🎯 Features

- **SMILES Input**: Enter reaction SMILES notation
- **Reaction Type Selection**: 60+ reaction types
- **Sample Library**: 77+ example reactions
- **Advanced Filtering**: Search and category filters
- **Dark Theme**: Professional interface
- **Cross-Platform**: Works on Windows, macOS, Linux

## 📊 Sample Library Statistics

- **Total Examples**: 77 reactions
- **Coupling Reactions**: 21 examples
- **Reduction Reactions**: 9 examples  
- **Oxidation Reactions**: 5 examples
- **Other Reactions**: 42 examples (Substitution, Elimination, etc.)

## 🔧 Development

### **File Dependencies:**

```python
simple_reaction_gui.py
    ├── imports PyQt6 (external)
    ├── imports reaction_types.get_reaction_types()
    └── imports sample_reactions.get_*_reactions()

reaction_types.py
    └── REACTION_TYPES list + helper functions

sample_reactions.py  
    └── SAMPLE_REACTIONS list + filtering functions
```

### **Adding New Features:**
1. **New Reaction Types**: Edit `REACTION_TYPES` in `reaction_types.py`
2. **New Sample Reactions**: Edit `SAMPLE_REACTIONS` in `sample_reactions.py`
3. **UI Changes**: Modify `simple_reaction_gui.py`

## 📦 Creating Executable

To create a standalone executable:

```bash
pip install pyinstaller
pyinstaller --onefile --windowed simple_reaction_gui.py
```

## 🚀 Next Steps

This is a basic interface ready for:
- ✅ **Prediction System Integration**
- ✅ **Database Connectivity** 
- ✅ **Advanced Analytics**
- ✅ **Export/Import Features**

## 📞 Support

For issues or questions:
1. Check that all 3 required files are present
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
