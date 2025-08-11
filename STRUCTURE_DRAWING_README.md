# Structure Drawing Integration

## Overview
The Reaction Predictor GUI now includes integrated molecular structure drawing capability using the JSME (JavaScript Molecular Editor). This allows users to draw chemical structures visually instead of manually typing SMILES strings.

## Features

### 🎨 Visual Structure Drawing
- Web-based molecular editor using JSME
- Support for drawing molecules and reactions
- Interactive drawing tools with templates
- Real-time SMILES generation

### 🔗 Seamless Integration
- "Draw Structure" button in the main GUI
- Automatic SMILES transfer to the prediction interface
- No manual copy-pasting required
- Thread-safe communication

### 🌐 Enhanced Web Interface
- Bootstrap-styled interface
- Clear instructions and status feedback
- Keyboard shortcuts for efficiency
- Error handling and validation

## How to Use

1. **Open the Main GUI**
   ```bash
   python simple_reaction_gui.py
   # or
   python demo_structure_drawing.py
   ```

2. **Access the Structure Drawer**
   - Look for the green "Draw Structure" button next to the SMILES input field
   - Click the button to open the molecular editor

3. **Draw Your Structure**
   - Use the web-based editor to draw molecules or reactions
   - Click "Get SMILES" to preview the SMILES string
   - Click "Save to Python" to transfer the structure

4. **Continue with Analysis**
   - Close the web browser
   - The SMILES will appear in the main GUI
   - Proceed with reaction condition prediction

## Technical Implementation

### Components
- **SmilesDrawingServer**: HTTP server handling GUI ↔ Web communication
- **SmilesDrawingDialog**: PyQt6 dialog managing the drawing interface
- **Enhanced HTML Interface**: Improved user experience with Bootstrap styling
- **JavaScript Integration**: JSME editor with custom save functionality

### Architecture
```
PyQt6 GUI ←→ HTTP Server ←→ Web Browser (JSME Editor)
     ↑                                        ↓
     └── SMILES Transfer ←── "Save to Python" ←┘
```

### File Structure
```
SMILES-Drawer/
├── index.html              # Original JSME interface
├── index_enhanced.html     # Enhanced interface with better UX
├── js/
│   ├── script.js           # Original JavaScript
│   └── enhanced_script.js  # Enhanced with status feedback
├── css/
│   └── styles.css          # Styling
└── jsme/
    └── jsme.nocache.js     # JSME library
```

## Requirements

### Python Dependencies
- PyQt6
- http.server (built-in)
- threading (built-in)
- webbrowser (built-in)

### Browser Requirements
- Any modern web browser
- JavaScript enabled
- Preferably Chrome for app-mode support

### JSME Library
- Included in the SMILES-Drawer/jsme/ directory
- No additional installation required

## Troubleshooting

### Common Issues

1. **"Structure editor failed to launch"**
   - Ensure SMILES-Drawer directory exists
   - Check that jsme/jsme.nocache.js is present
   - Verify browser accessibility

2. **"Could not connect to Python server"**
   - Restart the application
   - Check if port 8001 is available
   - Ensure firewall allows local connections

3. **Browser doesn't open**
   - Try manually opening http://localhost:8001
   - Check if Chrome is installed for app-mode
   - Fallback browsers should work

### Debug Mode
Enable console logging by opening browser developer tools (F12) to see detailed error messages.

## Example Usage

```python
# Basic integration example
from simple_reaction_gui import SmilesDrawingDialog

# Open structure drawing dialog
dialog = SmilesDrawingDialog(parent_window)
if dialog.exec() == dialog.DialogCode.Accepted:
    smiles = dialog.get_smiles()
    print(f"Drawn SMILES: {smiles}")
```

## Benefits

1. **User-Friendly**: Visual drawing is more intuitive than typing SMILES
2. **Error Reduction**: Reduces syntax errors in SMILES input
3. **Educational**: Helps users learn molecular structure representation
4. **Efficient**: Faster than manual SMILES construction for complex molecules
5. **Flexible**: Supports both molecules and reaction drawings

## Future Enhancements

- [ ] Support for loading existing SMILES into the editor
- [ ] Export to additional chemical formats (MOL, SDF)
- [ ] Integration with chemical databases for structure lookup
- [ ] Reaction template library
- [ ] Structure validation and cleanup tools
