"""
Simple Reaction Condition Predictor GUI - Basic Edition
======================================================

This GUI provides a basic interface for reaction condition prediction in organic chemistry:
- Simple SMILES input interface
- Reaction type selection
- Basic dark-themed UI
- Minimal functionality for initial testing

This is a simplified version focused on core input functionality only.
"""

from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QGridLayout, QLabel, QLineEdit, QPushButton, QTextEdit, QComboBox,
    QGroupBox, QMessageBox, QFrame, QDialog, QListWidget, QListWidgetItem,
    QSplitter, QCheckBox, QScrollArea, QProgressBar, QSizePolicy
)
from PyQt6.QtCore import Qt, QThread, pyqtSignal, QUrl, QTimer, QMutex
from PyQt6.QtGui import QFont, QPalette, QColor, QPixmap, QImage
import csv
import io
import base64
import http.server
import socketserver
import threading
import time
import json
import tempfile
import shutil
import webbrowser
import subprocess
import os
import mimetypes

# Import reaction types and sample reactions for dropdown
from reaction_types import get_reaction_types
from sample_reactions import (
    get_sample_reactions, get_coupling_reactions, get_reduction_reactions,
    get_oxidation_reactions, get_substitution_reactions, get_elimination_reactions,
    get_cycloaddition_reactions, get_buchwald_hartwig_reactions, search_reactions
)

def create_reaction_image(reaction_smiles: str, width: int = 600, height: int = 200) -> QPixmap:
    """Create reaction image from SMILES using RDKit"""
    # Avoid QPixmap/QPainter before QApplication exists
    try:
        from PyQt6.QtWidgets import QApplication  # local import to avoid circulars
        if QApplication.instance() is None:
            return QPixmap()  # null pixmap; caller should handle gracefully
    except Exception:
        return QPixmap()
    try:
        from rdkit import Chem
        from rdkit.Chem import rdDepictor, AllChem
        
        print(f"Creating reaction image for: '{reaction_smiles}'")
        
        # Parse reaction SMILES
        if ">>" not in reaction_smiles:
            print("No '>>' found in reaction SMILES")
            return create_placeholder_image(reaction_smiles, width, height)
        
        reactants_smiles, products_smiles = reaction_smiles.split(">>", 1)
        print(f"Reactants SMILES: '{reactants_smiles}'")
        print(f"Products SMILES: '{products_smiles}'")
        
        # Create molecules for reactants
        reactants = []
        for smi in reactants_smiles.split('.'):
            smi_clean = smi.strip()
            if smi_clean:  # Skip empty SMILES
                mol = Chem.MolFromSmiles(smi_clean)
                if mol is not None:
                    try:
                        rdDepictor.Compute2DCoords(mol)
                        reactants.append(mol)
                        print(f"Successfully parsed reactant: '{smi_clean}'")
                    except Exception as e:
                        print(f"Error computing 2D coords for reactant '{smi_clean}': {e}")
                else:
                    print(f"Could not parse reactant SMILES: '{smi_clean}'")
        
        # Create molecules for products
        products = []
        for smi in products_smiles.split('.'):
            smi_clean = smi.strip()
            if smi_clean:  # Skip empty SMILES
                mol = Chem.MolFromSmiles(smi_clean)
                if mol is not None:
                    try:
                        rdDepictor.Compute2DCoords(mol)
                        products.append(mol)
                        print(f"Successfully parsed product: '{smi_clean}'")
                    except Exception as e:
                        print(f"Error computing 2D coords for product '{smi_clean}': {e}")
                else:
                    print(f"Could not parse product SMILES: '{smi_clean}'")
        
        print(f"Successfully parsed {len(reactants)} reactants and {len(products)} products")
        
        if not reactants or not products:
            print(f"Failed to parse molecules. Reactants: {len(reactants)}, Products: {len(products)}")
            return create_placeholder_image(reaction_smiles, width, height)
        
        print(f"Successfully parsed {len(reactants)} reactants and {len(products)} products")
        
        # Create a complete reaction image showing all reactants and products
        return create_complete_reaction_image(reactants, products, width, height)
        
    except ImportError:
        # RDKit not available, create a placeholder
        return create_placeholder_image(reaction_smiles, width, height)
    except Exception as e:
        print(f"Error creating reaction image: {e}")
        return create_placeholder_image(reaction_smiles, width, height)

def create_complete_reaction_image(reactants, products, width, height, base_pixmap=None):
    """Create a complete reaction image with arrow showing all reactants and products"""
    # Avoid QPixmap/QPainter before QApplication exists
    try:
        from PyQt6.QtWidgets import QApplication
        if QApplication.instance() is None:
            return QPixmap()
    except Exception:
        return QPixmap()

    from PyQt6.QtGui import QPainter, QPen, QFont
    from PyQt6.QtCore import Qt

    # Use QImage for offscreen painting; convert to QPixmap at end
    image = QImage(width, height, QImage.Format.Format_ARGB32)
    image.fill(Qt.GlobalColor.white)
    painter = QPainter(image)
    painter.setRenderHint(QPainter.RenderHint.Antialiasing)

    # Layout
    arrow_width = 60
    mol_width = (width - arrow_width - 40) // 2
    reactant_x = 20
    reactant_y = 10

    try:
        from rdkit.Chem.Draw import rdMolDraw2D  # type: ignore
        rdkit_available = True
    except Exception:
        rdkit_available = False

    # Draw reactants
    try:
        if reactants and rdkit_available:
            if len(reactants) == 1:
                drawer = rdMolDraw2D.MolDraw2DCairo(mol_width, height - 20)
                drawer.DrawMolecule(reactants[0])
                drawer.FinishDrawing()
                img_data = drawer.GetDrawingText()
                mol_pixmap = QPixmap()
                if mol_pixmap.loadFromData(img_data):
                    painter.drawPixmap(reactant_x, reactant_y, mol_pixmap)
            else:
                reactant_width = mol_width // len(reactants) - 10
                plus_width = 15
                current_x = reactant_x
                for i, mol in enumerate(reactants):
                    actual_width = max(reactant_width, 80)
                    drawer = rdMolDraw2D.MolDraw2DCairo(actual_width, height - 40)
                    drawer.DrawMolecule(mol)
                    drawer.FinishDrawing()
                    img_data = drawer.GetDrawingText()
                    mol_pixmap = QPixmap()
                    if mol_pixmap.loadFromData(img_data):
                        y_centered = reactant_y + (height - 40 - mol_pixmap.height()) // 2
                        painter.drawPixmap(current_x, y_centered, mol_pixmap)
                    current_x += actual_width
                    if i < len(reactants) - 1:
                        painter.setPen(QPen(Qt.GlobalColor.black, 2))
                        painter.setFont(QFont("Arial", 12, QFont.Weight.Bold))
                        plus_y = height // 2
                        painter.drawText(current_x + 2, plus_y + 5, "+")
                        current_x += plus_width
        else:
            painter.setPen(QPen(Qt.GlobalColor.black))
            painter.setFont(QFont("Arial", 10))
            painter.drawText(reactant_x, height//2, f"{len(reactants)} Reactants")
    except Exception as e:
        print(f"Error drawing reactants: {e}")

    # Arrow
    arrow_x = reactant_x + mol_width + 10
    arrow_y = height // 2
    painter.setPen(QPen(Qt.GlobalColor.black, 3))
    painter.drawLine(arrow_x, arrow_y, arrow_x + arrow_width - 15, arrow_y)
    painter.drawLine(arrow_x + arrow_width - 15, arrow_y, arrow_x + arrow_width - 25, arrow_y - 8)
    painter.drawLine(arrow_x + arrow_width - 15, arrow_y, arrow_x + arrow_width - 25, arrow_y + 8)

    # Products
    product_x = arrow_x + arrow_width + 10
    try:
        if products and rdkit_available:
            if len(products) == 1:
                drawer = rdMolDraw2D.MolDraw2DCairo(mol_width, height - 20)
                drawer.DrawMolecule(products[0])
                drawer.FinishDrawing()
                img_data = drawer.GetDrawingText()
                mol_pixmap = QPixmap()
                if mol_pixmap.loadFromData(img_data):
                    painter.drawPixmap(product_x, reactant_y, mol_pixmap)
            else:
                product_width = mol_width // len(products) - 10
                plus_width = 15
                current_x = product_x
                for i, mol in enumerate(products):
                    actual_width = max(product_width, 80)
                    drawer = rdMolDraw2D.MolDraw2DCairo(actual_width, height - 40)
                    drawer.DrawMolecule(mol)
                    drawer.FinishDrawing()
                    img_data = drawer.GetDrawingText()
                    mol_pixmap = QPixmap()
                    if mol_pixmap.loadFromData(img_data):
                        y_centered = reactant_y + (height - 40 - mol_pixmap.height()) // 2
                        painter.drawPixmap(current_x, y_centered, mol_pixmap)
                    current_x += actual_width
                    if i < len(products) - 1:
                        painter.setPen(QPen(Qt.GlobalColor.black, 2))
                        painter.setFont(QFont("Arial", 12, QFont.Weight.Bold))
                        plus_y = height // 2
                        painter.drawText(current_x + 2, plus_y + 5, "+")
                        current_x += plus_width
        else:
            painter.setPen(QPen(Qt.GlobalColor.black))
            painter.setFont(QFont("Arial", 10))
            painter.drawText(product_x, height//2, f"{len(products)} Products")
    except Exception as e:
        print(f"Error drawing products: {e}")

    painter.end()
    return QPixmap.fromImage(image)

def create_placeholder_image(reaction_smiles: str, width: int = 600, height: int = 200) -> QPixmap:
    """Create a placeholder image when RDKit is not available"""
    # Avoid QPixmap/QPainter before QApplication exists
    try:
        from PyQt6.QtWidgets import QApplication
        if QApplication.instance() is None:
            return QPixmap()
    except Exception:
        return QPixmap()
    from PyQt6.QtGui import QPainter, QPen, QBrush, QFont
    from PyQt6.QtCore import Qt
    
    image = QImage(width, height, QImage.Format.Format_ARGB32)
    image.fill(Qt.GlobalColor.white)
    
    painter = QPainter(image)
    painter.setPen(QPen(Qt.GlobalColor.black, 2))
    
    # Draw border
    painter.drawRect(5, 5, width - 10, height - 10)
    
    # Draw arrow
    arrow_y = height // 2
    arrow_start = width // 3
    arrow_end = 2 * width // 3
    
    painter.drawLine(arrow_start, arrow_y, arrow_end, arrow_y)
    painter.drawLine(arrow_end - 10, arrow_y - 10, arrow_end, arrow_y)
    painter.drawLine(arrow_end - 10, arrow_y + 10, arrow_end, arrow_y)
    
    # Add text
    font = QFont("Arial", 10)
    painter.setFont(font)
    
    # Draw "Reactants" and "Products"
    painter.drawText(20, arrow_y - 20, "Reactants")
    painter.drawText(arrow_end + 20, arrow_y - 20, "Products")
    
    # Draw SMILES (truncated if too long)
    if ">>" in reaction_smiles:
        reactants, products = reaction_smiles.split(">>", 1)
        reactants = reactants[:30] + "..." if len(reactants) > 30 else reactants
        products = products[:30] + "..." if len(products) > 30 else products
        
        small_font = QFont("Arial", 8)
        painter.setFont(small_font)
        painter.drawText(20, arrow_y + 40, reactants)
        painter.drawText(arrow_end + 20, arrow_y + 40, products)
    
    # Add note about RDKit
    painter.setFont(QFont("Arial", 8))
    painter.setPen(QPen(Qt.GlobalColor.gray))
    painter.drawText(20, height - 20, "Install RDKit for structure visualization")
    
    painter.end()
    return QPixmap.fromImage(image)

class SimplePredictionWorker(QThread):
    """Enhanced prediction worker with recommendation engine"""
    
    finished = pyqtSignal(dict)
    error = pyqtSignal(str)
    progress = pyqtSignal(int)
    
    def __init__(self, reaction_smiles: str, reaction_type: str):
        super().__init__()
        self.reaction_smiles = reaction_smiles
        self.reaction_type = reaction_type
    
    def run(self):
        """Run prediction with enhanced recommendation engine"""
        try:
            self.progress.emit(20)
            
            # Basic validation
            if not self.reaction_smiles.strip():
                raise ValueError("SMILES input is required")
            
            self.progress.emit(40)
            
            # Try to use enhanced recommendation engine first
            try:
                from enhanced_recommendation_engine import EnhancedRecommendationEngine
                engine = EnhancedRecommendationEngine()
                self.progress.emit(60)
                
                # Get enhanced recommendations
                recommendations = engine.get_recommendations(
                    self.reaction_smiles, 
                    self.reaction_type
                )
                
                self.progress.emit(80)
                
                result = {
                    'reaction_smiles': self.reaction_smiles,
                    'reaction_type': self.reaction_type,
                    'status': 'Enhanced recommendations generated successfully',
                    'recommendations': recommendations,
                    'available_recommenders': engine.get_available_recommenders(),
                    'analysis_type': recommendations.get('analysis_type', 'enhanced')
                }
                
            except ImportError as e:
                # Fallback to original recommendation engine
                try:
                    from recommendation_engine import RecommendationEngine
                    engine = RecommendationEngine()
                    self.progress.emit(60)
                    
                    # Get original recommendations
                    recommendations = engine.get_recommendations(
                        self.reaction_smiles, 
                        self.reaction_type
                    )
                    
                    self.progress.emit(80)
                    
                    result = {
                        'reaction_smiles': self.reaction_smiles,
                        'reaction_type': self.reaction_type,
                        'status': 'Recommendations generated successfully',
                        'recommendations': recommendations,
                        'available_recommenders': engine.get_available_recommenders(),
                        'analysis_type': recommendations.get('analysis_type', 'general')
                    }
                    
                except Exception as rec_error:
                    # Fallback to basic analysis
                    result = {
                        'reaction_smiles': self.reaction_smiles,
                        'selected_reaction_type': self.reaction_type,
                        'status': 'Basic analysis completed (recommendation engines failed)',
                        'message': f'Recommendation error: {str(rec_error)}. Enhanced ligand/solvent systems may not be installed.',
                        'analysis_type': 'basic'
                    }
            except Exception as rec_error:
                # Error with enhanced engine
                result = {
                    'reaction_smiles': self.reaction_smiles,
                    'selected_reaction_type': self.reaction_type,
                    'status': 'Analysis completed with limitations',
                    'message': f'Enhanced recommendation error: {str(rec_error)}',
                    'analysis_type': 'basic'
                }
            
            self.progress.emit(100)
            self.finished.emit(result)
            
        except Exception as e:
            self.error.emit(f"Prediction Error: {str(e)}")

class SampleReactionsBrowser(QDialog):
    """Dialog for browsing and selecting from sample reactions"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.selected_reaction = None
        self.init_ui()
        self.load_reactions()
    
    def init_ui(self):
        """Initialize the sample reactions browser UI"""
        self.setWindowTitle("Sample Reactions Browser")
        # Removed fixed geometry - will be set dynamically by parent
        self.setModal(True)
        
        # Ensure the dialog has a proper title bar and is moveable
        self.setWindowFlags(Qt.WindowType.Dialog | Qt.WindowType.WindowTitleHint | Qt.WindowType.WindowCloseButtonHint | Qt.WindowType.WindowMinMaxButtonsHint)
        
        # Set minimum size to ensure usability
        self.setMinimumSize(1200, 700)
        
        # Main layout
        layout = QVBoxLayout(self)
        layout.setSpacing(10)
        layout.setContentsMargins(15, 15, 15, 15)
        
        # Filter section
        filter_group = QGroupBox("Filter Reactions")
        filter_layout = QVBoxLayout(filter_group)
        
        # Search box
        search_layout = QHBoxLayout()
        search_label = QLabel("Search:")
        search_label.setMinimumWidth(60)
        self.search_input = QLineEdit()
        self.search_input.setPlaceholderText("Search by reaction type, reagents, or description...")
        self.search_input.textChanged.connect(self.filter_reactions)
        search_layout.addWidget(search_label)
        search_layout.addWidget(self.search_input)
        filter_layout.addLayout(search_layout)
        
        # Category filters (Coupling-focused with subcategories)
        category_layout = QVBoxLayout()
        self.category_checkboxes = {}

        # Top-level quick toggles
        top_row = QHBoxLayout()
        for cat_name, cat_key in [("All", "all"), ("Coupling", "coupling_only"), ("Everything Else", "non_coupling")]:
            cb = QCheckBox(cat_name)
            # Default: "All" should be unchecked on start
            # (no explicit setChecked(True) here)
            cb.stateChanged.connect(self.filter_by_category)
            self.category_checkboxes[cat_key] = cb
            top_row.addWidget(cb)
        top_row.addStretch()
        category_layout.addLayout(top_row)

        # Coupling families
        family_row = QHBoxLayout()
        for cat_name, cat_key in [("C-C", "cc"), ("C-N", "cn"), ("C-O", "co"), ("C-S", "cs")]:
            cb = QCheckBox(cat_name)
            cb.stateChanged.connect(self.filter_by_category)
            self.category_checkboxes[cat_key] = cb
            family_row.addWidget(cb)
        family_row.addStretch()
        category_layout.addLayout(family_row)

        # Sub-categories per family (wrap across rows to avoid cutoff)
        sub_layout = QHBoxLayout()

        # Helper to populate a group box with a grid of checkboxes
        def add_checkboxes_grid(group_title, items):
            box = QGroupBox(group_title)
            grid = QGridLayout(box)
            cols = 2  # arrange in two columns to reduce horizontal overflow
            for idx, (name, key) in enumerate(items):
                cb = QCheckBox(name)
                cb.stateChanged.connect(self.filter_by_category)
                self.category_checkboxes[key] = cb
                row = idx // cols
                col = idx % cols
                grid.addWidget(cb, row, col)
            return box

        # C-C subcats
        cc_items = [
            ("Suzuki (Pd)", "cc_suzuki"),
            ("Stille (Pd)", "cc_stille"),
            ("Sonogashira (Pd)", "cc_sonogashira"),
            ("Heck (Pd)", "cc_heck"),
            ("Negishi (Pd)", "cc_negishi"),
            ("Kumada (Ni)", "cc_kumada"),
        ]
        sub_layout.addWidget(add_checkboxes_grid("C-C Couplings", cc_items))

        # C-N subcats
        cn_items = [
            ("Buchwald-Hartwig (Pd)", "cn_bh"),
            ("Ullmann (Cu)", "cn_ullmann"),
            ("Chan-Lam (Cu)", "cn_chanlam"),
        ]
        sub_layout.addWidget(add_checkboxes_grid("C-N Couplings", cn_items))

        # C-O subcats
        co_items = [
            ("Ullmann Ether (Cu)", "co_ullmann_ether"),
            ("Mitsunobu", "co_mitsunobu"),
        ]
        sub_layout.addWidget(add_checkboxes_grid("C-O Couplings", co_items))

        # C-S subcats
        cs_items = [
            ("Thioether Coupling (Pd)", "cs_thioether"),
        ]
        sub_layout.addWidget(add_checkboxes_grid("C-S Couplings", cs_items))

        category_layout.addLayout(sub_layout)

        filter_layout.addLayout(category_layout)
        layout.addWidget(filter_group)
        
        # Reactions list and details splitter
        content_splitter = QSplitter(Qt.Orientation.Horizontal)

        # Reactions list
        list_group = QGroupBox("Reactions List")
        list_layout = QVBoxLayout(list_group)

        self.reactions_list = QListWidget()
        self.reactions_list.setMinimumWidth(500)  # Increased from 400 to 500
        self.reactions_list.itemClicked.connect(self.on_reaction_selected)
        self.reactions_list.itemDoubleClicked.connect(self.on_reaction_double_clicked)
        list_layout.addWidget(self.reactions_list)

        # Stats label
        self.stats_label = QLabel("Total reactions: 0")
        self.stats_label.setStyleSheet("color: #CCCCCC; font-style: italic;")
        list_layout.addWidget(self.stats_label)

        content_splitter.addWidget(list_group)

        # Reaction details
        details_group = QGroupBox("Reaction Details")
        details_layout = QVBoxLayout(details_group)

        # Reaction image display
        self.details_image_label = QLabel()
        self.details_image_label.setMinimumHeight(200)
        self.details_image_label.setMaximumHeight(280)
        self.details_image_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.details_image_label.setScaledContents(True)
        self.details_image_label.setStyleSheet("""
            QLabel {
                background-color: #404040;
                border: 1px solid #555555;
                border-radius: 5px;
                color: #CCCCCC;
                font-style: italic;
            }
        """)
        self.details_image_label.setText("Select a reaction to view the structure")
        details_layout.addWidget(self.details_image_label)

        # Text details
        self.details_text = QTextEdit()
        self.details_text.setMinimumWidth(450)
        self.details_text.setMinimumHeight(200)
        self.details_text.setReadOnly(True)
        self.details_text.setPlaceholderText("Select a reaction to view details...")
        details_layout.addWidget(self.details_text)

        content_splitter.addWidget(details_group)
        content_splitter.setSizes([500, 600])

        layout.addWidget(content_splitter)

        # Buttons
        button_layout = QHBoxLayout()

        self.select_button = QPushButton("Select Reaction")
        self.select_button.setMinimumHeight(35)
        self.select_button.clicked.connect(self.select_reaction)
        self.select_button.setEnabled(False)

        cancel_button = QPushButton("Cancel")
        cancel_button.setMinimumHeight(35)
        cancel_button.clicked.connect(self.reject)

        button_layout.addStretch()
        button_layout.addWidget(self.select_button)
        button_layout.addWidget(cancel_button)

        layout.addLayout(button_layout)
        
        # Apply dark theme
        self.apply_dark_theme()
    
    def apply_dark_theme(self):
        """Apply dark theme to the dialog"""
        self.setStyleSheet("""
            QDialog {
                background-color: #2b2b2b;
                color: #ffffff;
            }
            QWidget {
                background-color: #2b2b2b;
                color: #ffffff;
            }
            QGroupBox {
                color: #ffffff;
                background-color: #3c3c3c;
                border: 2px solid #555555;
                border-radius: 5px;
                margin-top: 10px;
                padding-top: 10px;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                left: 10px;
                padding: 0 5px 0 5px;
            }
            QLabel {
                color: #ffffff;
            }
            QLineEdit {
                background-color: #404040;
                border: 1px solid #555555;
                color: #ffffff;
                padding: 5px;
                border-radius: 3px;
            }
            QLineEdit:focus {
                border: 2px solid #4CAF50;
            }
            QListWidget {
                background-color: #404040;
                border: 1px solid #555555;
                color: #ffffff;
                border-radius: 3px;
            }
            QListWidget::item {
                padding: 5px;
                border-bottom: 1px solid #555555;
            }
            QListWidget::item:selected {
                background-color: #4CAF50;
                color: white;
            }
            QListWidget::item:hover {
                background-color: #555555;
            }
            QTextEdit {
                background-color: #404040;
                border: 1px solid #555555;
                color: #ffffff;
                border-radius: 3px;
            }
            QCheckBox {
                color: #ffffff;
            }
            QCheckBox::indicator {
                width: 15px;
                height: 15px;
            }
            QCheckBox::indicator:unchecked {
                background-color: #404040;
                border: 1px solid #555555;
            }
            QCheckBox::indicator:checked {
                background-color: #4CAF50;
                border: 1px solid #4CAF50;
            }
            QPushButton {
                background-color: #404040;
                border: 1px solid #555555;
                color: #ffffff;
                padding: 8px;
                border-radius: 3px;
            }
            QPushButton:hover {
                background-color: #4CAF50;
            }
            QPushButton:pressed {
                background-color: #45a049;
            }
            QPushButton:disabled {
                background-color: #2b2b2b;
                color: #777777;
            }
        """)
    
    def load_reactions(self):
        """Load all sample reactions"""
        try:
            self.all_reactions = get_sample_reactions()
            self.filtered_reactions = self.all_reactions.copy()
            self.populate_list()
        except Exception as e:
            QMessageBox.warning(self, "Error", f"Could not load sample reactions: {e}")
            self.all_reactions = []
            self.filtered_reactions = []
    
    def populate_list(self):
        """Populate the reactions list widget"""
        self.reactions_list.clear()
        
        for reaction in self.filtered_reactions:
            if reaction.strip() and not reaction.startswith("Select a sample"):
                item = QListWidgetItem(reaction)
                # Make the display more readable
                display_text = self.format_reaction_display(reaction)
                item.setText(display_text)
                item.setData(Qt.ItemDataRole.UserRole, reaction)  # Store original data
                self.reactions_list.addItem(item)
        
        # Update stats
        count = self.reactions_list.count()
        self.stats_label.setText(f"Total reactions: {count}")
    
    def format_reaction_display(self, reaction):
        """Format reaction for better display in list"""
        if ">>" in reaction:
            parts = reaction.split(">>")
            if len(parts) == 2:
                reactants = parts[0].strip()
                if "(" in parts[1]:
                    products_and_desc = parts[1].split("(", 1)
                    products = products_and_desc[0].strip()
                    description = "(" + products_and_desc[1] if len(products_and_desc) > 1 else ""
                    
                    # For display purposes, show truncated SMILES but keep full SMILES for processing
                    display_reactants = reactants[:30] + "..." if len(reactants) > 30 else reactants
                    display_products = products[:30] + "..." if len(products) > 30 else products
                    
                    return f"{display_reactants} → {display_products} {description}"
        
        return reaction
    
    def filter_reactions(self):
        """Filter reactions based on search text"""
        search_text = self.search_input.text().lower()
        
        if not search_text:
            self.apply_category_filter()
            return
        
        # Get base reactions from category filter
        base_reactions = self.get_category_filtered_reactions()
        
        # Apply text search
        self.filtered_reactions = [
            reaction for reaction in base_reactions
            if search_text in reaction.lower()
        ]
        
        self.populate_list()
    
    def filter_by_category(self):
        """Filter reactions by selected categories"""
        # Handle "All" checkbox
        if self.category_checkboxes["all"].isChecked():
            # Uncheck other categories when "All" is checked
            for key, checkbox in self.category_checkboxes.items():
                if key != "all":
                    checkbox.setChecked(False)
        else:
            # If any specific category is checked, uncheck "All"
            any_specific_checked = any(
                checkbox.isChecked() for key, checkbox in self.category_checkboxes.items() 
                if key != "all"
            )
            if any_specific_checked:
                self.category_checkboxes["all"].setChecked(False)
        
        self.apply_category_filter()
    
    def apply_category_filter(self):
        """Apply category filtering"""
        if self.category_checkboxes["all"].isChecked():
            self.filtered_reactions = self.all_reactions.copy()
        else:
            self.filtered_reactions = self.get_category_filtered_reactions()
        
        # Reapply text search if any
        search_text = self.search_input.text().lower()
        if search_text:
            self.filtered_reactions = [
                reaction for reaction in self.filtered_reactions
                if search_text in reaction.lower()
            ]
        
        self.populate_list()
    
    def get_category_filtered_reactions(self):
        """Get reactions filtered by selected categories"""
        filtered = []
        
        try:
            # Quick toggles
            coupling_only = self.category_checkboxes["coupling_only"].isChecked()
            non_coupling = self.category_checkboxes["non_coupling"].isChecked()

            # Family selections
            cc = self.category_checkboxes["cc"].isChecked()
            cn = self.category_checkboxes["cn"].isChecked()
            co = self.category_checkboxes["co"].isChecked()
            cs = self.category_checkboxes["cs"].isChecked()

            # Subcategory selections
            cc_sub = [
                ("Suzuki", self.category_checkboxes["cc_suzuki"].isChecked()),
                ("Stille", self.category_checkboxes["cc_stille"].isChecked()),
                ("Sonogashira", self.category_checkboxes["cc_sonogashira"].isChecked()),
                ("Heck", self.category_checkboxes["cc_heck"].isChecked()),
                ("Negishi", self.category_checkboxes["cc_negishi"].isChecked()),
                ("Kumada", self.category_checkboxes["cc_kumada"].isChecked()),
            ]

            cn_sub = [
                ("Buchwald-Hartwig", self.category_checkboxes["cn_bh"].isChecked()),
                ("Ullmann C-N", self.category_checkboxes["cn_ullmann"].isChecked()),
                ("Chan-Lam", self.category_checkboxes["cn_chanlam"].isChecked()),
            ]

            co_sub = [
                ("Ullmann Ether", self.category_checkboxes["co_ullmann_ether"].isChecked()),
                ("Mitsunobu", self.category_checkboxes["co_mitsunobu"].isChecked()),
            ]

            cs_sub = [
                ("C-S Coupling", self.category_checkboxes["cs_thioether"].isChecked()),
                ("Thioether Formation", self.category_checkboxes["cs_thioether"].isChecked()),
            ]

            # Helper: add a reaction if it contains any of the tokens
            def add_by_tokens(tokens):
                for r in self.all_reactions:
                    if any(tok in r for tok in tokens):
                        filtered.append(r)

            # Coupling-only toggle
            if coupling_only:
                from sample_reactions import get_coupling_reactions
                filtered.extend(get_coupling_reactions())

            # Family level
            if cc:
                add_by_tokens(["Suzuki", "Stille", "Sonogashira", "Heck", "Negishi", "Kumada"])
            if cn:
                add_by_tokens(["Buchwald-Hartwig", "Ullmann C-N", "Chan-Lam"])
            if co:
                add_by_tokens(["Ullmann Ether", "Mitsunobu"])
            if cs:
                add_by_tokens(["C-S Coupling", "Thioether Formation"])

            # Subcategories
            if any(chk for _, chk in cc_sub):
                add_by_tokens([name for name, chk in cc_sub if chk])
            if any(chk for _, chk in cn_sub):
                add_by_tokens([name for name, chk in cn_sub if chk])
            if any(chk for _, chk in co_sub):
                add_by_tokens([name for name, chk in co_sub if chk])
            if any(chk for _, chk in cs_sub):
                add_by_tokens([name for name, chk in cs_sub if chk])

            # Non-coupling quick toggle (common other categories)
            if non_coupling:
                from sample_reactions import (
                    get_reduction_reactions, get_oxidation_reactions,
                    get_substitution_reactions, get_elimination_reactions,
                    get_cycloaddition_reactions,
                )
                filtered.extend(get_reduction_reactions())
                filtered.extend(get_oxidation_reactions())
                filtered.extend(get_substitution_reactions())
                filtered.extend(get_elimination_reactions())
                filtered.extend(get_cycloaddition_reactions())
        except Exception as e:
            print(f"Error filtering by category: {e}")
            return self.all_reactions.copy()
        
        # Remove duplicates while preserving order
        seen = set()
        unique_filtered = []
        for reaction in filtered:
            if reaction not in seen:
                seen.add(reaction)
                unique_filtered.append(reaction)
        
        return unique_filtered if unique_filtered else self.all_reactions.copy()
    
    def on_reaction_selected(self, item):
        """Handle reaction selection"""
        # Get the original reaction data (not the formatted display text)
        original_reaction = item.data(Qt.ItemDataRole.UserRole)
        self.show_reaction_details(original_reaction)
        self.select_button.setEnabled(True)
        self.selected_reaction = original_reaction
    
    def on_reaction_double_clicked(self, item):
        """Handle reaction double-click (select and close)"""
        self.on_reaction_selected(item)
        self.select_reaction()
    
    def show_reaction_details(self, reaction):
        """Show detailed information about the selected reaction"""
        if not reaction or reaction.startswith("Select a sample"):
            self.details_text.clear()
            self.details_image_label.clear()
            self.details_image_label.setText("Select a reaction to view the structure")
            return
        
        # Parse reaction to extract clean SMILES (before any description)
        # Important: Use the original reaction data, not the truncated display text
        original_reaction = reaction  # This should be the full, untruncated reaction
        smiles_part = ""
        description = ""
        reactants = ""
        products = ""
        
        if ">>" in original_reaction:
            # Split at the first occurrence of '>>' to separate reactants and products
            parts = original_reaction.split(">>", 1)
            reactants = parts[0].strip()
            
            # Now we need to carefully extract the products part before any description
            products_and_rest = parts[1].strip()
            
            # Find the first space that comes before a '(' - this indicates start of description
            if " (" in products_and_rest:
                # Split at the first occurrence of " (" to separate products from description
                product_parts = products_and_rest.split(" (", 1)
                products = product_parts[0].strip()
                description = " (" + product_parts[1] if len(product_parts) > 1 else ""
            else:
                # No description, the whole thing is products
                products = products_and_rest
                description = ""
            
            smiles_part = f"{reactants}>>{products}"
        else:
            # If no '>>' found, try to extract SMILES from beginning
            if " (" in original_reaction:
                smiles_part = original_reaction.split(" (")[0].strip()
            else:
                smiles_part = original_reaction.strip()
        
        # Debug: Print the SMILES being used for image generation
        print(f"Original reaction input: '{original_reaction}'")
        print(f"Extracted SMILES for image: '{smiles_part}'")
        print(f"Reactants: '{reactants}' | Products: '{products}'")
        
        # Update reaction image
        try:
            # Try full image first; on failure, display a placeholder
            if smiles_part:
                pixmap = None
                if ">>" in smiles_part:
                    pixmap = create_reaction_image(smiles_part, 580, 200)
                # If primary draw failed or not a strict reaction string, try placeholder
                if not pixmap or pixmap.isNull():
                    placeholder = create_placeholder_image(smiles_part, 580, 200)
                    if placeholder and not placeholder.isNull():
                        self.details_image_label.setPixmap(placeholder)
                    else:
                        self.details_image_label.clear()
                        self.details_image_label.setText("Could not generate reaction image")
                else:
                    self.details_image_label.setPixmap(pixmap)
            else:
                self.details_image_label.clear()
                self.details_image_label.setText("Invalid reaction format")
        except Exception as e:
            self.details_image_label.clear()
            self.details_image_label.setText(f"Error generating image: {str(e)}")
            print(f"Image generation error for '{smiles_part}': {e}")
        
        # Classify reaction type
        reaction_type = self.classify_reaction_type(original_reaction)
        
        # Count components
        reactant_count = 0
        if ">>" in smiles_part:
            reactants_smiles = smiles_part.split(">>")[0]
            reactant_count = len(reactants_smiles.split('.')) if reactants_smiles else 0
        
        details = f"""Reaction Details:

SMILES: {smiles_part}

Description: {description if description else 'No description available'}

Reaction Type: {reaction_type}

Components:
• Reactants: {reactant_count} component(s)
• Complexity: {self.assess_complexity(smiles_part)}

Notes:
• Double-click to select this reaction
• This will populate the main interface with this example
• You can then run predictions on this reaction
"""
        
        self.details_text.setPlainText(details)
    
    def classify_reaction_type(self, reaction):
        """Classify the reaction type based on description"""
        reaction_lower = reaction.lower()
        
        if "suzuki" in reaction_lower:
            return "C-C Coupling - Suzuki-Miyaura (Pd)"
        elif "stille" in reaction_lower:
            return "C-C Coupling - Stille (Pd)"
        elif "sonogashira" in reaction_lower:
            return "C-C Coupling - Sonogashira (Pd)"
        elif "heck" in reaction_lower:
            return "C-C Coupling - Heck (Pd)"
        elif "negishi" in reaction_lower:
            return "C-C Coupling - Negishi (Pd)"
        elif "buchwald" in reaction_lower or "hartwig" in reaction_lower:
            return "C-N Coupling - Buchwald-Hartwig (Pd)"
        elif "chan-lam" in reaction_lower:
            return "C-N Oxidative Coupling - Chan-Lam (Cu)"
        elif "ullmann ether" in reaction_lower:
            return "C-O Coupling - Ullmann (Cu)"
        elif "ullmann c-n" in reaction_lower or "ullmann" in reaction_lower:
            return "C-N Coupling - Ullmann (Cu)"
        elif "esterification" in reaction_lower:
            return "Esterification"
        elif "amidation" in reaction_lower:
            return "Amidation"
        elif "hydrogenation" in reaction_lower:
            return "Hydrogenation"
        elif "oxidation" in reaction_lower:
            return "Oxidation"
        elif "reduction" in reaction_lower:
            return "Reduction"
        elif "sn1" in reaction_lower:
            return "SN1 Substitution"
        elif "sn2" in reaction_lower:
            return "SN2 Substitution"
        elif "e1" in reaction_lower:
            return "E1 Elimination"
        elif "e2" in reaction_lower:
            return "E2 Elimination"
        elif "diels-alder" in reaction_lower:
            return "Diels-Alder Cycloaddition"
        elif "click" in reaction_lower:
            return "Click Chemistry"
        else:
            return "General Organic Reaction"
    
    def assess_complexity(self, smiles):
        """Assess reaction complexity"""
        if not smiles:
            return "Unknown"
        
        # Simple heuristics
        if len(smiles) < 50:
            return "Simple"
        elif len(smiles) < 100:
            return "Moderate"
        else:
            return "Complex"
    
    def select_reaction(self):
        """Select the current reaction and close dialog"""
        if self.selected_reaction:
            self.accept()
        else:
            QMessageBox.warning(self, "No Selection", "Please select a reaction first.")
    
    def get_selected_reaction(self):
        """Get the selected reaction"""
        return self.selected_reaction

class SmilesDrawingServer:
    """HTTP server for handling SMILES data from the web-based drawing interface"""
    
    def __init__(self, port=8001):
        self.port = port
        self.server = None
        self.server_thread = None
        self.received_smiles = None
        self.smiles_mutex = QMutex()
        
    def start_server(self):
        """Start the HTTP server in a background thread"""
        if self.server_thread and self.server_thread.is_alive():
            return
            
        handler = self.create_request_handler()
        try:
            self.server = socketserver.TCPServer(('', self.port), handler)
            self.server_thread = threading.Thread(target=self.server.serve_forever, daemon=True)
            self.server_thread.start()
            print(f"SMILES drawing server started on port {self.port}")
        except Exception as e:
            print(f"Failed to start SMILES server: {e}")
            raise
    
    def stop_server(self):
        """Stop the HTTP server"""
        if self.server:
            self.server.shutdown()
            self.server.server_close()
            self.server = None
        if self.server_thread:
            self.server_thread.join(timeout=1.0)
            self.server_thread = None
    
    def create_request_handler(self):
        """Create a request handler class with access to this server instance"""
        server_instance = self
        
        class SmilesRequestHandler(http.server.SimpleHTTPRequestHandler):
            def do_POST(self):
                if self.path == '/save_smiles':
                    try:
                        content_length = int(self.headers['Content-Length'])
                        post_data = self.rfile.read(content_length)
                        smiles_string = post_data.decode('utf-8').strip()
                        
                        # Store the received SMILES in a thread-safe way
                        server_instance.smiles_mutex.lock()
                        server_instance.received_smiles = smiles_string
                        server_instance.smiles_mutex.unlock()
                        
                        print(f"[SMILES Received]: {smiles_string}")
                        
                        self.send_response(200)
                        self.send_header('Content-type', 'text/plain')
                        self.send_header('Access-Control-Allow-Origin', '*')
                        self.end_headers()
                        self.wfile.write(b"OK")
                    except Exception as e:
                        print(f"Error processing POST request: {e}")
                        self.send_response(500)
                        self.end_headers()
                        self.wfile.write(b"ERROR")
                else:
                    self.send_response(404)
                    self.end_headers()
                    self.wfile.write(b"Not Found")
            
            def do_GET(self):
                # Serve files from the SMILES-Drawer directory
                try:
                    smiles_drawer_path = os.path.join(os.path.dirname(__file__), 'SMILES-Drawer')
                    
                    # Map the request path to the actual file
                    requested_path = self.path
                    
                    # Override the default index.html to use our enhanced version if available
                    if requested_path == '/' or requested_path == '/index.html':
                        enhanced_path = os.path.join(smiles_drawer_path, 'index_enhanced.html')
                        if os.path.exists(enhanced_path):
                            requested_path = '/index_enhanced.html'
                        else:
                            requested_path = '/index.html'
                    
                    # Handle compact interface request
                    if requested_path == '/compact.html':
                        compact_path = os.path.join(smiles_drawer_path, 'compact.html')
                        if os.path.exists(compact_path):
                            requested_path = '/compact.html'
                        else:
                            # Fallback to enhanced or regular interface
                            enhanced_path = os.path.join(smiles_drawer_path, 'index_enhanced.html')
                            if os.path.exists(enhanced_path):
                                requested_path = '/index_enhanced.html'
                            else:
                                requested_path = '/index.html'
                    
                    # Check for enhanced script
                    if requested_path == '/js/script.js':
                        enhanced_script_path = os.path.join(smiles_drawer_path, 'js', 'enhanced_script.js')
                        if os.path.exists(enhanced_script_path):
                            requested_path = '/js/enhanced_script.js'
                    
                    # Remove leading slash and construct full file path
                    relative_path = requested_path.lstrip('/')
                    full_path = os.path.join(smiles_drawer_path, relative_path)
                    
                    # Security check - ensure the path is within the SMILES-Drawer directory
                    if not os.path.commonpath([smiles_drawer_path, full_path]) == smiles_drawer_path:
                        self.send_error(403, "Forbidden")
                        return
                    
                    # Check if file exists
                    if not os.path.exists(full_path):
                        self.send_error(404, "File not found")
                        return
                    
                    # Serve the file
                    self.serve_file(full_path)
                    
                except Exception as e:
                    print(f"Error serving file: {e}")
                    self.send_error(500, f"Internal server error: {str(e)}")
            
            def serve_file(self, file_path):
                """Serve a specific file"""
                try:
                    # Determine content type
                    content_type = self.guess_type(file_path)
                    
                    # Read and serve the file
                    with open(file_path, 'rb') as f:
                        content = f.read()
                    
                    self.send_response(200)
                    self.send_header('Content-type', content_type)
                    self.send_header('Content-length', str(len(content)))
                    self.end_headers()
                    self.wfile.write(content)
                    
                except Exception as e:
                    print(f"Error reading file {file_path}: {e}")
                    self.send_error(500, f"Error reading file: {str(e)}")
            
            def guess_type(self, path):
                """Guess the content type based on file extension"""
                import mimetypes
                content_type, _ = mimetypes.guess_type(path)
                return content_type or 'application/octet-stream'
                    
            def log_message(self, format, *args):
                # Suppress server log messages to reduce console clutter
                pass
                
        return SmilesRequestHandler
    
    def get_smiles(self):
        """Get the most recently received SMILES string"""
        self.smiles_mutex.lock()
        smiles = self.received_smiles
        self.received_smiles = None  # Clear after reading
        self.smiles_mutex.unlock()
        return smiles

class SmilesDrawingDialog(QDialog):
    """Dialog for drawing chemical structures using the web-based JSME editor"""
    
    def __init__(self, parent=None, initial_smiles=""):
        super().__init__(parent)
        self.drawn_smiles = ""
        self.initial_smiles = initial_smiles
        self.server = None
        self.check_timer = None
        self.browser_process = None
        self.init_ui()
        
    def init_ui(self):
        """Initialize the dialog UI"""
        self.setWindowTitle("Draw Chemical Structure")
        self.setModal(True)
        self.resize(500, 300)

        layout = QVBoxLayout(self)

        # Instructions
        instructions = QLabel(
            """
<h3>Chemical Structure Drawing</h3>
<p>This opens a web-based molecular editor where you can:</p>
<ul>
<li>Draw molecules and reactions using the JSME editor</li>
<li>Use the "Get SMILES" button to see the SMILES string</li>
<li>Use the "Save to Python" button to send the structure back to this application</li>
<li>Close the browser window when you're finished</li>
</ul>
<p><b>Note:</b> Be sure to click "Save to Python" in the web editor before closing.</p>
"""
        )
        instructions.setWordWrap(True)
        instructions.setStyleSheet(
            "QLabel { color: #ffffff; background-color: #3c3c3c; padding: 15px; border-radius: 5px; }"
        )
        layout.addWidget(instructions)

        # Status label
        self.status_label = QLabel("Ready to launch the structure editor...")
        self.status_label.setStyleSheet("color: #4CAF50; font-weight: bold; padding: 10px;")
        layout.addWidget(self.status_label)

        # Progress bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        layout.addWidget(self.progress_bar)

        # Current SMILES display
        smiles_layout = QVBoxLayout()
        smiles_label = QLabel("Current SMILES:")
        smiles_label.setStyleSheet("font-weight: bold; color: #ffffff;")
        smiles_layout.addWidget(smiles_label)

        self.smiles_display = QTextEdit()
        self.smiles_display.setMaximumHeight(80)
        self.smiles_display.setPlaceholderText("No structure drawn yet...")
        self.smiles_display.setReadOnly(True)
        if self.initial_smiles:
            self.smiles_display.setPlainText(self.initial_smiles)
        smiles_layout.addWidget(self.smiles_display)
        layout.addLayout(smiles_layout)

        # Buttons
        button_layout = QHBoxLayout()

        self.launch_button = QPushButton("Launch Structure Editor")
        self.launch_button.clicked.connect(self.launch_editor)
        self.launch_button.setMinimumHeight(40)
        button_layout.addWidget(self.launch_button)

        self.ok_button = QPushButton("Use This Structure")
        self.ok_button.clicked.connect(self.accept)
        self.ok_button.setEnabled(False)
        self.ok_button.setMinimumHeight(40)
        button_layout.addWidget(self.ok_button)

        self.cancel_button = QPushButton("Cancel")
        self.cancel_button.clicked.connect(self.reject)
        self.cancel_button.setMinimumHeight(40)
        button_layout.addWidget(self.cancel_button)

        layout.addLayout(button_layout)

        # Apply dark theme
        self.setStyleSheet(
            """
            QDialog {
                background-color: #2b2b2b;
                color: #ffffff;
            }
            QPushButton {
                background-color: #404040;
                border: 1px solid #555555;
                color: #ffffff;
                padding: 8px;
                border-radius: 3px;
            }
            QPushButton:hover {
                background-color: #4CAF50;
            }
            QPushButton:disabled {
                background-color: #333333;
                color: #666666;
            }
            QTextEdit {
                background-color: #404040;
                border: 1px solid #555555;
                color: #ffffff;
                border-radius: 3px;
            }
            """
        )
    
    def launch_editor(self):
        """Launch the web-based structure editor"""
        try:
            self.launch_button.setEnabled(False)
            self.status_label.setText("Starting server...")
            self.progress_bar.setVisible(True)
            self.progress_bar.setRange(0, 0)  # Indeterminate progress
            
            # Check if SMILES-Drawer directory exists
            smiles_drawer_path = os.path.join(os.path.dirname(__file__), 'SMILES-Drawer')
            if not os.path.exists(smiles_drawer_path):
                raise FileNotFoundError(f"SMILES-Drawer directory not found at {smiles_drawer_path}")
            
            # Start the server
            self.server = SmilesDrawingServer(port=8001)
            self.server.start_server()
            
            # Give server time to start
            time.sleep(1)
            
            self.status_label.setText("Launching web browser...")
            
            # Launch browser with the drawing interface
            url = "http://localhost:8001"
            
            # Try to launch in a way that creates a separate window
            try:
                # Try with Chrome app mode first
                import subprocess
                chrome_paths = [
                    "C:/Program Files/Google/Chrome/Application/chrome.exe",
                    "C:/Program Files (x86)/Google/Chrome/Application/chrome.exe",
                    "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
                    "/usr/bin/google-chrome",
                    "/usr/bin/chromium-browser"
                ]
                
                chrome_found = False
                for chrome_path in chrome_paths:
                    if os.path.exists(chrome_path):
                        subprocess.Popen([
                            chrome_path,
                            f"--app={url}",
                            "--window-size=1024,768",
                            "--disable-web-security",
                            "--disable-features=VizDisplayCompositor"
                        ])
                        chrome_found = True
                        break
                
                if not chrome_found:
                    # Fallback to default browser
                    webbrowser.open(url)
                    
            except Exception as e:
                print(f"Error launching Chrome: {e}")
                # Fallback to default browser
                webbrowser.open(url)
            
            self.status_label.setText("Structure editor launched! Draw your structure and click 'Save to Python'")
            self.progress_bar.setVisible(False)
            
            # Start checking for received SMILES
            self.check_timer = QTimer()
            self.check_timer.timeout.connect(self.check_for_smiles)
            self.check_timer.start(1000)  # Check every second
            
            self.launch_button.setText("Editor Launched")
            
        except Exception as e:
            self.progress_bar.setVisible(False)
            self.launch_button.setEnabled(True)
            self.status_label.setText(f"Error: {str(e)}")
            QMessageBox.critical(self, "Error", f"Failed to launch structure editor:\n\n{str(e)}")
    
    def check_for_smiles(self):
        """Check if new SMILES data has been received from the web interface"""
        if self.server:
            smiles = self.server.get_smiles()
            if smiles:
                self.drawn_smiles = smiles
                self.smiles_display.setPlainText(smiles)
                self.ok_button.setEnabled(True)
                self.status_label.setText("Structure received! You can now use this structure or draw another one.")
                
                # Update button text to indicate success
                self.ok_button.setText("Use This Structure ✓")
                self.ok_button.setStyleSheet("""
                    QPushButton {
                        background-color: #4CAF50;
                        border: 1px solid #45a049;
                        color: white;
                        padding: 8px;
                        border-radius: 3px;
                        font-weight: bold;
                    }
                """)
    
    def get_smiles(self):
        """Get the drawn SMILES string"""
        return self.drawn_smiles
    
    def closeEvent(self, event):
        """Clean up when dialog is closed"""
        if self.check_timer:
            self.check_timer.stop()
        if self.server:
            self.server.stop_server()
        event.accept()

class SimpleReactionGUI(QMainWindow):
    """Simple Reaction Condition Predictor GUI"""
    
    def __init__(self):
        super().__init__()
        self.prediction_worker = None
        self.init_ui()
        self.apply_dark_theme()
    
    def init_ui(self):
        """Initialize the user interface"""
        self.setWindowTitle("Reaction Condition Predictor")
        # Removed setGeometry since we're maximizing anyway
        
        # Create central widget and main layout
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QVBoxLayout(central_widget)
        main_layout.setSpacing(8)
        main_layout.setContentsMargins(20, 12, 20, 12)
        
        # Input section
        self.create_input_section(main_layout)
        
        # Reaction scheme section with control buttons
        self.create_reaction_scheme_section(main_layout)
        
        # Results section
        self.create_results_section(main_layout)
        
        # Status bar
        self.statusBar().showMessage("Ready")
        self.statusBar().setStyleSheet("background-color: #424242; color: white;")
        
        # Maximize the window after all UI elements are created
        self.showMaximized()
    
    def create_input_section(self, parent_layout):
        """Create the input section"""
        # SMILES input layout - directly added to parent without groupbox
        smiles_layout = QHBoxLayout()
        smiles_layout.setSpacing(5)
        smiles_label = QLabel("Reaction SMILES:")
        smiles_label.setMinimumWidth(110)
        smiles_label.setMaximumWidth(110)
        self.smiles_input = QLineEdit()
        self.smiles_input.setPlaceholderText("Enter reaction SMILES (e.g., c1ccc(Br)cc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1)")
        self.smiles_input.setMinimumHeight(22)
        self.smiles_input.setMaximumHeight(22)
        
        smiles_layout.addWidget(smiles_label)
        smiles_layout.addWidget(self.smiles_input)
        
        # Add Draw Structure button
        self.draw_structure_button = QPushButton("Draw Structure")
        self.draw_structure_button.setMinimumHeight(22)
        self.draw_structure_button.setMaximumHeight(22)
        self.draw_structure_button.setMaximumWidth(110)
        self.draw_structure_button.setMinimumWidth(110)
        self.draw_structure_button.setToolTip("Open a web-based molecular editor to draw chemical structures")
        self.draw_structure_button.clicked.connect(self.open_structure_drawer)
        self.draw_structure_button.setStyleSheet("""
            QPushButton {
                background-color: #4CAF50;
                border: 1px solid #45a049;
                color: white;
                font-weight: bold;
                border-radius: 3px;
            }
            QPushButton:hover {
                background-color: #45a049;
            }
            QPushButton:pressed {
                background-color: #3d8b40;
            }
        """)
        smiles_layout.addWidget(self.draw_structure_button)
        
        # Add reaction type selection to the right of Draw Structure button
        self.reaction_type_combo = QComboBox()
        self.reaction_type_combo.setMinimumHeight(22)
        self.reaction_type_combo.setMaximumHeight(22)
        self.reaction_type_combo.setMinimumWidth(180)
        self.reaction_type_combo.setMaximumWidth(200)
        
        # Load reaction types
        try:
            reaction_types = get_reaction_types()
            self.reaction_type_combo.addItems(reaction_types)
        except Exception as e:
            # Fallback reaction types if import fails
            fallback_types = [
                "Auto-detect",
                "Suzuki-Miyaura Coupling",
                "Buchwald-Hartwig Amination",
                "Heck Coupling",
                "Sonogashira Coupling",
                "Stille Coupling",
                "Negishi Coupling",
                "Oxidation",
                "Reduction",
                "Esterification",
                "Amidation",
                "Substitution",
                "Elimination",
                "Cycloaddition",
                "General Organic Reaction"
            ]
            self.reaction_type_combo.addItems(fallback_types)
            print(f"Warning: Could not load reaction types, using fallback: {e}")
        
        smiles_layout.addWidget(self.reaction_type_combo)
        
        # Add the layout directly to parent
        parent_layout.addLayout(smiles_layout)
        
        # Connect SMILES input to update reaction scheme
        self.smiles_input.textChanged.connect(self.update_reaction_scheme)
    
    def create_reaction_scheme_section(self, parent_layout):
        """Create the reaction scheme visualization section with control buttons"""
        # Create horizontal layout for scheme and buttons
        scheme_buttons_layout = QHBoxLayout()
        
        # Create reaction image label directly without groupbox
        self.reaction_image_label = QLabel()
        self.reaction_image_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.reaction_image_label.setStyleSheet("""
            QLabel {
                background-color: #404040;
                border: 1px solid #555555;
                border-radius: 5px;
                color: #CCCCCC;
                font-style: italic;
                padding: 5px;
            }
        """)
        self.reaction_image_label.setText("Enter a reaction SMILES to see the reaction scheme")
        # Set size policy to make the label adaptive but with height constraints
        self.reaction_image_label.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.reaction_image_label.setWordWrap(True)
        self.reaction_image_label.setScaledContents(False)  # Don't scale content
        self.reaction_image_label.setFixedHeight(30)  # Fixed height when showing text only
        
        # Add image label directly to horizontal layout
        scheme_buttons_layout.addWidget(self.reaction_image_label, 3)  # Give scheme more space
        
        # Create control buttons section (vertical layout)
        buttons_widget = QWidget()
        buttons_layout = QVBoxLayout(buttons_widget)
        buttons_layout.setContentsMargins(10, 0, 0, 0)
        
        # Predict button
        self.predict_button = QPushButton("Predict Conditions")
        self.predict_button.setMinimumHeight(40)
        self.predict_button.setMaximumWidth(180)
        self.predict_button.clicked.connect(self.run_prediction)
        
        # Clear button
        self.clear_button = QPushButton("Clear")
        self.clear_button.setMinimumHeight(40)
        self.clear_button.setMaximumWidth(180)
        self.clear_button.clicked.connect(self.clear_inputs)
        
        # Browse samples button
        self.browse_samples_button = QPushButton("Browse Samples")
        self.browse_samples_button.setMinimumHeight(40)
        self.browse_samples_button.setMaximumWidth(180)
        self.browse_samples_button.clicked.connect(self.browse_sample_reactions)
        
        buttons_layout.addWidget(self.predict_button)
        buttons_layout.addWidget(self.clear_button)
        buttons_layout.addWidget(self.browse_samples_button)
        buttons_layout.addStretch()  # Push buttons to top
        
        # Add buttons widget to horizontal layout
        scheme_buttons_layout.addWidget(buttons_widget, 1)  # Give buttons less space
        
        parent_layout.addLayout(scheme_buttons_layout)
    
    def create_results_section(self, parent_layout):
        """Create results display section"""
        results_group = QGroupBox("Results")
        results_group.setStyleSheet("""
            QGroupBox {
                font-weight: bold;
                border: 2px solid #555555;
                border-radius: 5px;
                margin-top: 10px;
                padding-top: 10px;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                left: 10px;
                padding: 0 5px 0 5px;
            }
        """)
        results_layout = QVBoxLayout(results_group)
        
        # Create horizontal splitter for text and images
        results_splitter = QSplitter(Qt.Orientation.Horizontal)
        
        # Results text area (left side)
        text_widget = QWidget()
        text_layout = QVBoxLayout(text_widget)
        text_layout.setContentsMargins(0, 0, 5, 0)
        
        text_label = QLabel("Prediction Results:")
        text_label.setStyleSheet("font-weight: bold; color: #ffffff; margin-bottom: 5px;")
        text_layout.addWidget(text_label)
        
        self.results_text = QTextEdit()
        self.results_text.setMinimumHeight(600)  # Further increased from 400 to 600
        # Removed fixed width constraint to allow full expansion
        self.results_text.setPlaceholderText("Prediction results will appear here...")
        self.results_text.setReadOnly(True)
        text_layout.addWidget(self.results_text)
        
        results_splitter.addWidget(text_widget)
        
        # Related reactions images (right side)
        images_widget = QWidget()
        images_layout = QVBoxLayout(images_widget)
        images_layout.setContentsMargins(5, 0, 0, 0)
        
        images_label = QLabel("Top Related Reactions:")
        images_label.setStyleSheet("font-weight: bold; color: #ffffff; margin-bottom: 5px;")
        images_layout.addWidget(images_label)
        
        # Container for reaction images
        self.related_reactions_container = QWidget()
        self.related_reactions_layout = QVBoxLayout(self.related_reactions_container)
        self.related_reactions_layout.setSpacing(10)
        
        # Initially empty - will be populated with reaction images
        placeholder_label = QLabel("Related reactions will appear here after prediction")
        placeholder_label.setStyleSheet("""
            QLabel {
                background-color: #404040;
                border: 1px solid #555555;
                border-radius: 5px;
                color: #CCCCCC;
                font-style: italic;
                padding: 20px;
                text-align: center;
            }
        """)
        placeholder_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        placeholder_label.setMinimumHeight(450)  # Further increased from 300 to 450
        self.related_reactions_layout.addWidget(placeholder_label)
        
        # Scroll area for related reactions
        scroll_area = QScrollArea()
        scroll_area.setWidget(self.related_reactions_container)
        scroll_area.setWidgetResizable(True)
        # Removed width constraints to allow full expansion
        scroll_area.setStyleSheet("""
            QScrollArea {
                background-color: #404040;
                border: 1px solid #555555;
                border-radius: 3px;
            }
        """)
        images_layout.addWidget(scroll_area)
        
        results_splitter.addWidget(images_widget)
        
        # Set splitter proportions: 35% text, 65% images (changed from 60% text, 40% images)
        results_splitter.setSizes([350, 650])
        
        results_layout.addWidget(results_splitter)
        parent_layout.addWidget(results_group)
    
    def apply_dark_theme(self):
        """Apply dark theme to the application"""
        self.setStyleSheet("""
            QMainWindow {
                background-color: #2b2b2b;
                color: #ffffff;
            }
            QWidget {
                background-color: #2b2b2b;
                color: #ffffff;
            }
            QGroupBox {
                color: #ffffff;
                background-color: #3c3c3c;
            }
            QLabel {
                color: #ffffff;
            }
            QLineEdit {
                background-color: #404040;
                border: 1px solid #555555;
                color: #ffffff;
                padding: 2px;
                border-radius: 3px;
            }
            QLineEdit:focus {
                border: 2px solid #4CAF50;
            }
            QComboBox {
                background-color: #404040;
                border: 1px solid #555555;
                color: #ffffff;
                padding: 2px;
                border-radius: 3px;
            }
            QComboBox::drop-down {
                border: none;
            }
            QComboBox::down-arrow {
                image: none;
                border-left: 4px solid transparent;
                border-right: 4px solid transparent;
                border-top: 4px solid #ffffff;
            }
            QTextEdit {
                background-color: #404040;
                border: 1px solid #555555;
                color: #ffffff;
                border-radius: 3px;
            }
        """)
    
    def update_reaction_scheme(self):
        """Update the reaction scheme visualization when SMILES changes"""
        reaction_smiles = self.smiles_input.text().strip()
        
        if not reaction_smiles:
            self.reaction_image_label.clear()
            self.reaction_image_label.setText("Enter a reaction SMILES to see the reaction scheme")
            # Set fixed height for text-only display
            self.reaction_image_label.setFixedHeight(30)
            self.reaction_image_label.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
            return
        
        if ">>" not in reaction_smiles:
            self.reaction_image_label.clear()
            self.reaction_image_label.setText("Invalid reaction SMILES format (missing '>>')")
            # Set fixed height for text-only display
            self.reaction_image_label.setFixedHeight(30)
            self.reaction_image_label.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
            return
        
        try:
            # Create reaction image with compact height
            pixmap = create_reaction_image(reaction_smiles, 500, 70)  # Reduced height for compact display
            
            if pixmap and not pixmap.isNull():
                # Switch to adaptive sizing for image display
                self.reaction_image_label.setMaximumHeight(16777215)  # Remove height constraint
                self.reaction_image_label.setMinimumHeight(0)
                self.reaction_image_label.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum)
                self.reaction_image_label.setPixmap(pixmap)
                # Let the label adapt to the actual image size
                self.reaction_image_label.adjustSize()
            else:
                self.reaction_image_label.clear()
                self.reaction_image_label.setText("Could not generate reaction scheme")
                # Set fixed height for text-only display
                self.reaction_image_label.setFixedHeight(30)
                self.reaction_image_label.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
                
        except Exception as e:
            self.reaction_image_label.clear()
            self.reaction_image_label.setText("Error generating reaction scheme")
            # Set fixed height for text-only display
            self.reaction_image_label.setFixedHeight(30)
            self.reaction_image_label.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
    
    def run_prediction(self):
        """Run prediction with current inputs"""
        # Get inputs
        reaction_smiles = self.smiles_input.text().strip()
        reaction_type = self.reaction_type_combo.currentText()
        
        # Basic validation
        if not reaction_smiles:
            QMessageBox.warning(self, "Input Error", "Please enter a reaction SMILES.")
            return
        
        # Disable predict button during processing
        self.predict_button.setEnabled(False)
        self.predict_button.setText("Processing...")
        
        # Clear previous results
        self.results_text.clear()
        
        # Start prediction worker
        self.prediction_worker = SimplePredictionWorker(reaction_smiles, reaction_type)
        self.prediction_worker.finished.connect(self.on_prediction_finished)
        self.prediction_worker.error.connect(self.on_prediction_error)
        self.prediction_worker.progress.connect(self.on_prediction_progress)
        self.prediction_worker.start()
    
    def display_related_reactions(self, related_reactions):
        """Display images of related reactions from the dataset"""
        # Clear existing content
        for i in reversed(range(self.related_reactions_layout.count())):
            child = self.related_reactions_layout.itemAt(i).widget()
            if child:
                child.setParent(None)
        
        if not related_reactions:
            # Show placeholder if no related reactions
            placeholder_label = QLabel("No related reactions found in dataset")
            placeholder_label.setStyleSheet("""
                QLabel {
                    background-color: #404040;
                    border: 1px solid #555555;
                    border-radius: 5px;
                    color: #CCCCCC;
                    font-style: italic;
                    padding: 20px;
                    text-align: center;
                }
            """)
            placeholder_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
            placeholder_label.setMinimumHeight(100)
            self.related_reactions_layout.addWidget(placeholder_label)
            return
        
        # Display top 2 related reactions
        for i, reaction_data in enumerate(related_reactions[:2], 1):
            # Create container for each reaction
            reaction_container = QWidget()
            reaction_container.setStyleSheet("""
                QWidget {
                    background-color: #3c3c3c;
                    border: 1px solid #555555;
                    border-radius: 5px;
                    margin: 2px;
                }
            """)
            container_layout = QVBoxLayout(reaction_container)
            container_layout.setContentsMargins(8, 8, 8, 8)
            
            # Reaction title
            title_label = QLabel(f"Related Reaction #{i}")
            title_label.setStyleSheet("font-weight: bold; color: #4CAF50; margin-bottom: 5px;")
            container_layout.addWidget(title_label)
            
            # Extract SMILES and create image
            reaction_smiles = reaction_data.get('reaction_smiles', '')
            if not reaction_smiles:
                # Try alternative keys that might contain the SMILES
                reaction_smiles = reaction_data.get('smiles', reaction_data.get('reaction', ''))
            
            if reaction_smiles:
                try:
                    # Create reaction image with wider width to prevent overlap
                    pixmap = create_reaction_image(reaction_smiles, 480, 140)

                    image_label = QLabel()
                    image_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
                    image_label.setStyleSheet(
                        """
                            QLabel {
                                background-color: white;
                                border: 1px solid #777777;
                                border-radius: 3px;
                                padding: 5px;
                            }
                        """
                    )
                    image_label.setScaledContents(True)

                    if pixmap and not pixmap.isNull():
                        image_label.setPixmap(pixmap)
                        container_layout.addWidget(image_label)
                    else:
                        # Fallback to placeholder pixmap
                        placeholder = create_placeholder_image(reaction_smiles, 480, 140)
                        if placeholder and not placeholder.isNull():
                            image_label.setPixmap(placeholder)
                            container_layout.addWidget(image_label)
                        else:
                            # Final textual fallback
                            fallback_label = QLabel("Could not generate reaction image")
                            fallback_label.setStyleSheet("color: #CCCCCC; font-style: italic; text-align: center;")
                            fallback_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
                            container_layout.addWidget(fallback_label)
                except Exception as e:
                    print(f"Error creating image for related reaction {i}: {e}")
                    error_label = QLabel("Error generating image")
                    error_label.setStyleSheet("color: #FF6B6B; font-style: italic; text-align: center;")
                    error_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
                    container_layout.addWidget(error_label)
            
            # Additional info
            info_text = ""
            if 'yield' in reaction_data and reaction_data['yield']:
                try:
                    yield_val = float(reaction_data['yield'])
                    info_text += f"Yield: {yield_val:.1f}%\n"
                except (ValueError, TypeError):
                    pass
            
            if 'catalyst' in reaction_data and reaction_data['catalyst'] != 'N/A':
                info_text += f"Catalyst: {reaction_data['catalyst']}\n"
            
            if 'ligand' in reaction_data and reaction_data['ligand'] != 'N/A':
                info_text += f"Ligand: {reaction_data['ligand']}\n"
            
            if 'solvent' in reaction_data and reaction_data['solvent'] != 'N/A':
                info_text += f"Solvent: {reaction_data['solvent']}\n"
            
            if 'temperature' in reaction_data and reaction_data['temperature']:
                try:
                    temp_val = reaction_data['temperature']
                    if str(temp_val).strip() and str(temp_val).strip() != 'nan':
                        # Prefer numeric formatting with two decimals and space before unit
                        try:
                            temp_float = float(str(temp_val))
                            info_text += f"Temperature: {temp_float:.2f} °C\n"
                        except Exception:
                            # Fall back to string; avoid duplicating unit if already present
                            temp_str = str(temp_val)
                            if '°' in temp_str:
                                info_text += f"Temperature: {temp_str}\n"
                            else:
                                info_text += f"Temperature: {temp_str} °C\n"
                except:
                    pass
            
            if 'time' in reaction_data and reaction_data['time']:
                try:
                    time_val = reaction_data['time']
                    if str(time_val).strip() and str(time_val).strip() != 'nan':
                        # Prefer numeric formatting with two decimals and space before unit
                        try:
                            time_float = float(str(time_val))
                            info_text += f"Time: {time_float:.2f} h\n"
                        except Exception:
                            time_str = str(time_val)
                            if time_str.strip().endswith('h'):
                                info_text += f"Time: {time_str}\n"
                            else:
                                info_text += f"Time: {time_str} h\n"
                except:
                    pass
            
            if 'similarity' in reaction_data:
                info_text += f"Similarity: {reaction_data['similarity']:.1%}\n"
            
            if 'reaction_id' in reaction_data:
                info_text += f"ID: {reaction_data['reaction_id']}"
            
            if info_text:
                info_label = QLabel(info_text.strip())
                info_label.setStyleSheet("""
                    QLabel {
                        color: #CCCCCC;
                        font-size: 9px;
                        background-color: #2b2b2b;
                        border: 1px solid #555555;
                        border-radius: 3px;
                        padding: 4px;
                        margin-top: 5px;
                    }
                """)
                info_label.setWordWrap(True)
                container_layout.addWidget(info_label)
            
            self.related_reactions_layout.addWidget(reaction_container)
        
        # Add stretch to push content to top
        self.related_reactions_layout.addStretch()
    
    def clear_related_reactions(self):
        """Clear the related reactions display"""
        # Clear existing content
        for i in reversed(range(self.related_reactions_layout.count())):
            child = self.related_reactions_layout.itemAt(i).widget()
            if child:
                child.setParent(None)
        
        # Add placeholder
        placeholder_label = QLabel("Related reactions will appear here after prediction")
        placeholder_label.setStyleSheet("""
            QLabel {
                background-color: #404040;
                border: 1px solid #555555;
                border-radius: 5px;
                color: #CCCCCC;
                font-style: italic;
                padding: 20px;
                text-align: center;
            }
        """)
        placeholder_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        placeholder_label.setMinimumHeight(450)  # Further increased from 300 to 450
        self.related_reactions_layout.addWidget(placeholder_label)
    
    def get_related_reactions_from_dataset(self, input_smiles, top_k=2, reaction_type=None):
        """Find related reactions from actual dataset using chemical similarity.

        If reaction_type hints at Ullmann, prefer the Ullmann dataset; otherwise use Buchwald.
        Falls back gracefully when files are not present.
        """
        try:
            import pandas as pd
            import json
            import os
            import hashlib
            from dataset_registry import resolve_dataset_path
            
            # Try to import RDKit for similarity calculations
            try:
                from rdkit import Chem
                from rdkit.Chem import DataStructs, rdMolDescriptors
                rdkit_available = True
            except ImportError:
                rdkit_available = False
                print("RDKit not available - using fallback similarity")
            
            # Resolve dataset using centralized registry
            base_dir = os.path.dirname(__file__)
            dataset_path = resolve_dataset_path(reaction_type, base_dir)
            if not dataset_path:
                print("No reaction dataset found (checked Buchwald and Ullmann).")
                return self.get_mock_related_reactions(input_smiles)

            # Briefly inform which dataset is used
            try:
                self.statusBar().showMessage(f"Related reactions source: {os.path.basename(dataset_path)}")
            except Exception:
                pass
            
            # Read the dataset
            df = pd.read_csv(dataset_path)
            
            if df.empty:
                print("Dataset is empty")
                return self.get_mock_related_reactions(input_smiles)
            
            # Parse input reaction
            if ">>" not in input_smiles:
                print("Invalid input SMILES format")
                return self.get_mock_related_reactions(input_smiles)
            
            input_reactants, input_products = input_smiles.split(">>", 1)
            
            if not rdkit_available:
                # Fallback: simple text similarity
                return self._get_fallback_similar_reactions(df, input_smiles, top_k)
            
            # Calculate chemical similarity using RDKit
            input_reactant_mol = Chem.MolFromSmiles(input_reactants.strip())
            input_product_mol = Chem.MolFromSmiles(input_products.strip())
            
            if input_reactant_mol is None or input_product_mol is None:
                print("Failed to parse input molecules")
                return self._get_fallback_similar_reactions(df, input_smiles, top_k)
            
            # Generate fingerprints for input
            input_reactant_fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(input_reactant_mol, 2, nBits=2048)
            input_product_fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(input_product_mol, 2, nBits=2048)
            
            similarities = []
            
            for idx, row in df.iterrows():
                try:
                    # Get reactant and product SMILES from dataset
                    reactant_smiles = row.get('ReactantSMILES', '')
                    product_smiles = row.get('ProductSMILES', '')
                    
                    if not reactant_smiles or not product_smiles:
                        continue
                    
                    # Parse dataset molecules
                    dataset_reactant_mol = Chem.MolFromSmiles(reactant_smiles)
                    dataset_product_mol = Chem.MolFromSmiles(product_smiles)
                    
                    if dataset_reactant_mol is None or dataset_product_mol is None:
                        continue
                    
                    # Generate fingerprints for dataset reaction
                    dataset_reactant_fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(dataset_reactant_mol, 2, nBits=2048)
                    dataset_product_fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(dataset_product_mol, 2, nBits=2048)
                    
                    # Calculate Tanimoto similarity
                    reactant_similarity = DataStructs.TanimotoSimilarity(input_reactant_fp, dataset_reactant_fp)
                    product_similarity = DataStructs.TanimotoSimilarity(input_product_fp, dataset_product_fp)
                    
                    # Combined similarity (average of reactant and product similarity)
                    combined_similarity = (reactant_similarity + product_similarity) / 2.0
                    
                    # Create reaction entry
                    reaction_entry = {
                        'reaction_smiles': f"{reactant_smiles}>>{product_smiles}",
                        'similarity': combined_similarity,
                        'yield': row.get('Yield_%', 0.0),
                        'catalyst': self._parse_json_field(row.get('CoreDetail', '[]')),
                        'ligand': self._parse_json_field(row.get('Ligand', '[]')),
                        'solvent': self._parse_json_field(row.get('Solvent', '[]')),
                        'temperature': row.get('Temperature_C', ''),
                        'time': row.get('Time_h', ''),
                        'reference': row.get('Reference', ''),
                        'reaction_id': row.get('ReactionID', f'RXN_{idx}')
                    }
                    
                    similarities.append((combined_similarity, reaction_entry))
                    
                except Exception as e:
                    print(f"Error processing reaction {idx}: {e}")
                    continue
            
            if not similarities:
                print("No valid similarities calculated")
                return self._get_fallback_similar_reactions(df, input_smiles, top_k)
            
            # Sort by similarity (highest first) and return top_k
            similarities.sort(key=lambda x: x[0], reverse=True)
            top_reactions = [reaction for _, reaction in similarities[:top_k]]
            
            print(f"Found {len(top_reactions)} similar reactions from dataset")
            return top_reactions
            
        except Exception as e:
            print(f"Error in dataset similarity calculation: {e}")
            # Fallback to mock data if anything goes wrong
            return self.get_mock_related_reactions(input_smiles)
    
    def _parse_json_field(self, json_str):
        """Parse JSON field from dataset, handling various formats"""
        import pandas as pd
        import json
        
        try:
            if pd.isna(json_str) or json_str == '':
                return 'N/A'
            
            # Handle JSON array format
            if json_str.startswith('[') and json_str.endswith(']'):
                parsed = json.loads(json_str)
                if isinstance(parsed, list) and len(parsed) > 0:
                    # Heuristic: Ullmann dataset sometimes tokenizes a single long name
                    # into many short fragments. If many short tokens, join with spaces.
                    try:
                        items = [str(x).strip() for x in parsed if str(x).strip()]
                        if not items:
                            return 'N/A'
                        avg_len = sum(len(it) for it in items) / max(1, len(items))
                        if len(items) > 3 and avg_len < 8:
                            return " ".join(items)
                        # Otherwise, join multiple entries with comma for readability
                        return ", ".join(items)
                    except Exception:
                        # Fallback to first item
                        return str(parsed[0])
                return 'N/A'
            
            # Handle plain string
            return str(json_str)
            
        except (json.JSONDecodeError, ValueError):
            return str(json_str) if json_str else 'N/A'
    
    def _get_fallback_similar_reactions(self, df, input_smiles, top_k):
        """Fallback similarity when RDKit is not available"""
        import json
        
        try:
            # Simple approach: randomly select reactions but make it deterministic
            import hashlib
            hash_obj = hashlib.md5(input_smiles.encode())
            hash_int = int(hash_obj.hexdigest()[:8], 16)
            
            # Select random but deterministic subset
            n_reactions = min(len(df), 20)  # Work with subset for efficiency
            start_idx = hash_int % max(1, len(df) - n_reactions)
            subset_df = df.iloc[start_idx:start_idx + n_reactions]
            
            selected_reactions = []
            
            for idx, row in subset_df.head(top_k).iterrows():
                try:
                    reactant_smiles = row.get('ReactantSMILES', '')
                    product_smiles = row.get('ProductSMILES', '')
                    
                    if not reactant_smiles or not product_smiles:
                        continue
                    
                    reaction_entry = {
                        'reaction_smiles': f"{reactant_smiles}>>{product_smiles}",
                        'similarity': 0.5 + (hash_int % 30) / 100.0,  # Fake similarity 0.5-0.8
                        'yield': row.get('Yield_%', 0.0),
                        'catalyst': self._parse_json_field(row.get('CoreDetail', '[]')),
                        'ligand': self._parse_json_field(row.get('Ligand', '[]')),
                        'solvent': self._parse_json_field(row.get('Solvent', '[]')),
                        'temperature': row.get('Temperature_C', ''),
                        'time': row.get('Time_h', ''),
                        'reference': row.get('Reference', ''),
                        'reaction_id': row.get('ReactionID', f'RXN_{idx}')
                    }
                    
                    selected_reactions.append(reaction_entry)
                    
                except Exception as e:
                    print(f"Error processing fallback reaction {idx}: {e}")
                    continue
            
            return selected_reactions
            
        except Exception as e:
            print(f"Error in fallback similarity: {e}")
            return self.get_mock_related_reactions(input_smiles)

    def get_mock_related_reactions(self, input_smiles):
        """Generate mock related reactions for demonstration purposes"""
        # This is a temporary function to demonstrate the feature
        # In a real implementation, this would come from the recommendation engine
        
        # Diverse set of Buchwald-Hartwig mock reactions
        all_mock_reactions = [
            {
                'reaction_smiles': 'Brc1ccc(C)cc1.CN(C)C>>CN(C)c1ccc(C)cc1',
                'yield': 92.5,
                'catalyst': 'Pd2(dba)3',
                'ligand': 'XPhos',
                'similarity': 0.87
            },
            {
                'reaction_smiles': 'Brc1ccc(OC)cc1.NCc1ccccc1>>COc1ccc(NCc2ccccc2)cc1',
                'yield': 88.3,
                'catalyst': 'Pd(OAc)2',
                'ligand': 'SPhos',
                'similarity': 0.82
            },
            {
                'reaction_smiles': 'Ic1cccc(C(F)(F)F)c1.Nc1ccccc1>>FC(F)(F)c1cccc(Nc2ccccc2)c1',
                'yield': 76.8,
                'catalyst': 'Pd(OAc)2',
                'ligand': 'XPhos',
                'similarity': 0.79
            },
            {
                'reaction_smiles': 'Clc1ccc(cc1)C(=O)OC.N1CCOCC1>>O=C(c1ccc(N2CCOCC2)cc1)OC',
                'yield': 84.2,
                'catalyst': 'Pd2(dba)3',
                'ligand': 'RuPhos',
                'similarity': 0.73
            },
            {
                'reaction_smiles': 'Brc1cncc(C)c1.CNc1ccccc1>>CN(c1ccccc1)c1cncc(C)c1',
                'yield': 71.5,
                'catalyst': 'Pd(OAc)2',
                'ligand': 'DavePhos',
                'similarity': 0.68
            },
            {
                'reaction_smiles': 'Brc1ccc(CN)cc1.c1ccc(N)cc1>>c1ccc(Nc2ccc(CN)cc2)cc1',
                'yield': 89.7,
                'catalyst': 'Pd2(dba)3',
                'ligand': 'BINAP',
                'similarity': 0.91
            },
            {
                'reaction_smiles': 'Clc1nc(Cl)nc(N)n1.NCCc1ccccc1>>NCCc1ccccc1.Nc1nc(NCCc2ccccc2)nc(N)n1',
                'yield': 65.4,
                'catalyst': 'Pd(OAc)2',
                'ligand': 'XantPhos',
                'similarity': 0.64
            },
            {
                'reaction_smiles': 'Brc1ccc(F)cc1.N1CCNCC1>>Fc1ccc(N2CCNCC2)cc1',
                'yield': 93.1,
                'catalyst': 'Pd2(dba)3',
                'ligand': 'XPhos',
                'similarity': 0.85
            }
        ]
        
        # Simple hash-based selection to ensure variety based on input
        import hashlib
        hash_obj = hashlib.md5(input_smiles.encode())
        hash_int = int(hash_obj.hexdigest()[:8], 16)
        
        # Select 2 different reactions based on input hash
        start_idx = hash_int % (len(all_mock_reactions) - 1)
        selected_reactions = []
        
        # First reaction
        selected_reactions.append(all_mock_reactions[start_idx])
        
        # Second reaction (ensure it's different)
        second_idx = (start_idx + 2) % len(all_mock_reactions)
        selected_reactions.append(all_mock_reactions[second_idx])
        
        return selected_reactions
    
    def on_prediction_finished(self, result):
        """Handle prediction completion with enhanced results"""
        self.predict_button.setEnabled(True)
        self.predict_button.setText("Predict Conditions")

        # Choose formatter based on analysis type
        analysis_type = result.get('analysis_type')
        if analysis_type in ('enhanced', 'comprehensive'):
            results_text = self._format_enhanced_recommendations(result)
        elif analysis_type == 'buchwald_hartwig':
            results_text = self._format_buchwald_recommendations(result)
        elif analysis_type == 'general':
            results_text = self._format_general_recommendations(result)
        else:
            results_text = self._format_basic_results(result)

        self.results_text.setPlainText(results_text)

        # Related reactions (compute before export so we can include in files)
        related_reactions = result.get('recommendations', {}).get('related_reactions', [])
        if not related_reactions:
            reaction_smiles = result.get('reaction_smiles', '')
            if reaction_smiles:
                rt_hint = (
                    result.get('recommendations', {}).get('reaction_type')
                    or result.get('reaction_type')
                    or self.reaction_type_combo.currentText()
                )
                related_reactions = self.get_related_reactions_from_dataset(reaction_smiles, reaction_type=rt_hint)
                if not related_reactions:
                    detected_type = result.get('recommendations', {}).get('reaction_type', '')
                    if (analysis_type == 'buchwald_hartwig' or 
                        'buchwald' in detected_type.lower() or 
                        'cross' in detected_type.lower() or
                        'coupling' in detected_type.lower()):
                        related_reactions = self.get_mock_related_reactions(reaction_smiles)

        # Export a well-structured JSON (no separate HTML export)
        try:
            export_dir = os.path.join(os.getcwd(), 'exports')
            os.makedirs(export_dir, exist_ok=True)
            ts = time.strftime('%Y%m%d-%H%M%S')

            # Use shared export builder
            try:
                from prediction_export import build_export_payload  # type: ignore
                export = build_export_payload(result, related_reactions)
            except Exception:
                export = self._build_export_payload(result, related_reactions=related_reactions)

            out_path = os.path.join(export_dir, f'prediction_{ts}.json')
            with open(out_path, 'w', encoding='utf-8') as f:
                json.dump(export, f, ensure_ascii=False, indent=2)
            self.statusBar().showMessage(f"Prediction completed • Exported JSON: {os.path.basename(out_path)}")
        except Exception as e:
            self.statusBar().showMessage(f"Prediction completed • JSON export failed: {e}")

        # Related reactions display
        self.display_related_reactions(related_reactions)

    def _build_export_payload(self, result: dict, related_reactions=None) -> dict:
        """Convert internal result into a clean, stable JSON payload.

                Schema:
                    {
                        meta: { generated_at, analysis_type, status },
                        input: { reaction_smiles, selected_reaction_type },
                        detection: { reaction_type },
                        dataset: { ligands_available, solvents_available, reaction_types_supported },
                        top_conditions: [ for top 3, a single 'chemicals' list with all required components (starting materials, metal precursor, ligand, base, solvent) and 'conditions' with time/temperature ],
                        analytics: { source, top: { ligands[], solvents[], bases[] }, cooccurrence: { best_ligand_solvent, best_base_solvent }, numeric_stats: { temperature_c, time_h, yield_pct }, typical_catalyst_loading },
                        related_reactions: [ { reaction_smiles, yield, catalyst, ligand, solvent, temperature, time, similarity, reaction_id, reference } ]
                    }
        """
        recs = result.get('recommendations', {}) or {}

        combined = []
        for c in recs.get('combined_conditions', []) or []:
            combined.append({
                'ligand': c.get('ligand'),
                'ligand_compatibility': c.get('ligand_compatibility'),
                'solvent': c.get('solvent'),
                'solvent_abbreviation': c.get('solvent_abbreviation'),
                'solvent_compatibility': c.get('solvent_compatibility'),
                'combined_score': c.get('combined_score'),
                'recommendation_confidence': c.get('recommendation_confidence'),
                'typical_conditions': c.get('typical_conditions', {}),
                'synergy_bonus': c.get('synergy_bonus', 0),
                'suggested_base': c.get('suggested_base')
            })

        ligands = []
        for l in recs.get('ligand_recommendations', []) or []:
            ligands.append({
                'ligand': l.get('ligand'),
                'compatibility_score': l.get('compatibility_score'),
                'applications': l.get('applications'),
                'reaction_suitability': l.get('reaction_suitability'),
            })

        solvents = []
        for s in recs.get('solvent_recommendations', []) or []:
            solvents.append({
                'solvent': s.get('solvent'),
                'abbreviation': s.get('abbreviation'),
                'compatibility_score': s.get('compatibility_score'),
                'applications': s.get('applications'),
                'reaction_suitability': s.get('reaction_suitability'),
            })

        alternatives = {}
        alt = recs.get('property_based_alternatives', {}) or {}
        for key in ['budget_friendly_ligands', 'low_boiling_solvents', 'green_solvents']:
            if key in alt:
                alternatives[key] = alt[key]

        related = related_reactions or recs.get('related_reactions', []) or []

    # Build top 3 detailed conditions as a flat chemicals list
        def _split_smiles(sm: str):
            try:
                lhs, rhs = (sm or '').split('>>', 1)
                reactants = [t for t in lhs.split('.') if t]
                product = rhs
                return reactants, product
            except Exception:
                return [], sm

        def _solvent_cas(name: str | None):
            try:
                from reagents.solvent import create_solvent_dataframe  # type: ignore
                df = create_solvent_dataframe()
                if name and 'Solvent' in df.columns and 'CAS Number' in df.columns:
                    row = df[df['Solvent'] == name]
                    if not row.empty:
                        return row.iloc[0].get('CAS Number')
            except Exception:
                pass
            return None

        rxn_smiles = result.get('reaction_smiles', '')
        reactants, _ = _split_smiles(rxn_smiles)
        detected_type = recs.get('reaction_type') or result.get('reaction_type') or ''

        def _default_metal_precursor(rt: str):
            if isinstance(rt, str) and rt.lower() == 'ullmann':
                return {'name': 'CuI', 'cas': '7681-65-4', 'smiles': None, 'equivalents': None}
            else:
                return {'name': 'Pd(OAc)2', 'cas': '3375-31-3', 'smiles': None, 'equivalents': None}

        # Lightweight CAS lookup from data/cas_dictionary.csv for ligands/bases
        def _load_cas_map():
            cas_by_token = {}
            cas_by_name = {}
            try:
                data_dir = os.path.join(os.path.dirname(__file__), 'data')
                path = os.path.join(data_dir, 'cas_dictionary.csv')
                with open(path, 'r', encoding='utf-8') as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        cas = (row.get('CAS') or '').strip()
                        name = (row.get('Name') or '').strip()
                        token = (row.get('Token') or '').strip()
                        if cas:
                            if token:
                                cas_by_token[token.lower()] = cas
                            if name:
                                cas_by_name[name.lower()] = cas
            except Exception:
                pass
            # Common synonyms/aliases
            aliases = {
                'k2co3': 'K2CO3',
                'cs2co3': 'Cs2CO3',
                'k3po4': 'K3PO4',
                'kotbu': 'KOtBu',
                'potassium tert-butoxide': 'KOtBu',
                'naotbu': 'NaOtBu',
                'dipea': 'DIPEA',
                'dbu': 'DBU',
                'pyridine': 'Pyridine',
                'xphos': 'XPhos',
                'sphos': 'SPhos',
                'ruphos': 'RuPhos',
                'brettphos': 'BrettPhos',
                'tbuxphos': 'tBuXPhos',
                'johnphos': 'JohnPhos',
                'xantphos': 'XantPhos',
                'dppe': 'DPPE',
                'dppf': 'DPPF',
                'binap': 'BINAP',
            }
            return cas_by_token, cas_by_name, aliases

        cas_by_token, cas_by_name, alias_map = _load_cas_map()

        def _lookup_cas(name: str | None):
            if not name:
                return None
            key = name.strip().lower()
            # direct token hit
            if key in cas_by_token:
                return cas_by_token[key]
            # direct name hit
            if key in cas_by_name:
                return cas_by_name[key]
            # alias to token/name
            if key in alias_map:
                alias = alias_map[key]
                ak = alias.lower()
                if ak in cas_by_token:
                    return cas_by_token[ak]
                if ak in cas_by_name:
                    return cas_by_name[ak]
            return None

        top_conditions = []
        for c in (recs.get('combined_conditions', []) or [])[:3]:
            conditions = c.get('typical_conditions', {}) or {}
            base_name = c.get('suggested_base') or conditions.get('base')
            chemicals = []
            # Starting materials from reactants
            for smi in reactants:
                chemicals.append({'name': None, 'cas': None, 'smiles': smi, 'equivalents': None, 'role': 'starting_material'})
            # Metal precursor (based on detected type)
            mp = _default_metal_precursor(detected_type)
            chemicals.append({**mp, 'role': 'metal_precursor'})
            # Ligand with CAS lookup
            lig_name = c.get('ligand')
            chemicals.append({'name': lig_name, 'cas': _lookup_cas(lig_name), 'smiles': None, 'equivalents': None, 'role': 'ligand'})
            # Base
            if base_name:
                # Default 2.0 eq for bases if not specified
                chemicals.append({'name': base_name, 'cas': _lookup_cas(base_name), 'smiles': None, 'equivalents': 2.0, 'role': 'base'})
            # Solvent
            chemicals.append({'name': c.get('solvent'), 'abbreviation': c.get('solvent_abbreviation'), 'cas': _solvent_cas(c.get('solvent')), 'smiles': None, 'equivalents': None, 'role': 'solvent'})

            top_conditions.append({
                'reaction': {'smiles': rxn_smiles},
                'chemicals': chemicals,
                'conditions': {
                    'temperature': conditions.get('temperature'),
                    'time': conditions.get('time'),
                    'atmosphere': conditions.get('atmosphere'),
                }
            })

        # Optional analytics snippet (Ullmann) if available
        def _load_analytics_snippet(rt: str | None):
            try:
                if not isinstance(rt, str) or rt.strip().lower() != 'ullmann':
                    return None
                base = os.path.join(os.path.dirname(__file__), 'data', 'analytics', 'Ullmann')
                latest = os.path.join(base, 'latest.json')
                if not os.path.exists(latest):
                    return None
                with open(latest, 'r', encoding='utf-8') as f:
                    summ = json.load(f)
                top = (summ.get('top') or {})
                co = (summ.get('cooccurrence') or {})
                nums = (summ.get('numeric_stats') or {})
                def _simple(items, n=3):
                    out = []
                    for it in (items or [])[:n]:
                        out.append({
                            'name': it.get('name'),
                            'pct': it.get('pct'),
                            'count': it.get('count')
                        })
                    return out
                def _best(pair_list):
                    if not pair_list:
                        return None
                    first = pair_list[0]
                    return {
                        'a': first.get('a'),
                        'b': first.get('b'),
                        'pct': first.get('pct'),
                        'count': first.get('count')
                    }
                snippet = {
                    'source': 'Ullmann',
                    'top': {
                        'ligands': _simple(top.get('ligands')),
                        'solvents': _simple(top.get('solvents')),
                        'bases': _simple(top.get('bases')),
                    },
                    'cooccurrence': {
                        'best_ligand_solvent': _best(co.get('ligand_solvent')),
                        'best_base_solvent': _best(co.get('base_solvent')),
                    },
                    'numeric_stats': {
                        'temperature_c': nums.get('temperature_c'),
                        'time_h': nums.get('time_h'),
                        'yield_pct': nums.get('yield_pct'),
                    }
                }
                try:
                    for c in (recs.get('combined_conditions') or []):
                        tc = c.get('typical_conditions') or {}
                        if tc.get('catalyst_loading'):
                            snippet['typical_catalyst_loading'] = tc.get('catalyst_loading')
                            break
                except Exception:
                    pass
                return snippet
            except Exception:
                return None

        analytics_snippet = _load_analytics_snippet(detected_type)

        # Embed analytics inside dataset for simplified schema
        dataset_block = recs.get('dataset_info', {}) or {}
        if not isinstance(dataset_block, dict):
            dataset_block = {}
        if analytics_snippet:
            dataset_block = dict(dataset_block)
            dataset_block['analytics'] = analytics_snippet

        payload = {
            'meta': {
                'generated_at': time.strftime('%Y-%m-%dT%H:%M:%S'),
                'analysis_type': result.get('analysis_type', 'unknown'),
                'status': result.get('status'),
            },
            'input': {
                'reaction_smiles': result.get('reaction_smiles'),
                'selected_reaction_type': result.get('reaction_type') or result.get('selected_reaction_type'),
            },
            'detection': {
                'reaction_type': recs.get('reaction_type'),
            },
            'dataset': dataset_block,
            'top_conditions': top_conditions,
            'related_reactions': related,
        }
        return payload
    
    def _format_buchwald_recommendations(self, result):
        """Format Buchwald-Hartwig specific recommendations"""
        
        text = f"""Buchwald-Hartwig Amination Condition Recommendations

Reaction: {result['reaction_smiles']}
Status: {result['status']}

"""
        
        if 'error' in result:
            text += f"Error: {result['error']}\n"
            return text
        
        recommendations_data = result.get('recommendations', {})
        
        if 'error' in recommendations_data:
            text += f"Recommendation Error: {recommendations_data['error']}\n"
            return text
        
        dataset_size = recommendations_data.get('dataset_size', 'Unknown')
        text += f"Dataset Size: {dataset_size} reactions analyzed\n\n"
        
        recs = recommendations_data.get('recommendations', [])
        
        if not recs:
            text += """No specific recommendations found.

This could be due to:
• Limited structural similarity to known reactions
• Insufficient data in the current dataset
• Novel substrate combination

Suggestions:
• Check literature for similar substrates
• Start with common Buchwald-Hartwig conditions (Pd2(dba)3/XPhos)
• Consider substrate electronic properties when choosing ligands
"""
            return text
        
        text += f"Found {len(recs)} catalyst family recommendations:\n"
        text += "=" * 70 + "\n\n"
        
        for i, rec in enumerate(recs, 1):
            text += f"""RANK #{i}: {rec['family_name']}

Recommended Catalyst System:
  • Metal Precursor: {rec['recommended_catalyst']}
  • Ligand: {rec['recommended_ligand']}

Performance Statistics:
  • Average Yield: {rec['performance']['avg_yield']:.1f}%
  • Best Yield: {rec['performance']['max_yield']:.1f}%
  • Success Rate (≥80%): {rec['performance']['success_rate']:.1%}
  • Literature Examples: {rec['performance']['total_reactions']} reactions
  • Confidence Score: {rec['confidence_score']:.1%}

Alternative Catalysts in this Family:"""
            
            if rec['alternatives']:
                for alt in rec['alternatives']:
                    text += f"""
  → {alt['catalyst']} / {alt['ligand']} (Avg: {alt['avg_yield']:.1f}%, Success: {alt['success_rate']:.1%})"""
            else:
                text += "\n  → No alternatives found in dataset"
            
            text += "\n\n" + "-" * 70 + "\n\n"
        
        text += """INTERPRETATION GUIDE:

• Higher ranked families show better overall performance
• Confidence scores reflect the amount of supporting literature data
• Alternative catalysts within each family are likely to give similar results
• Consider substrate sterics when choosing between families:
  - Bulky substrates: Use bulky phosphine ligands (XPhos, SPhos)
  - Electron-poor substrates: Use electron-rich ligands
  - Standard substrates: Moderate phosphines often work well

• Metal precursor selection tips:
  - Pd2(dba)3: Pre-reduced, good for sensitive substrates
  - Pd(OAc)2: Most common, requires in-situ reduction
  - PdCl2: Alternative Pd(II) source, similar to Pd(OAc)2
"""
        
        return text
    
    def _format_general_recommendations(self, result):
        """Format general recommendations"""
        recommendations = result.get('recommendations', {})
        
        # Build a minimal display type with metal if we can infer from selected type
        def _display_general_type():
            selected = result.get('selected_reaction_type') or ''
            low = selected.lower()
            if any(k in low for k in ["suzuki", "stille", "sonogashira", "heck", "negishi", "buchwald", "hartwig"]):
                return f"{selected} (Pd)" if selected else "Pd"
            if any(k in low for k in ["chan-lam", "ullmann"]):
                return f"{selected} (Cu)" if selected else "Cu"
            if "kumada" in low:
                return f"{selected} (Ni)" if selected else "Ni"
            return selected or "Unknown"

        text = f"""General Reaction Analysis

Reaction: {result['reaction_smiles']}
Status: {result['status']}

{recommendations.get('message', 'No specific recommendations available.')}

"""
        
        suggestions = recommendations.get('suggestions', [])
        if suggestions:
            text += "Suggestions:\n"
            for i, suggestion in enumerate(suggestions, 1):
                text += f"  {i}. {suggestion}\n"
        
        text += f"""

Available Recommenders: {', '.join(result.get('available_recommenders', []))}

Note: For specific reaction types like Buchwald-Hartwig amination, 
the system can provide detailed catalyst recommendations when the 
appropriate dataset is available.
"""
        
        return text
    
    def _format_enhanced_recommendations(self, result):
        """Format enhanced ligand and solvent recommendations"""
        recommendations = result.get('recommendations', {})
        
        # Build a clearer detected reaction type with metal annotation when possible
        def _display_reaction_type():
            rt = recommendations.get('reaction_type', 'Unknown') or 'Unknown'
            origin = recommendations.get('detected_from') or result.get('selected_reaction_type') or ''
            def _metal_for(name: str | None):
                if not name:
                    return None
                low = name.lower()
                if 'ullmann' in low:
                    return 'Cu'
                if 'chan-lam' in low:
                    return 'Cu'
                if 'kumada' in low:
                    return 'Ni'
                if 'buchwald' in low or 'hartwig' in low:
                    return 'Pd'
                if any(k in low for k in ['suzuki', 'stille', 'sonogashira', 'heck', 'negishi']):
                    return 'Pd'
                if 'c-o' in low and 'ullmann' in low:
                    return 'Cu'
                if 'c-s' in low:
                    return 'Pd'
                return None
            metal = _metal_for(rt) or _metal_for(origin) or ('Pd' if (rt or '').lower() == 'cross-coupling' else None)
            if metal:
                name_part = rt if (rt and rt.lower() not in ['cross-coupling', 'general organic reaction']) else ''
                return f"{name_part} ({metal})" if name_part else f"{metal}"
            return rt

        text = f"""🧪 ENHANCED REACTION CONDITION RECOMMENDATIONS

Reaction: {result['reaction_smiles']}
Detected Type: {_display_reaction_type()}
Status: {result['status']}

"""
        
        if 'error' in recommendations:
            text += f"❌ Error: {recommendations['error']}\n\n"
            return text
        
        # Dataset information
        dataset_info = recommendations.get('dataset_info', {})
        if dataset_info:
            text += f"""📊 Database Coverage:
• Available Ligands: {dataset_info.get('ligands_available', 'N/A')}
• Available Solvents: {dataset_info.get('solvents_available', 'N/A')}
• Supported Reactions: {', '.join(dataset_info.get('reaction_types_supported', []))}

"""

        # Optional analytics snippet (Ullmann) if available
        def _load_analytics_snippet(rt: str | None):
            try:
                if not isinstance(rt, str) or rt.strip().lower() != 'ullmann':
                    return None
                base = os.path.join(os.path.dirname(__file__), 'data', 'analytics', 'Ullmann')
                latest = os.path.join(base, 'latest.json')
                if not os.path.exists(latest):
                    return None
                with open(latest, 'r', encoding='utf-8') as f:
                    summ = json.load(f)
                top = (summ.get('top') or {})
                co = (summ.get('cooccurrence') or {})
                return {
                    'ligands': top.get('ligands') or [],
                    'solvents': top.get('solvents') or [],
                    'bases': top.get('bases') or [],
                    'ligand_solvent': co.get('ligand_solvent') or [],
                    'base_solvent': co.get('base_solvent') or [],
                }
            except Exception:
                return None

        if dataset_info.get('analytics_loaded'):
            snip = _load_analytics_snippet(recommendations.get('reaction_type'))
            if snip:
                def _fmt_top(items):
                    out = []
                    for it in (items or [])[:3]:
                        name = it.get('name')
                        pct = it.get('pct')
                        if pct is not None:
                            pct = f"{pct*100:.1f}%"
                        out.append(f"{name} ({pct})")
                    return ", ".join([s for s in out if s])
                def _fmt_pair(items):
                    if not items:
                        return None
                    first = items[0]
                    pct = first.get('pct')
                    if pct is not None:
                        pct = f"{pct*100:.1f}%"
                    return f"{first.get('a')} + {first.get('b')} ({pct})"

                text += "📈 Analytics (Ullmann dataset):\n"
                tl = _fmt_top(snip.get('ligands'))
                ts = _fmt_top(snip.get('solvents'))
                tb = _fmt_top(snip.get('bases'))
                if tl:
                    text += f"• Top Ligands: {tl}\n"
                if ts:
                    text += f"• Top Solvents: {ts}\n"
                if tb:
                    text += f"• Top Bases: {tb}\n"
                pls = _fmt_pair(snip.get('ligand_solvent'))
                pbs = _fmt_pair(snip.get('base_solvent'))
                if pls:
                    text += f"• Best Ligand–Solvent Pair: {pls}\n"
                if pbs:
                    text += f"• Best Base–Solvent Pair: {pbs}\n"
                text += "\n"

        # Combined condition recommendations (top section)
        combined_conditions = recommendations.get('combined_conditions', [])
        if combined_conditions:
            text += "🎯 TOP RECOMMENDED REACTION CONDITIONS:\n"
            text += "=" * 75 + "\n\n"
            
            for i, combo in enumerate(combined_conditions[:3], 1):
                confidence_emoji = "🟢" if combo['recommendation_confidence'] == 'High' else "🟡" if combo['recommendation_confidence'] == 'Medium' else "🔴"
                
                text += f"""RECOMMENDATION #{i} {confidence_emoji} {combo['recommendation_confidence']} Confidence

🧬 Catalyst System:
  • Ligand: {combo['ligand']} (Score: {combo['ligand_compatibility']})
  • Solvent: {combo['solvent']} ({combo['solvent_abbreviation']}) (Score: {combo['solvent_compatibility']})
  • Combined Score: {combo['combined_score']}"""
                
                if combo['synergy_bonus'] > 0:
                    text += f" (includes +{combo['synergy_bonus']} synergy bonus)"
                
                # Add typical conditions
                conditions = combo.get('typical_conditions', {})
                if conditions:
                    text += f"""

⚙️ Typical Conditions:
  • Temperature: {conditions.get('temperature', 'N/A')}
  • Time: {conditions.get('time', 'N/A')}
  • Atmosphere: {conditions.get('atmosphere', 'N/A')}"""
                    
                    if 'catalyst_loading' in conditions:
                        text += f"\n  • Catalyst Loading: {conditions['catalyst_loading']}"
                    if 'base' in conditions:
                        text += f"\n  • Base: {conditions['base']}"
                    if 'additives' in conditions:
                        text += f"\n  • Notes: {conditions['additives']}"
                
                text += "\n\n" + "-" * 75 + "\n\n"
        
        # Individual ligand recommendations
        ligand_recs = recommendations.get('ligand_recommendations', [])
        if ligand_recs:
            text += "🧬 TOP LIGAND OPTIONS:\n"
            text += "-" * 40 + "\n"
            
            for i, lig in enumerate(ligand_recs[:5], 1):
                text += f"""  {i}. {lig['ligand']} (Score: {lig['compatibility_score']})
     • Applications: {lig['applications']}
     • Suitability: {lig['reaction_suitability']}
"""
        
        # Individual solvent recommendations  
        solvent_recs = recommendations.get('solvent_recommendations', [])
        if solvent_recs:
            text += "\n🧪 TOP SOLVENT OPTIONS:\n"
            text += "-" * 40 + "\n"
            
            for i, sol in enumerate(solvent_recs[:5], 1):
                text += f"""  {i}. {sol['solvent']} ({sol['abbreviation']}) (Score: {sol['compatibility_score']})
     • Applications: {sol['applications']}
     • Suitability: {sol['reaction_suitability']}
"""
        
        # Property-based alternatives
        alternatives = recommendations.get('property_based_alternatives', {})
        if alternatives and len(alternatives) > 1:  # More than just error key
            text += "\n🎛️ SPECIALIZED OPTIONS:\n"
            text += "-" * 40 + "\n"
            
            if 'budget_friendly_ligands' in alternatives:
                text += "\n💰 Budget-Friendly Ligands:\n"
                for lig in alternatives['budget_friendly_ligands']:
                    text += f"  • {lig['ligand']} (Score: {lig['compatibility_score']})\n"
            
            if 'low_boiling_solvents' in alternatives:
                text += "\n🌡️ Low-Boiling Solvents (Easy Removal):\n"
                for sol in alternatives['low_boiling_solvents']:
                    text += f"  • {sol['solvent']} ({sol['abbreviation']}) (Score: {sol['compatibility_score']})\n"
            
            if 'green_solvents' in alternatives:
                text += "\n🌱 Green/Sustainable Solvents:\n"
                for sol in alternatives['green_solvents']:
                    text += f"  • {sol['solvent']} ({sol['abbreviation']}) (Score: {sol['compatibility_score']})\n"
        
        # Reaction-specific guidance
        guidance = recommendations.get('reaction_specific_notes', '')
        if guidance:
            text += f"\n{guidance}\n"
        
        # Add footer
        text += f"""
📚 REFERENCES & NOTES:
• Compatibility scores: 0.9-1.0 (Excellent), 0.7-0.8 (Good), 0.5-0.6 (Fair)
• Combined scores include synergy bonuses for proven ligand-solvent pairs
• Start with highest-scoring recommendations and optimize from there
• Consider substrate sterics and electronics when making final selection

🔬 For best results: Use inert atmosphere, dry solvents, and degassed conditions.
"""
        
        return text
    
    def _format_basic_results(self, result):
        """Format basic results for fallback cases"""
        text = f"""Prediction Completed

Input Information:
• Reaction SMILES: {result['reaction_smiles']}
• Selected Reaction Type: {result.get('selected_reaction_type', 'Not specified')}

Status: {result['status']}
"""
        
        if 'message' in result:
            text += f"Message: {result['message']}\n"
        
        text += """
Note: The recommendation engine is being enhanced to provide 
specific condition recommendations for different reaction types.
"""
        
        return text
    
    def on_prediction_error(self, error_message):
        """Handle prediction error"""
        self.predict_button.setEnabled(True)
        self.predict_button.setText("Predict Conditions")
        
        self.results_text.setPlainText(f"Error occurred during prediction:\n\n{error_message}")
        self.statusBar().showMessage("Prediction failed")
        
        QMessageBox.critical(self, "Prediction Error", error_message)
    
    def on_prediction_progress(self, value):
        """Handle prediction progress updates"""
        self.statusBar().showMessage(f"Processing... {value}%")
    
    def clear_inputs(self):
        """Clear all input fields"""
        self.smiles_input.clear()
        self.reaction_type_combo.setCurrentIndex(0)
        self.results_text.clear()
        
        # Clear reaction scheme and reset to compact text display
        self.reaction_image_label.clear()
        self.reaction_image_label.setText("Enter a reaction SMILES to see the reaction scheme")
        self.reaction_image_label.setFixedHeight(30)
        self.reaction_image_label.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        
        # Clear related reactions display
        self.clear_related_reactions()
        
        self.statusBar().showMessage("Inputs cleared")
    
    def open_structure_drawer(self):
        """Open the structure drawing interface directly"""
        try:
            # Start the HTTP server for SMILES communication
            if not hasattr(self, 'drawing_server') or self.drawing_server is None:
                self.drawing_server = SmilesDrawingServer(port=8001)
                self.drawing_server.start_server()
                
                # Start a timer to check for received SMILES
                self.smiles_check_timer = QTimer()
                self.smiles_check_timer.timeout.connect(self.check_for_received_smiles)
                self.smiles_check_timer.start(1000)  # Check every second
            
            # Check if SMILES-Drawer directory exists
            smiles_drawer_path = os.path.join(os.path.dirname(__file__), 'SMILES-Drawer')
            if not os.path.exists(smiles_drawer_path):
                raise FileNotFoundError(f"SMILES-Drawer directory not found at {smiles_drawer_path}")
            
            self.statusBar().showMessage("Launching structure editor...")
            
            # Launch browser with the drawing interface
            url = "http://localhost:8001/compact.html"  # Use compact interface
            
            # Try to launch in a way that creates a separate window
            try:
                # Try with Chrome app mode first
                chrome_paths = [
                    "C:/Program Files/Google/Chrome/Application/chrome.exe",
                    "C:/Program Files (x86)/Google/Chrome/Application/chrome.exe",
                    "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
                    "/usr/bin/google-chrome",
                    "/usr/bin/chromium-browser"
                ]
                
                chrome_found = False
                for chrome_path in chrome_paths:
                    if os.path.exists(chrome_path):
                        subprocess.Popen([
                            chrome_path,
                            f"--app={url}",
                            "--window-size=900,600",
                            "--disable-web-security",
                            "--disable-features=VizDisplayCompositor"
                        ])
                        chrome_found = True
                        break
                
                if not chrome_found:
                    # Fallback to default browser
                    webbrowser.open(url)
                    
            except Exception as e:
                print(f"Error launching Chrome: {e}")
                # Fallback to default browser
                webbrowser.open(url)
            
            self.statusBar().showMessage("Structure editor opened. Draw your structure and click 'Send to GUI'")
                
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to open structure drawer:\n\n{str(e)}")
            self.statusBar().showMessage("Error opening structure drawer")
    
    def check_for_received_smiles(self):
        """Check if new SMILES data has been received from the web interface"""
        if hasattr(self, 'drawing_server') and self.drawing_server:
            smiles = self.drawing_server.get_smiles()
            if smiles:
                # Update the SMILES input field
                self.smiles_input.setText(smiles)
                self.statusBar().showMessage("Structure imported successfully!")
                
                # Auto-update the reaction scheme
                self.update_reaction_scheme()
                
                # Show a brief confirmation
                QMessageBox.information(self, "Structure Imported", 
                                      f"Structure imported successfully!\n\nSMILES: {smiles}")
    
    def closeEvent(self, event):
        """Clean up when main window is closed"""
        # Stop the SMILES check timer
        if hasattr(self, 'smiles_check_timer'):
            self.smiles_check_timer.stop()
        
        # Stop the drawing server
        if hasattr(self, 'drawing_server') and self.drawing_server:
            self.drawing_server.stop_server()
        
        event.accept()
    
    def load_example(self):
        """Load an example reaction"""
        # Example Buchwald-Hartwig amination
        example_smiles = "c1ccc(Br)cc1.CN>>c1ccc(N(C)C)cc1"
        example_type = "Buchwald-Hartwig Amination"
        
        self.smiles_input.setText(example_smiles)
        
        # Find and set the reaction type
        index = self.reaction_type_combo.findText(example_type)
        if index >= 0:
            self.reaction_type_combo.setCurrentIndex(index)
        
        # Update reaction scheme
        self.update_reaction_scheme()
        
        self.statusBar().showMessage("Example Buchwald-Hartwig reaction loaded")
    
    def browse_sample_reactions(self):
        """Open the sample reactions browser dialog"""
        try:
            browser = SampleReactionsBrowser(self)
            
            # Center the dialog relative to the parent window
            if self.isMaximized():
                # If parent is maximized, position dialog in a good visible location
                screen_geometry = self.screen().geometry()
                dialog_width = 1400
                dialog_height = 800
                x = (screen_geometry.width() - dialog_width) // 2
                y = max(50, (screen_geometry.height() - dialog_height) // 2)  # Ensure title bar is visible
                browser.setGeometry(x, y, dialog_width, dialog_height)
            else:
                # Center relative to parent window
                parent_geometry = self.geometry()
                dialog_width = 1400
                dialog_height = 800
                x = parent_geometry.x() + (parent_geometry.width() - dialog_width) // 2
                y = parent_geometry.y() + (parent_geometry.height() - dialog_height) // 2
                # Ensure the dialog is not positioned off-screen
                x = max(50, min(x, self.screen().geometry().width() - dialog_width - 50))
                y = max(50, min(y, self.screen().geometry().height() - dialog_height - 50))
                browser.setGeometry(x, y, dialog_width, dialog_height)
            
            if browser.exec() == QDialog.DialogCode.Accepted:
                selected_reaction = browser.get_selected_reaction()
                if selected_reaction:
                    self.load_selected_sample_reaction(selected_reaction)
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Could not open sample reactions browser:\n{e}")
    
    def load_selected_sample_reaction(self, reaction_text):
        """Load the selected sample reaction into the interface"""
        if not reaction_text or reaction_text.startswith("Select a sample"):
            return
        
        # Extract SMILES from the reaction text
        smiles = self.extract_smiles_from_sample(reaction_text)
        reaction_type = self.extract_reaction_type_from_sample(reaction_text)
        
        if smiles:
            self.smiles_input.setText(smiles)
            
            # Update reaction scheme
            self.update_reaction_scheme()
            
            # Try to set the appropriate reaction type
            if reaction_type:
                # Find matching reaction type in combo box
                for i in range(self.reaction_type_combo.count()):
                    combo_text = self.reaction_type_combo.itemText(i)
                    if reaction_type.lower() in combo_text.lower():
                        self.reaction_type_combo.setCurrentIndex(i)
                        break
            
            self.statusBar().showMessage(f"Sample reaction loaded: {reaction_type}")
            
            # Show brief info about the loaded reaction
            description = self.extract_description_from_sample(reaction_text)
            if description:
                QMessageBox.information(
                    self, 
                    "Sample Reaction Loaded", 
                    f"Loaded: {description}\n\nReaction Type: {reaction_type}\n\nYou can now run predictions on this reaction."
                )
        else:
            QMessageBox.warning(self, "Error", "Could not extract SMILES from the selected reaction.")
    
    def extract_smiles_from_sample(self, reaction_text):
        """Extract SMILES string from sample reaction text"""
        if ">>" in reaction_text:
            # Split by '>>' and take everything before any parentheses in the product part
            parts = reaction_text.split(">>", 1)  # Only split on first occurrence
            if len(parts) >= 2:
                reactants = parts[0].strip()
                products_and_desc = parts[1].strip()
                
                # Remove description part - look for " (" which indicates start of description
                if " (" in products_and_desc:
                    products = products_and_desc.split(" (", 1)[0].strip()
                else:
                    # If no description with " (", the whole thing is products
                    products = products_and_desc
                
                return f"{reactants}>>{products}"
        
        return reaction_text.strip()
    
    def extract_reaction_type_from_sample(self, reaction_text):
        """Extract reaction type from sample reaction description"""
        if "(" in reaction_text:
            description = reaction_text.split("(", 1)[1]
            if ")" in description:
                description = description.split(")")[0]
            
            # Map common descriptions to reaction types
            desc_lower = description.lower()
            
            if "suzuki" in desc_lower:
                return "C-C Coupling - Suzuki-Miyaura"
            elif "stille" in desc_lower:
                return "C-C Coupling - Stille"
            elif "sonogashira" in desc_lower:
                return "C-C Coupling - Sonogashira"
            elif "heck" in desc_lower:
                return "C-C Coupling - Heck"
            elif "negishi" in desc_lower:
                return "C-C Coupling - Negishi"
            elif "buchwald" in desc_lower or "hartwig" in desc_lower or "b-h" in desc_lower:
                return "C-N Coupling - Buchwald-Hartwig"
            elif "chan-lam" in desc_lower:
                return "C-N Oxidative Coupling - Chan-Lam"
            elif "ullmann ether" in desc_lower:
                return "C-O Coupling - Ullmann"
            elif "ullmann" in desc_lower:
                return "C-N Coupling - Ullmann"
            elif "esterification" in desc_lower:
                return "Esterification"
            elif "amidation" in desc_lower:
                return "Amidation"
            elif "hydrogenation" in desc_lower:
                return "Hydrogenation"
            elif "oxidation" in desc_lower:
                return "Oxidation"
            elif "reduction" in desc_lower or "nabh4" in desc_lower or "lialh4" in desc_lower:
                return "Reduction"
            elif "sn1" in desc_lower:
                return "SN1 - Nucleophilic Substitution (1st order)"
            elif "sn2" in desc_lower:
                return "SN2 - Nucleophilic Substitution (2nd order)"
            elif "e1" in desc_lower:
                return "E1 - Elimination (1st order)"
            elif "e2" in desc_lower:
                return "E2 - Elimination (2nd order)"
            elif "diels-alder" in desc_lower:
                return "Diels-Alder Cycloaddition"
            elif "click" in desc_lower:
                return "Click Chemistry"
        
        return "Auto-detect"
    
    def extract_description_from_sample(self, reaction_text):
        """Extract human-readable description from sample reaction text"""
        if "(" in reaction_text:
            description = reaction_text.split("(", 1)[1]
            if ")" in description:
                description = description.split(")")[0]
            return description
        return ""

def main():
    """Main function to run the application"""
    app = QApplication([])
    
    # Set application properties
    app.setApplicationName("Simple Reaction Predictor")
    app.setApplicationVersion("1.0")
    
    # Create and show the main window
    window = SimpleReactionGUI()
    window.show()
    
    # Run the application
    app.exec()

if __name__ == "__main__":
    main()
