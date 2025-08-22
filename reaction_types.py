"""
Reaction Types Data
==================

This file contains all the reaction type categories used in the GUI
for reaction type selection and specialized analysis.
"""

# Comprehensive reaction type categories for the dropdown
REACTION_TYPES = [
    "Auto detect reaction type",
    "───────── Specific Reaction Types ─────────",
    "C-C Coupling - Suzuki-Miyaura (Pd)",
    "C-C Coupling - Stille (Pd)",
    "C-C Coupling - Sonogashira (Pd)", 
    "C-C Coupling - Heck (Pd)",
    "C-C Coupling - Negishi (Pd)",
    "C-C Coupling - Kumada (Ni)",
    "C-N Coupling - Buchwald-Hartwig (Pd)",
    "C-N Coupling - Ullmann (Cu)",
    "C-N Coupling - Chan-Lam (Cu)",
    "C-O Coupling - Ullmann Ether (Cu)",
    "C-O Coupling - Mitsunobu",
    "C-S Coupling - Thioether Formation (Pd)",
    "───────── Substitution Reactions ─────────",
    "SN1 - Nucleophilic Substitution (1st order)",
    "SN2 - Nucleophilic Substitution (2nd order)", 
    "SNAr - Nucleophilic Aromatic Substitution",
    "SE - Electrophilic Substitution",
    "───────── Elimination Reactions ─────────",
    "E1 - Elimination (1st order)",
    "E2 - Elimination (2nd order)",
    "Dehydration - Alcohol to Alkene",
    "Dehydrohalogenation",
    "───────── Addition Reactions ─────────",
    "Hydrogenation - Alkene/Alkyne",
    "Hydrogenation - Carbonyl", 
    "Hydrogenation - Aromatic",
    "Hydroboration-Oxidation",
    "Oxymercuration-Demercuration",
    "Epoxidation",
    "Dihydroxylation",
    "───────── Carbonyl Chemistry ─────────",
    "Aldol Condensation",
    "Aldol Addition",
    "Wittig Reaction",
    "Mannich Reaction",
    "Claisen Condensation",
    "Michael Addition",
    "Robinson Annulation",
    "───────── Cycloaddition Reactions ─────────",
    "Diels-Alder [4+2]",
    "1,3-Dipolar Cycloaddition",
    "Click Chemistry (CuAAC)",
    "[2+2] Cycloaddition",
    "Pauson-Khand Reaction",
    "───────── Oxidation Reactions ─────────",
    "Alcohol Oxidation (1° → Aldehyde)",
    "Alcohol Oxidation (2° → Ketone)",
    "Swern Oxidation",
    "PCC/PDC Oxidation",
    "TEMPO Oxidation",
    "Epoxidation",
    "Ozonolysis",
    "───────── Reduction Reactions ─────────",
    "NaBH4 - Carbonyl Reduction",
    "LiAlH4 - Strong Reduction",
    "DIBAL-H - Selective Reduction",
    "Catalytic Hydrogenation",
    "Transfer Hydrogenation",
    "Wolff-Kishner Reduction",
    "Clemmensen Reduction",
    "───────── Acylation/Esterification ─────────",
    "Esterification - Acid + Alcohol",
    "Transesterification",
    "Amidation - Acid + Amine",
    "Friedel-Crafts Acylation",
    "Steglich Esterification",
    "───────── Organometallic Reactions ─────────",
    "Grignard Reaction",
    "Organolithium Addition",
    "Organozinc Coupling",
    "Organocopper Reaction",
    "───────── Metathesis Reactions ─────────",
    "Ring-Closing Metathesis (RCM)",
    "Cross-Metathesis (CM)",
    "Ring-Opening Metathesis (ROM)",
    "───────── Rearrangement Reactions ─────────",
    "Claisen Rearrangement",
    "Cope Rearrangement",
    "Pinacol Rearrangement",
    "Beckmann Rearrangement",
    "Curtius Rearrangement",
    "───────── Heterocycle Synthesis ─────────",
    "Pyridine Synthesis",
    "Quinoline Synthesis (Skraup)",
    "Indole Synthesis (Fischer)",
    "Benzimidazole Synthesis",
    "Thiazole Synthesis",
    "Oxazole Synthesis",
    "Triazole Synthesis",
    "───────── Multicomponent Reactions ─────────",
    "Ugi Reaction (4-component)",
    "Passerini Reaction (3-component)",
    "Biginelli Reaction",
    "Strecker Reaction",
    "───────── Protection/Deprotection ─────────",
    "Alcohol Protection",
    "Amine Protection",
    "Carbonyl Protection",
    "Carboxylic Acid Protection",
    "───────── Named Reactions ─────────",
    "Williamson Ether Synthesis",
    "Gabriel Amine Synthesis",
    "Malonic Ester Synthesis", 
    "Acetoacetic Ester Synthesis",
    "Hell-Volhard-Zelinsky",
    "Sandmeyer Reaction",
    "Birch Reduction",
    "Baeyer-Villiger Oxidation"
]

def get_reaction_types():
    """Get the list of all reaction types"""
    return REACTION_TYPES

def get_coupling_reactions():
    """Get only coupling reaction types"""
    return [r for r in REACTION_TYPES if "Coupling" in r]

def get_oxidation_reactions():
    """Get only oxidation reaction types"""
    return [r for r in REACTION_TYPES if "Oxidation" in r]

def get_reduction_reactions():
    """Get only reduction reaction types"""
    return [r for r in REACTION_TYPES if "Reduction" in r or "Hydrogenation" in r]
