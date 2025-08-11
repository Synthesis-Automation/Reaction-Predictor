# 🧪 How to Use the Enhanced Reaction Condition Recommender

## Overview
Your reaction predictor now includes an **enhanced ligand and solvent recommendation system** that provides intelligent, reaction-specific condition suggestions for major organometallic and organic reactions.

## ✨ New Features Added

### 1. **Comprehensive Reagent Database**
- **130+ ligands** across all major categories
- **52+ solvents** with detailed properties  
- **Reaction compatibility scoring** for each reagent

### 2. **Intelligent Recommendations**
- **Reaction type detection** from SMILES patterns
- **Combined condition scoring** (ligand + solvent synergy)
- **Property-based filtering** (cost, boiling point, etc.)
- **Confidence scoring** for reliability assessment

### 3. **Supported Reaction Types**
- ✅ **Cross-Coupling** (Suzuki, Buchwald-Hartwig, Heck, etc.)
- ✅ **Hydrogenation** (asymmetric and achiral)
- ✅ **Metathesis** (olefin metathesis)
- ✅ **C-H Activation** 
- ✅ **Carbonylation**

## 🚀 How to Use

### Step 1: Enter Your Reaction
```
Input: c1ccc(Br)cc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1
Type: Suzuki-Miyaura Coupling
```

### Step 2: Get Enhanced Recommendations
Click **"Predict Conditions"** to receive:

#### 🎯 **Top Combined Recommendations**
```
RECOMMENDATION #1 🟢 High Confidence
🧬 Catalyst System:
  • Ligand: XPhos (Score: 0.9)
  • Solvent: THF (Score: 0.9)  
  • Combined Score: 1.0 (includes +0.1 synergy bonus)

⚙️ Typical Conditions:
  • Temperature: 80-120°C
  • Time: 4-24 hours
  • Atmosphere: Inert (N₂ or Ar)
  • Catalyst Loading: 1-5 mol%
  • Base: K₂CO₃ or Cs₂CO₃
```

#### 🧬 **Individual Reagent Options**
```
TOP LIGAND OPTIONS:
  1. SPhos (Score: 0.9)
     • Applications: Buchwald-Hartwig amination
     • Suitability: Cross-Coupling

  2. XPhos (Score: 0.9)
     • Applications: Suzuki-Miyaura coupling
     • Suitability: Cross-Coupling

TOP SOLVENT OPTIONS:
  1. DMF (Score: 0.9)
     • Applications: Cross-coupling, amide formation
     • Suitability: Cross-Coupling
```

#### 🎛️ **Specialized Options**
```
💰 Budget-Friendly Ligands:
  • PPh3 (Score: 0.8)
  • PMe3 (Score: 0.7)

🌡️ Low-Boiling Solvents (Easy Removal):
  • THF (Score: 0.9)
  • Acetone (Score: 0.7)

🌱 Green/Sustainable Solvents:
  • MeTHF (Score: 0.9)
  • CPME (Score: 0.8)
```

## 💡 **Usage Tips**

### For Best Results:
1. **Start with highest-scoring combinations** - they have the most literature support
2. **Check confidence levels** - 🟢 High > 🟡 Medium > 🔴 Low
3. **Consider your constraints** - use specialized options for budget/sustainability
4. **Read reaction-specific guidance** - each reaction type has optimization tips

### Example Workflow:
1. **Primary choice**: Highest-scoring combination
2. **Backup option**: Second-highest combination  
3. **Budget alternative**: Check budget-friendly options
4. **Green chemistry**: Consider sustainable solvent alternatives
5. **Troubleshooting**: Use guidance notes if reactions don't work

## 🔬 **Sample Test Cases**

### Cross-Coupling Reactions
```
Suzuki: c1ccc(Br)cc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1
Buchwald-Hartwig: c1ccc(Br)cc1.NC>>c1ccc(NC)cc1
```

### Hydrogenation Reactions  
```
Alkene: C=C>>CC
Ketone: CC(=O)C>>CC(O)C
```

### Other Reactions
```
Metathesis: C=C.C=C>>C=C (ring-closing)
C-H Activation: c1ccccc1.Br>>c1ccc(Br)cc1
Carbonylation: c1ccc(Br)cc1>>c1ccc(C=O)cc1
```

## 📊 **Understanding the Scoring System**

### Compatibility Scores:
- **0.9-1.0**: Excellent - Highly optimized for this reaction type
- **0.7-0.8**: Good - Generally reliable 
- **0.5-0.6**: Fair - May work with optimization
- **<0.5**: Poor - Not recommended unless specialized need

### Combined Scores:
- Include **synergy bonuses** for proven ligand-solvent pairs
- Account for **reaction-specific weighting** of properties
- Provide **confidence assessment** based on literature data

## 🛠️ **Troubleshooting**

### If You Get Basic Results:
```
Status: Basic analysis completed (recommendation engines failed)
```
This means the enhanced system isn't available. Check:
1. Are `ligand.py` and `solvent.py` in the `reagents/` folder?
2. Are required Python packages installed? (`pandas`, `numpy`, `scikit-learn`)
3. Try running the test files to verify system functionality

### If Recommendations Seem Off:
1. **Verify your SMILES** - incorrect SMILES lead to wrong reaction type detection
2. **Check reaction type selection** - manual selection overrides auto-detection  
3. **Consider substrate specifics** - the system gives general recommendations, fine-tuning may be needed

## 🎉 **Benefits of the Enhanced System**

1. **Time Saving**: No more literature searching for similar reactions
2. **Comprehensive Coverage**: 5 major reaction types, 130+ ligands, 52+ solvents
3. **Intelligent Scoring**: Accounts for synergy and reaction-specific properties
4. **Practical Guidance**: Includes typical conditions and optimization tips
5. **Flexibility**: Budget-friendly and green chemistry alternatives

---

**Happy Synthesizing! 🧪⚗️**

*For questions or issues, refer to the test files or check the technical documentation.*
