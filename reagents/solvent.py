import numpy as np
import pandas as pd
import os
import json

# Optional dependencies (graceful degradation if missing)
try:  # SciPy distances
    from scipy.spatial.distance import cdist  # type: ignore
except Exception:
    cdist = None  # type: ignore

try:  # SciPy clustering
    from scipy.cluster.hierarchy import dendrogram, linkage  # type: ignore
except Exception:
    dendrogram = None  # type: ignore
    linkage = None  # type: ignore

try:  # Matplotlib
    import matplotlib.pyplot as plt  # type: ignore
except Exception:
    plt = None  # type: ignore

try:  # NetworkX
    import networkx as nx  # type: ignore
except Exception:
    nx = None  # type: ignore

try:  # scikit-learn scaler
    from sklearn.preprocessing import MinMaxScaler  # type: ignore
except Exception:
    MinMaxScaler = None  # type: ignore

# Verify data lengths
solvents = [
    "Acetic Acid",
    "Acetone",
    "Acetonitrile",
    "Benzene",
    "1-Butanol",
    "2-Butanol",
    "tert-Butanol",
    "Carbon Tetrachloride",
    "Chloroform",
    "Cyclohexane",
    "Dichloromethane",
    "Diethyl Ether",
    "Dimethylformamide",
    "Dimethyl Sulfoxide",
    "Ethanol",
    "Ethyl Acetate",
    "Hexane",
    "Methanol",
    "1-Propanol",
    "2-Propanol",
    "Pyridine",
    "THF",
    "Toluene",
    "Water",
    "m-Xylene",
    "o-Xylene",
    "p-Xylene",
    "1,2-Dichloroethane",
    "1,4-Dioxane",
    "Chlorobenzene",
    "Cyclohexanone",
    "Dibutyl Ether",
    "Diisopropyl Ether",
    "Ethylene Glycol",
    "Heptane",
    "Isobutanol",
    "Isopropyl Acetate",
    "Methyl Ethyl Ketone",
    "N-Methyl-2-pyrrolidone",
    "Nitromethane",
    "Pentane",
    "Propionic Acid",
    "Tetrachloroethylene",
    "Triethylamine",
    "2-Methoxyethanol",
    "Butyl Acetate",
    "Diethylene Glycol",
    "Hexafluoroisopropanol",
    "MeTHF",
    "Cyclopentyl methyl ether",
    "Methyl isobutyl ketone",
    "t-Amyl-OH",
]

# Number of solvents
n_solvents = len(solvents)
print(f"Number of solvents: {n_solvents}")


# Verify all property arrays have the correct length
def verify_length(arr, name):
    if len(arr) != n_solvents:
        print(f"Warning: {name} has length {len(arr)}, expected {n_solvents}")
    return arr


# Create expanded solvent properties database
solvent_data = pd.DataFrame(
    {
        "Solvent": solvents,
        "CAS Number": [
            "64-19-7",
            "67-64-1",
            "75-05-8",
            "71-43-2",
            "71-36-3",
            "78-92-2",
            "75-65-0",
            "56-23-5",
            "67-66-3",
            "110-82-7",
            "75-09-2",
            "60-29-7",
            "68-12-2",
            "67-68-5",
            "64-17-5",
            "141-78-6",
            "110-54-3",
            "67-56-1",
            "71-23-8",
            "67-63-0",
            "110-86-1",
            "109-99-9",
            "108-88-3",
            "7732-18-5",
            "108-38-3",
            "95-47-6",
            "106-42-3",
            "107-06-2",
            "123-91-1",
            "108-90-7",
            "108-94-1",
            "142-96-1",
            "108-20-3",
            "107-21-1",
            "142-82-5",
            "78-83-1",
            "108-21-4",
            "78-93-3",
            "872-50-4",
            "75-52-5",
            "109-66-0",
            "79-09-4",
            "127-18-4",
            "121-44-8",
            "109-86-4",
            "123-86-4",
            "111-46-6",
            "920-66-1",
            "96-47-9",
            "5614-37-9",
            "108-10-1",
            "75-85-4",  # New solvents
        ],
        "Abbreviation": [
            "AcOH",
            "Me2CO",
            "MeCN",
            "PhH",
            "BuOH",
            "sBuOH",
            "tBuOH",
            "CCl4",
            "CHCl3",
            "cHex",
            "DCM",
            "Et2O",
            "DMF",
            "DMSO",
            "EtOH",
            "EtOAc",
            "Hex",
            "MeOH",
            "nPrOH",
            "iPrOH",
            "Py",
            "THF",
            "PhMe",
            "H2O",
            "m-Xyl",
            "o-Xyl",
            "p-Xyl",
            "DCE",
            "Diox",
            "PhCl",
            "cHexone",
            "Bu2O",
            "iPr2O",
            "EG",
            "Hept",
            "iBuOH",
            "iPrOAc",
            "MEK",
            "NMP",
            "MeNO2",
            "Pent",
            "PrCOOH",
            "PCE",
            "NEt3",
            "MeOEtOH",
            "BuOAc",
            "DEG",
            "HFIP",
            "MeTHF",
            "CPME",
            "MIBK",
            "tAmOH",  # New solvents
        ],
        "Dielectric Constant": verify_length(
            [
                6.2,
                21.01,
                36.64,
                2.28,
                17.8,
                16.1,
                12.5,
                2.24,
                4.81,
                2.02,
                8.93,
                4.33,
                36.7,
                47.2,
                24.3,
                6.02,
                1.88,
                32.6,
                20.1,
                18.3,
                12.3,
                7.6,
                2.38,
                80.1,
                2.27,
                2.57,
                2.27,
                10.36,
                2.21,
                5.62,
                18.3,
                3.08,
                3.88,
                37.7,
                1.92,
                17.93,
                18.3,
                18.5,
                32.2,
                35.87,
                1.84,
                3.44,
                2.5,
                2.44,
                16.93,
                5.01,
                31.69,
                16.7,
                7.0,
                4.76,
                13.1,
                12.0,  # Values for new solvents
            ],
            "Dielectric Constant",
        ),
        "Polarity Index": verify_length(
            [
                6.2,
                5.1,
                5.8,
                2.7,
                3.9,
                3.9,
                4.0,
                1.6,
                4.1,
                0.2,
                3.1,
                2.8,
                6.4,
                7.2,
                4.3,
                4.4,
                0.1,
                5.1,
                4.0,
                3.9,
                5.3,
                4.0,
                2.4,
                9.0,
                2.5,
                2.5,
                2.5,
                3.5,
                4.8,
                2.7,
                4.5,
                2.9,
                2.4,
                6.9,
                0.2,
                3.9,
                4.3,
                4.7,
                6.7,
                6.7,
                0.0,
                3.0,
                2.9,
                2.9,
                5.5,
                4.0,
                6.0,
                6.6,
                4.0,
                3.4,
                4.2,
                4.1,  # Values for new solvents
            ],
            "Polarity Index",
        ),
        "Boiling Point (°C)": verify_length(
            [
                118,
                56.05,
                81.65,
                80.1,
                117.7,
                99.5,
                82.2,
                76.72,
                61.2,
                80.7,
                39.7,
                34.5,
                153,
                189,
                78.37,
                77.1,
                68.7,
                64.7,
                97.2,
                82.6,
                115.2,
                66,
                110.6,
                100,
                139,
                144.4,
                138.4,
                83.5,
                101.1,
                131.7,
                155.6,
                141,
                68.5,
                197.3,
                98.4,
                108,
                85,
                79.6,
                202,
                101.2,
                36.1,
                141,
                121.2,
                89.5,
                124.5,
                126.1,
                245,
                58.2,
                80.2,
                106,
                116.5,
                102.0,  # Values for new solvents
            ],
            "Boiling Point (°C)",
        ),
        "Density (g/mL)": verify_length(
            [
                1.0446,
                0.7845,
                0.7857,
                0.8765,
                0.8095,
                0.808,
                0.786,
                1.594,
                1.4892,
                0.7781,
                1.3255,
                0.7134,
                0.9445,
                1.1004,
                0.7893,
                0.897,
                0.6548,
                0.7918,
                0.8039,
                0.7855,
                0.9819,
                0.8892,
                0.8669,
                0.9982,
                0.864,
                0.8802,
                0.8611,
                1.253,
                1.033,
                1.106,
                0.947,
                0.764,
                0.72,
                1.113,
                0.684,
                0.802,
                0.872,
                0.805,
                1.028,
                1.138,
                0.626,
                0.993,
                1.622,
                0.726,
                0.965,
                0.882,
                1.118,
                1.596,
                0.854,
                0.860,
                0.802,
                0.810,  # Values for new solvents
            ],
            "Density (g/mL)",
        ),
        "Dipole Moment (D)": verify_length(
            [
                1.74,
                2.69,
                3.92,
                0,
                1.66,
                1.63,
                1.67,
                0,
                1.04,
                0,
                1.60,
                1.15,
                3.82,
                3.96,
                1.69,
                1.78,
                0.08,
                1.70,
                1.68,
                1.66,
                2.19,
                1.75,
                0.31,
                1.85,
                0.37,
                0.45,
                0.3,
                1.86,
                0.45,
                1.69,
                2.87,
                1.18,
                1.13,
                2.28,
                0,
                1.64,
                1.8,
                2.76,
                4.09,
                3.46,
                0,
                1.75,
                0,
                0.87,
                2.18,
                1.84,
                2.31,
                2.05,
                1.60,
                1.23,
                2.79,
                1.70,  # Values for new solvents
            ],
            "Dipole Moment (D)",
        ),
        "Donor Number (DN)": verify_length(
            [
                12.6,
                17.0,
                14.1,
                0.1,
                23.0,
                22.0,
                21.5,
                0.0,
                4.0,
                0.0,
                1.0,
                19.2,
                26.6,
                29.8,
                19.6,
                17.1,
                0.0,
                19.0,
                20.0,
                21.1,
                33.1,
                20.0,
                0.1,
                18.0,
                0.0,
                0.0,
                0.0,
                0.0,
                14.8,
                3.3,
                16.1,
                20.2,
                19.5,
                19.0,
                0.0,
                21.5,
                16.8,
                17.4,
                27.3,
                2.7,
                0.0,
                12.8,
                0.0,
                61.0,
                19.0,
                16.9,
                20.0,
                16.0,
                18.0,
                16.5,
                14.3,
                20.5,  # Values for new solvents
            ],
            "Donor Number (DN)",
        ),
        "Hydrogen Bond Donor": verify_length(
            [
                1.12,
                0.08,
                0.19,
                0.00,
                0.84,
                0.83,
                0.68,
                0.00,
                0.20,
                0.00,
                0.13,
                0.00,
                0.00,
                0.00,
                0.86,
                0.00,
                0.00,
                0.98,
                0.84,
                0.76,
                0.00,
                0.00,
                0.00,
                1.17,
                0.00,
                0.00,
                0.00,
                0.10,
                0.00,
                0.00,
                0.00,
                0.00,
                0.00,
                0.90,
                0.00,
                0.79,
                0.00,
                0.00,
                0.00,
                0.22,
                0.00,
                1.12,
                0.00,
                0.00,
                0.90,
                0.00,
                0.91,
                1.96,
                0.00,
                0.00,
                0.00,
                0.65,  # Values for new solvents
            ],
            "Hydrogen Bond Donor",
        ),
        "Reaction_Compatibility": [
            # Compatibility scores for different reaction types (comma-separated string)
            # Format: "Cross-Coupling,Hydrogenation,Metathesis,C-H_Activation,Carbonylation"
            # Values 0-1 (0: incompatible, 1: highly compatible)
            # Standard Solvents
            "0.4,0.3,0.6,0.5,0.9",
            "0.7,0.4,0.3,0.6,0.8",
            "0.8,0.7,0.3,0.7,0.8",  # AcOH, Acetone, MeCN
            "0.7,0.5,0.2,0.8,0.7",
            "0.5,0.8,0.2,0.5,0.6",
            "0.5,0.8,0.2,0.5,0.6",  # Benzene, 1-BuOH, 2-BuOH
            "0.5,0.8,0.2,0.5,0.6",
            "0.3,0.3,0.8,0.4,0.4",
            "0.6,0.3,0.5,0.6,0.7",  # tBuOH, CCl4, CHCl3
            "0.4,0.6,0.6,0.5,0.5",
            "0.8,0.4,0.3,0.7,0.8",
            "0.7,0.5,0.5,0.6,0.6",  # Cyclohexane, DCM, Et2O
            "0.9,0.7,0.2,0.8,0.9",
            "0.8,0.3,0.2,0.9,0.8",
            "0.6,0.9,0.3,0.6,0.7",  # DMF, DMSO, EtOH
            "0.7,0.4,0.4,0.6,0.8",
            "0.3,0.5,0.6,0.4,0.4",
            "0.5,0.9,0.3,0.5,0.6",  # EtOAc, Hexane, MeOH
            "0.5,0.9,0.2,0.5,0.6",
            "0.6,0.8,0.2,0.6,0.6",
            "0.7,0.6,0.2,0.8,0.7",  # 1-PrOH, 2-PrOH, Pyridine
            "0.9,0.6,0.3,0.8,0.8",
            "0.7,0.5,0.2,0.8,0.7",
            "0.2,0.7,0.1,0.3,0.4",  # THF, Toluene, Water
            "0.7,0.5,0.2,0.8,0.7",
            "0.7,0.5,0.2,0.8,0.7",
            "0.7,0.5,0.2,0.8,0.7",  # Xylenes
            "0.7,0.4,0.3,0.6,0.7",
            "0.8,0.6,0.3,0.7,0.8",
            "0.6,0.5,0.2,0.7,0.6",  # DCE, Dioxane, PhCl
            "0.6,0.3,0.4,0.5,0.7",
            "0.6,0.5,0.5,0.5,0.6",
            "0.7,0.5,0.5,0.6,0.6",  # Cyclohexanone, Bu2O, iPr2O
            "0.4,0.8,0.2,0.4,0.5",
            "0.3,0.5,0.6,0.4,0.4",
            "0.5,0.8,0.2,0.5,0.6",  # EG, Heptane, iBuOH
            "0.7,0.4,0.4,0.6,0.8",
            "0.7,0.3,0.4,0.6,0.8",
            "0.9,0.6,0.2,0.8,0.9",  # iPrOAc, MEK, NMP
            "0.6,0.3,0.2,0.7,0.7",
            "0.3,0.5,0.6,0.4,0.4",
            "0.4,0.3,0.6,0.5,0.8",  # MeNO2, Pentane, PrCOOH
            "0.4,0.3,0.7,0.5,0.5",
            "0.8,0.2,0.1,0.9,0.6",
            "0.5,0.8,0.2,0.5,0.6",  # PCE, NEt3, MeOEtOH
            "0.7,0.4,0.4,0.6,0.8",
            "0.5,0.7,0.2,0.5,0.6",
            "0.3,0.2,0.1,0.4,0.5",  # BuOAc, DEG, HFIP
            "0.9,0.6,0.3,0.8,0.8",
            "0.8,0.5,0.4,0.7,0.7",
            "0.7,0.3,0.4,0.6,0.8",  # MeTHF, CPME, MIBK
            "0.5,0.8,0.2,0.5,0.6",  # tAmOH
        ],
        "Typical_Applications": [
            # Typical applications for each solvent
            "Carbonylation, esterification",
            "Cross-coupling, extraction",
            "Cross-coupling, C-H activation",
            "C-H activation, aromatic chemistry",
            "Hydrogenation, extraction",
            "Hydrogenation, extraction",
            "Hydrogenation, bulky substrates",
            "Metathesis, chlorination",
            "Extraction, cross-coupling",
            "Non-polar reactions, extraction",
            "Cross-coupling, extraction",
            "Extraction, hydrogenation",
            "Cross-coupling, amide formation",
            "Oxidation, polar reactions",
            "Hydrogenation, reduction",
            "Extraction, esterification",
            "Non-polar extraction",
            "Hydrogenation, reduction",
            "Hydrogenation, reduction",
            "Hydrogenation, reduction",
            "Cross-coupling, base reactions",
            "Cross-coupling, organometallic",
            "C-H activation, aromatic",
            "Hydrophilic reactions",
            "C-H activation, aromatic",
            "C-H activation, aromatic",
            "C-H activation, aromatic",
            "Cross-coupling, extraction",
            "Cross-coupling, coordination",
            "C-H activation, electrophilic",
            "Ketone chemistry, oxidation",
            "Extraction, non-polar",
            "Extraction, organometallic",
            "Polymerization, high bp",
            "Extraction, non-polar",
            "Hydrogenation, reduction",
            "Esterification, extraction",
            "Ketone chemistry, oxidation",
            "Polar aprotic reactions",
            "Nitration, polar reactions",
            "Non-polar extraction",
            "Esterification, carbonylation",
            "Metathesis, dehydration",
            "Base reactions, amination",
            "Polyol chemistry, reduction",
            "Esterification, extraction",
            "High bp polar reactions",
            "Fluorinated chemistry",
            "Green THF alternative",
            "Green ether alternative",
            "Ketone chemistry, extraction",
            "Bulky substrate reduction",
        ],
    }
)


def create_solvent_dataframe(json_path: str | None = None, json_dir: str | None = None, prefer_builtin: bool = False):
    """
    Create the solvent DataFrame.

    If a JSON file or directory is present, load from there; otherwise fall back to the
    built-in table. Supports either:
      - data/solvents.json (array of solvent objects), or
      - data/solvents/ (directory of *.json files, one solvent per file)
    """
    # Resolve defaults relative to repo root (../data from this file)
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, 'data'))
    default_file = os.path.join(base_dir, 'solvents.json')
    default_dir = os.path.join(base_dir, 'solvents')

    json_path = json_path or default_file
    json_dir = json_dir or default_dir

    # If explicitly requested, return built-in table
    if prefer_builtin:
        return solvent_data

    def _normalize_entry(entry: dict) -> dict:
        name = entry.get('solvent') or entry.get('name') or entry.get('Solvent')
        abbr = entry.get('abbreviation') or entry.get('Abbreviation')
        cas = entry.get('cas') or entry.get('CAS Number')
        rc = entry.get('reaction_compatibility') or entry.get('Reaction_Compatibility')
        if isinstance(rc, dict):
            order = ["Cross-Coupling", "Hydrogenation", "Metathesis", "C-H_Activation", "Carbonylation"]
            rc_str = ",".join(str(float(rc.get(k, 0.5))) for k in order)
        elif isinstance(rc, (list, tuple)):
            rc_str = ",".join(str(float(x)) for x in rc)
        else:
            rc_str = rc if isinstance(rc, str) else "0.5,0.5,0.5,0.5,0.5"

        apps = entry.get('typical_applications') or entry.get('Typical_Applications') or ''
        if isinstance(apps, (list, tuple)):
            apps_str = ", ".join(map(str, apps))
        else:
            apps_str = str(apps)

        return {
            "Solvent": name,
            "CAS Number": cas,
            "Abbreviation": abbr,
            "Dielectric Constant": entry.get('dielectric_constant') or entry.get('Dielectric Constant'),
            "Polarity Index": entry.get('polarity_index') or entry.get('Polarity Index'),
            "Boiling Point (°C)": entry.get('boiling_point_c') or entry.get('Boiling Point (°C)'),
            "Density (g/mL)": entry.get('density_g_ml') or entry.get('Density (g/mL)'),
            "Dipole Moment (D)": entry.get('dipole_moment_d') or entry.get('Dipole Moment (D)'),
            "Donor Number (DN)": entry.get('donor_number_dn') or entry.get('Donor Number (DN)'),
            "Hydrogen Bond Donor": entry.get('hydrogen_bond_donor') or entry.get('Hydrogen Bond Donor'),
            "Reaction_Compatibility": rc_str,
            "Typical_Applications": apps_str,
        }

    records: list[dict] = []
    try:
        if os.path.isfile(json_path):
            with open(json_path, 'r', encoding='utf-8') as f:
                data = json.load(f)
            if isinstance(data, dict) and 'solvents' in data:
                data = data['solvents']
            if isinstance(data, list):
                records = [_normalize_entry(d) for d in data]
        elif os.path.isdir(json_dir):
            for fn in os.listdir(json_dir):
                if fn.lower().endswith('.json'):
                    with open(os.path.join(json_dir, fn), 'r', encoding='utf-8') as f:
                        d = json.load(f)
                    if isinstance(d, dict) and 'solvent' in d and isinstance(d['solvent'], dict):
                        d = d['solvent']
                    records.append(_normalize_entry(d))
    except Exception as e:
        print(f"Warning: failed to load solvents JSON: {e}")
        records = []

    if records:
        df = pd.DataFrame.from_records(records)
        df = df[df['Solvent'].notna() & (df['Solvent'].astype(str).str.len() > 0)]
        return df.reset_index(drop=True)

    return solvent_data


def create_solvent_feature_matrix():
    """Create a normalized feature matrix for machine learning operations"""
    df = create_solvent_dataframe()

    # Select numerical features (using actual column names from DataFrame)
    feature_columns = [
        "Dielectric Constant",
        "Polarity Index",
        "Boiling Point (°C)",
        "Density (g/mL)",
        "Dipole Moment (D)",
        "Donor Number (DN)",
        "Hydrogen Bond Donor",
    ]

    X = []
    for _, row in df.iterrows():
        features = []
        for col in feature_columns:
            if col in df.columns:
                features.append(row[col])
            else:
                features.append(0)  # Default value if column missing
        X.append(features)

    # Normalize features
    X = np.array(X)
    scaler = MinMaxScaler()
    X_normalized = scaler.fit_transform(X)

    return X_normalized


# Reaction-specific property weights for solvents
SOLVENT_REACTION_WEIGHTS = {
    "Cross-Coupling": {
        "Dielectric Constant": 0.15,
        "Polarity Index": 0.20,
        "Boiling Point (°C)": 0.10,
        "Density (g/mL)": 0.05,
        "Dipole Moment (D)": 0.15,
        "Donor Number (DN)": 0.25,
        "Hydrogen Bond Donor": 0.10,
    },
    "Hydrogenation": {
        "Dielectric Constant": 0.10,
        "Polarity Index": 0.15,
        "Boiling Point (°C)": 0.15,
        "Density (g/mL)": 0.05,
        "Dipole Moment (D)": 0.10,
        "Donor Number (DN)": 0.35,
        "Hydrogen Bond Donor": 0.10,
    },
    "Metathesis": {
        "Dielectric Constant": 0.20,
        "Polarity Index": 0.25,
        "Boiling Point (°C)": 0.15,
        "Density (g/mL)": 0.05,
        "Dipole Moment (D)": 0.10,
        "Donor Number (DN)": 0.05,
        "Hydrogen Bond Donor": 0.20,
    },
    "C-H_Activation": {
        "Dielectric Constant": 0.15,
        "Polarity Index": 0.25,
        "Boiling Point (°C)": 0.20,
        "Density (g/mL)": 0.05,
        "Dipole Moment (D)": 0.15,
        "Donor Number (DN)": 0.10,
        "Hydrogen Bond Donor": 0.10,
    },
    "Carbonylation": {
        "Dielectric Constant": 0.15,
        "Polarity Index": 0.20,
        "Boiling Point (°C)": 0.10,
        "Density (g/mL)": 0.05,
        "Dipole Moment (D)": 0.20,
        "Donor Number (DN)": 0.20,
        "Hydrogen Bond Donor": 0.10,
    },
}


def parse_solvent_reaction_compatibility(compatibility_str, reaction_type):
    """Parse reaction compatibility string and return score for specific reaction type"""
    reaction_types = [
        "Cross-Coupling",
        "Hydrogenation",
        "Metathesis",
        "C-H_Activation",
        "Carbonylation",
    ]

    if reaction_type not in reaction_types:
        return 0.5  # Default compatibility for unknown reactions

    try:
        scores = [float(x) for x in compatibility_str.split(",")]
        reaction_idx = reaction_types.index(reaction_type)
        return scores[reaction_idx] if reaction_idx < len(scores) else 0.5
    except:
        return 0.5


def calculate_solvent_weighted_similarity(
    solvent1_features, solvent2_features, weights
):
    """Calculate weighted similarity between two solvents based on reaction-specific weights"""
    similarity = 0.0

    feature_names = [
        "Dielectric Constant",
        "Polarity Index",
        "Boiling Point (°C)",
        "Density (g/mL)",
        "Dipole Moment (D)",
        "Donor Number (DN)",
        "Hydrogen Bond Donor",
    ]

    for i, feature in enumerate(feature_names):
        if feature in weights:
            # Normalized difference (0 = identical, 1 = maximum difference)
            max_val = max(abs(solvent1_features[i]), abs(solvent2_features[i]), 1)
            diff = abs(solvent1_features[i] - solvent2_features[i]) / max_val
            similarity += weights[feature] * (1 - diff)

    return similarity


def recommend_solvents_for_reaction(
    target_solvent=None, reaction_type="Cross-Coupling", top_n=5, min_compatibility=0.3
):
    """
    Recommend solvents for a specific reaction type.

    Parameters:
    - target_solvent: Name of reference solvent (optional)
    - reaction_type: Type of reaction ('Cross-Coupling', 'Hydrogenation', 'Metathesis',
                    'C-H_Activation', 'Carbonylation')
    - top_n: Number of recommendations to return
    - min_compatibility: Minimum compatibility score for the reaction type

    Returns:
    - List of dictionaries with solvent recommendations
    """
    df = create_solvent_dataframe()

    # Filter solvents by reaction compatibility
    compatible_solvents = []
    for idx, row in df.iterrows():
        compatibility = parse_solvent_reaction_compatibility(
            row["Reaction_Compatibility"], reaction_type
        )
        if compatibility >= min_compatibility:
            compatible_solvents.append(
                {
                    "index": idx,
                    "name": row.get("Solvent", row.get("name", "")),
                    "compatibility": compatibility,
                    "applications": row.get("Typical_Applications", ""),
                    "abbreviation": row.get("Abbreviation", ""),
                }
            )

    # Sort by compatibility score
    compatible_solvents.sort(key=lambda x: x["compatibility"], reverse=True)

    # If target solvent is specified, calculate similarity scores
    if target_solvent:
        target_idx = None
        for idx, name in enumerate(df["Solvent"]):
            if name.lower() == target_solvent.lower():
                target_idx = idx
                break

        if target_idx is not None:
            X = create_solvent_feature_matrix()
            target_features = X[target_idx]
            weights = SOLVENT_REACTION_WEIGHTS.get(
                reaction_type, SOLVENT_REACTION_WEIGHTS["Cross-Coupling"]
            )

            # Calculate similarity scores
            for solvent in compatible_solvents:
                if solvent["index"] != target_idx:  # Don't include the target itself
                    solvent_features = X[solvent["index"]]
                    similarity = calculate_solvent_weighted_similarity(
                        target_features, solvent_features, weights
                    )
                    solvent["similarity"] = similarity
                else:
                    solvent["similarity"] = 0  # Remove target from recommendations

            # Sort by similarity score (keeping only compatible solvents)
            compatible_solvents = [
                s for s in compatible_solvents if s["index"] != target_idx
            ]
            compatible_solvents.sort(
                key=lambda x: (x["compatibility"] * 0.6 + x["similarity"] * 0.4),
                reverse=True,
            )

    # Return top recommendations
    recommendations = []
    for i, solvent in enumerate(compatible_solvents[:top_n]):
        rec = {
            "rank": i + 1,
            "solvent": solvent["name"],
            "abbreviation": solvent["abbreviation"],
            "compatibility_score": round(solvent["compatibility"], 3),
            "applications": solvent["applications"],
            "reaction_suitability": reaction_type,
        }

        if "similarity" in solvent:
            rec["similarity_score"] = round(solvent["similarity"], 3)
            rec["combined_score"] = round(
                solvent["compatibility"] * 0.6 + solvent["similarity"] * 0.4, 3
            )

        recommendations.append(rec)

    return recommendations


def get_reaction_specific_solvents(reaction_type, property_preferences=None):
    """
    Get solvents optimized for a specific reaction type with optional property preferences.

    Parameters:
    - reaction_type: Type of reaction
    - property_preferences: Dict with preferred ranges for properties

    Returns:
    - Filtered solvent recommendations
    """
    recommendations = recommend_solvents_for_reaction(
        reaction_type=reaction_type, top_n=10, min_compatibility=0.4
    )

    if property_preferences:
        df = create_solvent_dataframe()
        filtered_recs = []

        for rec in recommendations:
            solvent_name = rec["solvent"]
            solvent_idx = df[df["Solvent"] == solvent_name].index[0]

            # Check property preferences
            meets_criteria = True
            if "bp_max" in property_preferences:
                if (
                    df.loc[solvent_idx, "Boiling Point (°C)"]
                    > property_preferences["bp_max"]
                ):
                    meets_criteria = False

            if "bp_min" in property_preferences:
                if (
                    df.loc[solvent_idx, "Boiling Point (°C)"]
                    < property_preferences["bp_min"]
                ):
                    meets_criteria = False

            if "polarity_max" in property_preferences:
                if (
                    df.loc[solvent_idx, "Polarity Index"]
                    > property_preferences["polarity_max"]
                ):
                    meets_criteria = False

            if "polarity_min" in property_preferences:
                if (
                    df.loc[solvent_idx, "Polarity Index"]
                    < property_preferences["polarity_min"]
                ):
                    meets_criteria = False

            if "protic" in property_preferences:
                is_protic = df.loc[solvent_idx, "Hydrogen Bond Donor"] > 0.5
                if property_preferences["protic"] != is_protic:
                    meets_criteria = False

            if meets_criteria:
                filtered_recs.append(rec)

        return filtered_recs[:5]  # Return top 5 that meet criteria

    return recommendations


# Define property weights
property_weights = {
    "CAS Number": 0.00,  # Zero weight - not used for similarity
    "Abbreviation": 0.00,  # Zero weight - not used for similarity
    "Dielectric Constant": 0.08,
    "Polarity Index": 0.27,
    "Boiling Point (°C)": 0.02,
    "Density (g/mL)": 0.02,
    "Dipole Moment (D)": 0.11,
    "Donor Number (DN)": 0.20,
    "Hydrogen Bond Donor": 0.30,
}

# Export solvent database to Excel with the new properties
excel_path = "solvent_database.xlsx"
with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
    # Main data sheet
    solvent_data.to_excel(writer, sheet_name="Solvent Properties", index=False)

    # Additional information sheet with aligned arrays
    properties = ["CAS Number", "Abbreviation"] + list(property_weights.keys())
    weights = [0.00, 0.00] + list(property_weights.values())
    descriptions = [
        "Chemical Abstracts Service registry number",
        "Common abbreviation used in literature",
        "Relative permittivity - measure of polarity",
        "Empirical measure of solvent polarity",
        "Temperature at normal pressure where liquid->gas",
        "Mass per unit volume at 20°C",
        "Measure of charge separation",
        "Lewis basicity measure",
        "Ability to donate hydrogen bonds",
    ]

    # Ensure all lists have the same length
    while len(descriptions) < len(properties):
        descriptions.append("")

    info_data = pd.DataFrame(
        {
            "Property": properties,
            "Weight": weights,
            "Description": descriptions[: len(properties)],  # Trim to match length
        }
    )
    info_data.to_excel(writer, sheet_name="Property Information", index=False)

print(f"Solvent database has been exported to {excel_path}")


# Function to recommend similar solvents
def recommend_solvents(selected_solvent, num_recommendations=3):
    if selected_solvent not in solvent_data["Solvent"].values:
        print(f"Solvent {selected_solvent} not found in the database.")
        return []

    # Get only numeric columns for normalization (exclude 'Solvent', 'CAS Number', and 'Abbreviation')
    numeric_columns = [
        col
        for col in property_weights.keys()
        if col not in ["CAS Number", "Abbreviation"]
    ]

    # Normalize the properties to 0-1 scale
    properties_normalized = solvent_data.copy()
    for column in numeric_columns:
        min_val = properties_normalized[column].min()
        max_val = properties_normalized[column].max()
        properties_normalized[column] = (properties_normalized[column] - min_val) / (
            max_val - min_val
        )

    # Apply weights to properties
    properties_weighted = properties_normalized.copy()
    for column in numeric_columns:
        properties_weighted[column] = (
            properties_weighted[column] * property_weights[column]
        )

    # Extract weighted properties for similarity calculation
    solvent_properties = properties_weighted[numeric_columns].values
    solvent_names = solvent_data["Solvent"].values

    # Find the index of the selected solvent
    selected_idx = np.where(solvent_names == selected_solvent)[0][0]
    selected_vector = solvent_properties[selected_idx].reshape(1, -1)

    # Compute weighted Euclidean distances
    distances = cdist(selected_vector, solvent_properties, metric="euclidean")[0]

    # Sort and get the closest solvents (excluding itself)
    closest_indices = np.argsort(distances)[1 : num_recommendations + 1]
    recommended_solvents = solvent_names[closest_indices]

    # Add similarity scores
    similarity_scores = 1 / (1 + distances[closest_indices])
    recommendations = list(zip(recommended_solvents, similarity_scores))

    return recommendations


def analyze_solvent_clusters():
    # Get only numeric columns for normalization
    numeric_columns = [
        col
        for col in property_weights.keys()
        if col not in ["CAS Number", "Abbreviation"]
    ]

    # Normalize all properties for clustering
    properties_normalized = solvent_data.copy()
    for column in numeric_columns:
        min_val = properties_normalized[column].min()
        max_val = properties_normalized[column].max()
        properties_normalized[column] = (properties_normalized[column] - min_val) / (
            max_val - min_val
        )

    # Apply weights to properties
    properties_weighted = properties_normalized.copy()
    for column in numeric_columns:
        properties_weighted[column] = (
            properties_weighted[column] * property_weights[column]
        )

    # Prepare data for clustering
    X = properties_weighted[numeric_columns].values

    # Rest of clustering code remains the same
    Z = linkage(X, method="ward")
    plt.figure(figsize=(15, 10))
    plt.title("Hierarchical Clustering of Solvents")
    dendrogram(
        Z, labels=solvent_data["Solvent"].values, leaf_rotation=90, leaf_font_size=8
    )
    plt.xlabel("Solvents")
    plt.ylabel("Distance")
    plt.tight_layout()
    plt.savefig("solvent_clusters.png", dpi=300, bbox_inches="tight")
    plt.close()
    print(
        "Clustering analysis completed. The dendrogram has been saved as 'solvent_clusters.png'"
    )


def create_solvent_network(threshold=0.7):
    # Get only numeric columns
    numeric_columns = [
        col
        for col in property_weights.keys()
        if col not in ["CAS Number", "Abbreviation"]
    ]

    # Normalize all properties
    properties_normalized = solvent_data[numeric_columns].copy()
    scaler = MinMaxScaler()
    properties = scaler.fit_transform(properties_normalized)

    # Rest of network code remains the same
    similarity_matrix = 1 / (1 + cdist(properties, properties, metric="euclidean"))
    G = nx.Graph()

    # Add nodes (solvents)
    for solvent in solvent_data["Solvent"]:
        G.add_node(solvent)

    # Add edges based on similarity threshold
    for i in range(len(solvent_data)):
        for j in range(i + 1, len(solvent_data)):
            if similarity_matrix[i, j] > threshold:
                G.add_edge(
                    solvent_data["Solvent"].iloc[i],
                    solvent_data["Solvent"].iloc[j],
                    weight=similarity_matrix[i, j],
                )

    # Set up the plot
    plt.figure(figsize=(20, 20))

    # Create layout
    pos = nx.spring_layout(G, k=1, iterations=50)

    # Draw network
    nx.draw_networkx_nodes(G, pos, node_size=1000, node_color="lightblue")
    nx.draw_networkx_edges(G, pos, alpha=0.2)
    nx.draw_networkx_labels(G, pos, font_size=8)

    plt.title("Solvent Similarity Network\n(Edges show similarities above threshold)")
    plt.axis("off")

    # Save the network plot
    plt.savefig("solvent_network.png", dpi=300, bbox_inches="tight")
    plt.close()

    print(
        "Network analysis completed. The graph has been saved as 'solvent_network.png'"
    )
    return G


if __name__ == "__main__":
    # Example Usage
    selected_solvent = "Acetone"
    recommended = recommend_solvents(selected_solvent, num_recommendations=5)
    print(f"\nRecommendations for {selected_solvent}:")
    for solvent, score in recommended:
        print(f"{solvent}: {score:.2f} similarity score")

    # Perform clustering analysis
    analyze_solvent_clusters()
    create_solvent_network(threshold=0.7)  # Adjust threshold to control network density
