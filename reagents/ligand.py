import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from scipy.cluster.hierarchy import dendrogram, linkage
import matplotlib.pyplot as plt
import networkx as nx
from sklearn.preprocessing import MinMaxScaler

# Define expanded ligand database
ligands = [
    # Monodentate Phosphines - Simple
    "PPh3",
    "PCy3",
    "PtBu3",
    "P(o-tol)3",
    "P(p-tol)3",
    "P(p-F-Ph)3",
    "P(p-Cl-Ph)3",
    "P(p-CF3-Ph)3",
    "PMePh2",
    "PMe2Ph",
    "PMe3",
    "PEt3",
    "P(nBu)3",
    "P(iPr)3",
    "P(2-furyl)3",
    "P(C6F5)3",
    "P(p-MeO-Ph)3",
    "P(o-MeO-Ph)3",
    "P(mes)3",
    # Buchwald-type Ligands
    "SPhos",
    "XPhos",
    "RuPhos",
    "BrettPhos",
    "tBuXPhos",
    "JohnPhos",
    "DavePhos",
    "MePhos",
    "AmplPhos",
    "QPhos",
    "CyJohnPhos",
    "AlPhos",
    "Me4tBuXPhos",
    "AdBrettPhos",
    "tBuBrettPhos",
    "Ph-XPhos",
    "MorDalPhos",
    # Bidentate Phosphines - BINAP Family
    "BINAP",
    "Tol-BINAP",
    "H8-BINAP",
    "BIPHEP",
    "MeO-BIPHEP",
    "SEGPHOS",
    "DM-SEGPHOS",
    "DTBM-SEGPHOS",
    "C3-TunePhos",
    "DifluorPhos",
    "SynPhos",
    # Other Bidentate Phosphines
    "DPPE",
    "DPPP",
    "DPPB",
    "DPPF",
    "XantPhos",
    "DPEPhos",
    "JOSIPHOS",
    "DIOP",
    "DUPHOS",
    "TangPhos",
    "BenzP*",
    "FerroTANE",
    "JosiPhos-1",
    "CyPF-tBu",
    "QPhos",
    "BIBOP",
    "MeO-F12-BIPHEP",
    "CTH-JORPHOS",
    # NHC Ligands (removed IPr* and IMes*(2) variants)
    "IPr",
    "IMes",
    "SIPr",
    "SIMes",
    "IPrCl",
    "IPent",
    "IHept",
    "IAd",
    "ICy",
    "ItBu",
    "IBox",
    "IBiox",
    "IDD",
    "IAd*",
    # Nitrogen Ligands
    "TMEDA",
    "Pyridine",
    "Phenanthroline",
    "Bipyridine",
    "4-DMAP",
    "DBU",
    "DABCO",
    "Terpy",
    "Neocuproine",
    "Bathophenanthroline",
    "Bathocuproine",
    "PyBOX",
    "Box",
    "Quinox",
    "TMCDA",
    "PMETA",
    "DiPIC",
    # Phosphoramidites
    "MonoPhos",
    "PipPhos",
    "SIPHOS",
    "QUINAP",
    "TADDOL-based",
    "PhosphorAmidite",
    # Carbene Precursors
    "IMes·HCl",
    "IPr·HCl",
    "SIMes·HCl",
    "SIPr·HCl",
    "ICy·HCl",
    "ItBu·HCl",
    # Other
    "dba",
    "Acac",
    "COD",
    "Cp",
    "Cp*",
    "tBuCp",
    "IndCp",
    "TMS-Cp",
    "Allyl",
    "cod",
    "nbd",
    "PhCN",
    "4-pic",
    "TMHD",
    "hfacac",
]

# Create expanded ligand properties database
ligand_data = pd.DataFrame(
    {
        "Ligand": ligands,
        "Cone Angle (°)": [
            # Monodentate Phosphines
            145,
            170,
            182,
            150,
            145,
            145,
            145,
            148,
            136,
            122,
            118,
            132,
            132,
            160,
            133,
            184,
            145,
            150,
            212,
            # Buchwald-type
            155,
            160,
            158,
            160,
            160,
            155,
            155,
            145,
            158,
            160,
            157,
            158,
            165,
            163,
            162,
            159,
            156,
            # BINAP Family
            165,
            165,
            166,
            160,
            165,
            165,
            168,
            175,
            163,
            164,
            166,
            # Other Bidentate
            125,
            127,
            130,
            125,
            155,
            140,
            170,
            125,
            125,
            128,
            165,
            124,
            172,
            168,
            160,
            158,
            167,
            164,
            # NHC Ligands (removed corresponding values)
            175,
            170,
            175,
            170,
            175,
            178,
            180,
            175,
            168,
            182,
            165,
            167,
            179,
            177,
            # Nitrogen Ligands
            80,
            100,
            95,
            95,
            105,
            110,
            105,
            98,
            96,
            97,
            97,
            102,
            94,
            93,
            82,
            85,
            88,
            # Phosphoramidites
            125,
            128,
            127,
            158,
            135,
            130,
            # Carbene Precursors
            170,
            175,
            170,
            175,
            168,
            182,
            # Other
            90,
            90,
            80,
            85,
            85,
            90,
            88,
            86,
            82,
            80,
            80,
            95,
            98,
            92,
            88,
        ],
        "Electronic Parameter (cm⁻¹)": [
            # Monodentate Phosphines
            2068.9,
            2056.4,
            2056.1,
            2066.7,
            2066.7,
            2071.3,
            2072.8,
            2074.3,
            2067.2,
            2065.3,
            2064.1,
            2061.7,
            2060.3,
            2059.2,
            2068.5,
            2084.3,
            2066.1,
            2066.3,
            2064.8,
            # Buchwald-type
            2065.7,
            2064.8,
            2065.7,
            2061.8,
            2060.3,
            2064.5,
            2063.7,
            2066.2,
            2063.5,
            2062.8,
            2064.7,
            2063.9,
            2061.2,
            2060.5,
            2060.8,
            2065.1,
            2064.2,
            # BINAP Family
            2067.8,
            2067.5,
            2067.9,
            2067.4,
            2066.5,
            2064.5,
            2064.2,
            2063.8,
            2065.7,
            2068.2,
            2066.8,
            # Other Bidentate
            2073.6,
            2073.0,
            2073.3,
            2072.1,
            2066.5,
            2066.8,
            2065.4,
            2071.8,
            2071.5,
            2070.8,
            2069.5,
            2072.4,
            2065.2,
            2064.8,
            2062.8,
            2068.5,
            2067.2,
            2066.9,
            # NHC Ligands (removed corresponding values)
            2051.5,
            2051.2,
            2051.5,
            2051.2,
            2052.0,
            2050.8,
            2050.5,
            2051.8,
            2052.5,
            2053.2,
            2052.8,
            2052.5,
            2051.6,
            2051.9,
            # Nitrogen Ligands
            2083.2,
            2083.2,
            2085.6,
            2085.2,
            2081.5,
            2080.8,
            2082.5,
            2084.8,
            2085.4,
            2085.8,
            2085.7,
            2084.2,
            2084.5,
            2085.1,
            2083.5,
            2083.8,
            2084.1,
            # Phosphoramidites
            2069.5,
            2069.8,
            2069.2,
            2068.5,
            2069.0,
            2069.4,
            # Carbene Precursors
            2051.2,
            2051.5,
            2051.2,
            2051.5,
            2052.5,
            2053.2,
            # Other
            2090.0,
            2090.0,
            2089.5,
            2088.0,
            2087.5,
            2088.2,
            2088.5,
            2088.3,
            2089.2,
            2089.5,
            2089.5,
            2082.5,
            2083.0,
            2089.8,
            2091.2,
        ],
        "Bite Angle (°)": [
            # Monodentate Phosphines (all 0)
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            # Buchwald-type (mostly monodentate)
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            # BINAP Family
            93,
            93,
            93,
            92,
            92,
            95,
            95,
            95,
            94,
            93,
            94,
            # Other Bidentate
            85,
            91,
            98,
            96,
            110,
            102,
            96,
            98,
            99,
            85,
            92,
            95,
            96,
            97,
            102,
            95,
            92,
            94,
            # NHC Ligands
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            # Nitrogen Ligands
            85,
            0,
            82,
            81,
            0,
            0,
            0,
            86,
            82,
            82,
            82,
            88,
            87,
            83,
            85,
            85,
            84,
            # Phosphoramidites
            0,
            0,
            0,
            72,
            0,
            0,
            # Carbene Precursors
            0,
            0,
            0,
            0,
            0,
            0,
            # Other
            0,
            90,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            88,
            90,
        ],
        "Steric Bulk (Å³)": [
            # Monodentate Phosphines
            245,
            310,
            350,
            285,
            265,
            260,
            255,
            280,
            220,
            190,
            160,
            180,
            240,
            270,
            210,
            290,
            275,
            280,
            320,
            # Buchwald-type
            390,
            385,
            390,
            445,
            425,
            365,
            365,
            355,
            395,
            380,
            375,
            385,
            435,
            450,
            440,
            395,
            380,
            # BINAP Family
            520,
            525,
            515,
            510,
            520,
            530,
            535,
            545,
            525,
            515,
            530,
            # Other Bidentate
            360,
            370,
            380,
            425,
            485,
            445,
            420,
            410,
            415,
            390,
            420,
            415,
            425,
            430,
            445,
            440,
            435,
            430,
            # NHC Ligands
            420,
            400,
            420,
            400,
            425,
            440,
            445,
            415,
            380,
            390,
            385,
            380,
            430,
            410,
            # Nitrogen Ligands
            165,
            85,
            165,
            160,
            95,
            110,
            105,
            175,
            170,
            180,
            185,
            190,
            170,
            165,
            170,
            165,
            170,
            # Phosphoramidites
            280,
            285,
            290,
            310,
            295,
            285,
            # Carbene Precursors
            400,
            420,
            400,
            420,
            380,
            390,
            # Other
            120,
            120,
            110,
            85,
            95,
            105,
            95,
            100,
            80,
            110,
            105,
            95,
            90,
            125,
            130,
        ],
        "Donor Strength (pKa)": [
            # Monodentate Phosphines
            2.73,
            9.70,
            11.40,
            3.10,
            3.84,
            1.97,
            2.40,
            1.80,
            4.57,
            6.50,
            8.65,
            8.69,
            8.43,
            8.64,
            2.41,
            1.30,
            4.57,
            4.59,
            6.43,
            # Buchwald-type
            4.50,
            4.30,
            4.50,
            4.80,
            7.00,
            5.20,
            5.20,
            4.40,
            4.70,
            4.50,
            5.10,
            4.60,
            7.20,
            4.90,
            7.10,
            4.40,
            5.00,
            # BINAP Family
            2.75,
            2.77,
            2.76,
            2.74,
            2.75,
            2.77,
            2.78,
            2.79,
            2.76,
            2.73,
            2.77,
            # Other Bidentate
            2.76,
            2.77,
            2.78,
            2.80,
            2.85,
            2.83,
            2.84,
            2.79,
            2.81,
            2.78,
            2.82,
            2.81,
            2.85,
            2.83,
            2.84,
            2.80,
            2.78,
            2.79,
            # NHC Ligands
            8.50,
            8.30,
            8.50,
            8.30,
            8.45,
            8.70,
            8.75,
            8.40,
            8.15,
            8.20,
            8.10,
            8.05,
            8.45,
            8.35,
            # Nitrogen Ligands
            5.20,
            5.20,
            4.80,
            4.85,
            9.70,
            12.0,
            8.80,
            4.70,
            4.85,
            4.75,
            4.80,
            4.90,
            4.85,
            4.80,
            5.30,
            4.75,
            4.85,
            # Phosphoramidites
            3.50,
            3.55,
            3.45,
            3.40,
            3.50,
            3.45,
            # Carbene Precursors
            8.30,
            8.50,
            8.30,
            8.50,
            8.15,
            8.20,
            # Other
            8.80,
            8.80,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            5.20,
            5.25,
            8.70,
            8.85,
        ],
        "Price Category": [
            # 1: <$100/g, 2: $100-500/g, 3: $500-1000/g, 4: >$1000/g
            # Monodentate Phosphines
            1,
            2,
            3,
            1,
            1,
            2,
            1,
            2,
            1,
            1,
            1,
            1,
            1,
            2,
            2,
            3,
            2,
            2,
            2,
            # Buchwald-type
            3,
            3,
            3,
            4,
            4,
            3,
            3,
            2,
            4,
            3,
            3,
            3,
            4,
            4,
            4,
            3,
            3,
            # BINAP Family
            3,
            3,
            3,
            3,
            3,
            4,
            4,
            4,
            4,
            4,
            4,
            # Other Bidentate
            1,
            1,
            1,
            2,
            3,
            3,
            4,
            3,
            4,
            4,
            4,
            4,
            4,
            4,
            3,
            4,
            4,
            4,
            # NHC Ligands
            3,
            2,
            3,
            2,
            3,
            4,
            4,
            3,
            2,
            2,
            3,
            3,
            3,
            3,
            # Nitrogen Ligands
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            2,
            2,
            2,
            2,
            3,
            2,
            2,
            1,
            2,
            2,
            # Phosphoramidites
            3,
            3,
            3,
            4,
            3,
            3,
            # Carbene Precursors
            2,
            3,
            2,
            3,
            2,
            2,
            # Other
            1,
            1,
            1,
            1,
            2,
            2,
            2,
            2,
            1,
            1,
            1,
            1,
            1,
            1,
            2,
        ],
        "Coordination Mode": [
            # 1: Monodentate, 2: Bidentate, 3: Tridentate, 4: Other
            # Monodentate Phosphines (all 1)
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            # Buchwald-type (all 1)
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            # BINAP Family (all 2)
            2,
            2,
            2,
            2,
            2,
            2,
            2,
            2,
            2,
            2,
            2,
            # Other Bidentate (all 2)
            2,
            2,
            2,
            2,
            2,
            2,
            2,
            2,
            2,
            2,
            2,
            2,
            2,
            2,
            2,
            2,
            2,
            2,
            # NHC Ligands (all 1)
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            # Nitrogen Ligands
            2,
            1,
            2,
            2,
            1,
            1,
            1,
            3,
            2,
            2,
            2,
            3,
            2,
            2,
            2,
            2,
            2,
            # Phosphoramidites (all 1)
            1,
            1,
            1,
            2,
            1,
            1,
            # Carbene Precursors (all 1)
            1,
            1,
            1,
            1,
            1,
            1,
            # Other
            2,
            2,
            2,
            1,
            1,
            1,
            1,
            1,
            1,
            2,
            2,
            1,
            1,
            2,
            2,
        ],
        "Reaction_Compatibility": [
            # Compatibility scores for different reaction types (comma-separated string)
            # Format: "Cross-Coupling,Hydrogenation,Metathesis,C-H_Activation,Carbonylation"
            # Values 0-1 (0: incompatible, 1: highly compatible)
            # Monodentate Phosphines
            "0.8,0.9,0.3,0.6,0.8",
            "0.7,0.9,0.2,0.5,0.7",
            "0.8,0.8,0.3,0.7,0.8",
            "0.7,0.8,0.3,0.5,0.7",
            "0.7,0.8,0.3,0.5,0.7",
            "0.6,0.7,0.3,0.4,0.6",
            "0.6,0.7,0.3,0.4,0.6",
            "0.5,0.6,0.3,0.4,0.5",
            "0.7,0.8,0.3,0.5,0.7",
            "0.6,0.8,0.3,0.4,0.6",
            "0.5,0.7,0.3,0.3,0.5",
            "0.5,0.7,0.3,0.3,0.5",
            "0.6,0.8,0.3,0.4,0.6",
            "0.6,0.8,0.3,0.4,0.6",
            "0.5,0.7,0.3,0.3,0.5",
            "0.4,0.5,0.3,0.3,0.4",
            "0.6,0.7,0.3,0.4,0.6",
            "0.6,0.7,0.3,0.4,0.6",
            "0.5,0.6,0.3,0.3,0.5",
            # Buchwald-type Ligands (optimized for cross-coupling)
            "0.9,0.7,0.2,0.6,0.8",
            "0.9,0.7,0.2,0.6,0.8",
            "0.9,0.7,0.2,0.6,0.8",
            "0.9,0.6,0.2,0.7,0.8",
            "0.9,0.6,0.2,0.7,0.8",
            "0.8,0.7,0.2,0.6,0.8",
            "0.8,0.7,0.2,0.6,0.8",
            "0.8,0.7,0.2,0.5,0.7",
            "0.9,0.6,0.2,0.6,0.8",
            "0.8,0.7,0.2,0.6,0.8",
            "0.8,0.7,0.2,0.6,0.8",
            "0.8,0.7,0.2,0.6,0.8",
            "0.9,0.6,0.2,0.7,0.8",
            "0.9,0.6,0.2,0.7,0.8",
            "0.9,0.6,0.2,0.7,0.8",
            "0.8,0.7,0.2,0.6,0.8",
            "0.8,0.7,0.2,0.6,0.8",
            # BINAP Family (excellent for asymmetric hydrogenation)
            "0.8,0.95,0.3,0.6,0.8",
            "0.8,0.95,0.3,0.6,0.8",
            "0.8,0.95,0.3,0.6,0.8",
            "0.8,0.9,0.3,0.6,0.8",
            "0.8,0.9,0.3,0.6,0.8",
            "0.8,0.95,0.3,0.6,0.8",
            "0.8,0.95,0.3,0.6,0.8",
            "0.8,0.95,0.3,0.6,0.8",
            "0.8,0.9,0.3,0.6,0.8",
            "0.8,0.9,0.3,0.6,0.8",
            "0.8,0.9,0.3,0.6,0.8",
            # Other Bidentate Phosphines
            "0.7,0.8,0.3,0.5,0.7",
            "0.7,0.8,0.3,0.5,0.7",
            "0.7,0.8,0.3,0.5,0.7",
            "0.8,0.8,0.3,0.6,0.8",
            "0.8,0.8,0.4,0.6,0.8",
            "0.8,0.8,0.4,0.6,0.8",
            "0.8,0.9,0.3,0.6,0.8",
            "0.7,0.8,0.3,0.5,0.7",
            "0.8,0.9,0.3,0.6,0.8",
            "0.7,0.8,0.3,0.5,0.7",
            "0.8,0.8,0.3,0.6,0.8",
            "0.8,0.9,0.3,0.6,0.8",
            "0.8,0.8,0.4,0.6,0.8",
            "0.8,0.8,0.4,0.6,0.8",
            "0.8,0.8,0.4,0.6,0.8",
            "0.8,0.8,0.4,0.6,0.8",
            "0.8,0.8,0.4,0.6,0.8",
            "0.8,0.8,0.4,0.6,0.8",
            # NHC Ligands (excellent for metathesis and some cross-coupling)
            "0.7,0.4,0.9,0.8,0.6",
            "0.7,0.4,0.9,0.8,0.6",
            "0.7,0.4,0.9,0.8,0.6",
            "0.7,0.4,0.9,0.8,0.6",
            "0.7,0.4,0.9,0.8,0.6",
            "0.8,0.4,0.9,0.8,0.6",
            "0.8,0.4,0.9,0.8,0.6",
            "0.8,0.4,0.9,0.8,0.6",
            "0.6,0.4,0.8,0.7,0.5",
            "0.7,0.4,0.9,0.8,0.6",
            "0.6,0.4,0.8,0.7,0.5",
            "0.6,0.4,0.8,0.7,0.5",
            "0.7,0.4,0.9,0.8,0.6",
            "0.8,0.4,0.9,0.8,0.6",
            # Nitrogen Ligands
            "0.3,0.6,0.2,0.5,0.4",
            "0.5,0.7,0.2,0.6,0.5",
            "0.4,0.8,0.2,0.5,0.4",
            "0.4,0.8,0.2,0.5,0.4",
            "0.6,0.7,0.2,0.7,0.6",
            "0.2,0.3,0.2,0.4,0.3",
            "0.2,0.3,0.2,0.4,0.3",
            "0.4,0.8,0.2,0.5,0.4",
            "0.4,0.8,0.2,0.5,0.4",
            "0.4,0.8,0.2,0.5,0.4",
            "0.4,0.8,0.2,0.5,0.4",
            "0.5,0.8,0.2,0.6,0.5",
            "0.5,0.8,0.2,0.6,0.5",
            "0.4,0.8,0.2,0.5,0.4",
            "0.3,0.6,0.2,0.5,0.4",
            "0.4,0.7,0.2,0.5,0.4",
            "0.4,0.7,0.2,0.5,0.4",
            # Phosphoramidites (good for asymmetric reactions)
            "0.7,0.9,0.3,0.5,0.7",
            "0.7,0.9,0.3,0.5,0.7",
            "0.7,0.9,0.3,0.5,0.7",
            "0.8,0.9,0.3,0.6,0.8",
            "0.7,0.9,0.3,0.5,0.7",
            "0.7,0.9,0.3,0.5,0.7",
            # Carbene Precursors
            "0.7,0.4,0.9,0.8,0.6",
            "0.7,0.4,0.9,0.8,0.6",
            "0.7,0.4,0.9,0.8,0.6",
            "0.7,0.4,0.9,0.8,0.6",
            "0.6,0.4,0.8,0.7,0.5",
            "0.7,0.4,0.9,0.8,0.6",
            # Other/Specialty
            "0.3,0.3,0.5,0.4,0.3",
            "0.5,0.5,0.3,0.4,0.5",
            "0.6,0.6,0.7,0.5,0.6",
            "0.4,0.4,0.6,0.5,0.4",
            "0.4,0.4,0.6,0.5,0.4",
            "0.4,0.4,0.6,0.5,0.4",
            "0.4,0.4,0.6,0.5,0.4",
            "0.4,0.4,0.6,0.5,0.4",
            "0.5,0.5,0.7,0.6,0.5",
            "0.6,0.6,0.7,0.5,0.6",
            "0.6,0.6,0.7,0.5,0.6",
            "0.4,0.6,0.3,0.4,0.5",
            "0.4,0.6,0.3,0.4,0.4",
            "0.5,0.5,0.6,0.5,0.5",
            "0.5,0.5,0.6,0.5,0.5",
        ],
        "Typical_Applications": [
            # Typical applications for each ligand
            # Monodentate Phosphines
            "General cross-coupling, hydrogenation",
            "Hydrogenation, cross-coupling",
            "Cross-coupling, C-H activation",
            "Cross-coupling",
            "Cross-coupling",
            "Electron-poor substrates",
            "Electron-poor substrates",
            "Specialized cross-coupling",
            "General catalysis",
            "Hydrogenation",
            "Small substrate hydrogenation",
            "Hydrogenation",
            "Cross-coupling",
            "Cross-coupling",
            "Specialized applications",
            "Electron-poor substrates",
            "General catalysis",
            "General catalysis",
            "Bulky substrate catalysis",
            # Buchwald-type Ligands
            "Buchwald-Hartwig amination",
            "Suzuki-Miyaura coupling",
            "Buchwald-Hartwig amination",
            "Challenging cross-coupling",
            "Sterically hindered substrates",
            "General cross-coupling",
            "General cross-coupling",
            "Cross-coupling",
            "Specialty cross-coupling",
            "Cross-coupling",
            "Cross-coupling",
            "Cross-coupling",
            "Challenging substrates",
            "Bulky substrates",
            "Bulky substrates",
            "Specialty coupling",
            "Specialty coupling",
            # BINAP Family
            "Asymmetric hydrogenation",
            "Asymmetric hydrogenation",
            "Asymmetric hydrogenation",
            "Asymmetric reactions",
            "Asymmetric reactions",
            "Asymmetric hydrogenation",
            "Asymmetric hydrogenation",
            "Asymmetric hydrogenation",
            "Asymmetric reactions",
            "Asymmetric reactions",
            "Asymmetric reactions",
            # Other Bidentate Phosphines
            "Cross-coupling",
            "Cross-coupling",
            "Cross-coupling",
            "Cross-coupling, carbonylation",
            "Cross-coupling, carbonylation",
            "Cross-coupling, carbonylation",
            "Asymmetric catalysis",
            "Cross-coupling",
            "Asymmetric hydrogenation",
            "Cross-coupling",
            "Cross-coupling",
            "Asymmetric catalysis",
            "Cross-coupling, carbonylation",
            "Cross-coupling, carbonylation",
            "Cross-coupling, carbonylation",
            "Cross-coupling, carbonylation",
            "Cross-coupling, carbonylation",
            "Specialty coupling",
            # NHC Ligands
            "Metathesis, cross-coupling",
            "Metathesis, cross-coupling",
            "Metathesis, cross-coupling",
            "Metathesis, cross-coupling",
            "Cross-coupling",
            "Metathesis",
            "Metathesis",
            "Metathesis",
            "Cross-coupling",
            "Cross-coupling",
            "Specialty applications",
            "Specialty applications",
            "Metathesis",
            "Metathesis",
            # Nitrogen Ligands
            "Coordination",
            "Coordination, C-H activation",
            "Coordination",
            "Coordination",
            "Nucleophilic catalysis",
            "Base catalysis",
            "Base catalysis",
            "Coordination",
            "Coordination",
            "Coordination",
            "Coordination",
            "Asymmetric catalysis",
            "Asymmetric catalysis",
            "Coordination",
            "Coordination",
            "Coordination",
            "Coordination",
            # Phosphoramidites
            "Asymmetric hydrogenation",
            "Asymmetric hydrogenation",
            "Asymmetric catalysis",
            "Asymmetric catalysis",
            "Asymmetric catalysis",
            "Asymmetric catalysis",
            # Carbene Precursors
            "Metathesis precursor",
            "Metathesis precursor",
            "Metathesis precursor",
            "Metathesis precursor",
            "Cross-coupling precursor",
            "Cross-coupling precursor",
            # Other/Specialty
            "Stabilizing ligand",
            "Stabilizing ligand",
            "Diene ligand",
            "Cyclic ligand",
            "Cyclic ligand",
            "Cyclic ligand",
            "Cyclic ligand",
            "Cyclic ligand",
            "π-Allyl ligand",
            "Diene ligand",
            "Diene ligand",
            "Nitrile ligand",
            "Pyridine derivative",
            "β-Diketonate",
            "β-Diketonate",
        ],
    }
)


def create_ligand_dataframe():
    """Create the ligand DataFrame from the ligand_data dictionary"""
    return pd.DataFrame(ligand_data)


def create_feature_matrix():
    """Create a normalized feature matrix for machine learning operations"""
    df = create_ligand_dataframe()

    # Select numerical features (using actual column names from DataFrame)
    feature_columns = [
        "Cone Angle (°)",
        "Electronic Parameter (cm⁻¹)",
        "Bite Angle (°)",
        "Steric Bulk (Å³)",
        "Donor Strength (pKa)",
        "Price Category",
        "Coordination Mode",
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


# Reaction-specific property weights
REACTION_WEIGHTS = {
    "Cross-Coupling": {
        "Cone Angle (°)": 0.25,
        "Electronic Parameter (cm⁻¹)": 0.20,
        "Bite Angle (°)": 0.15,
        "Steric Bulk (Å³)": 0.15,
        "Donor Strength (pKa)": 0.10,
        "Price Category": 0.05,
        "Coordination Mode": 0.10,
    },
    "Hydrogenation": {
        "Cone Angle (°)": 0.15,
        "Electronic Parameter (cm⁻¹)": 0.30,
        "Bite Angle (°)": 0.15,
        "Steric Bulk (Å³)": 0.10,
        "Donor Strength (pKa)": 0.15,
        "Price Category": 0.05,
        "Coordination Mode": 0.10,
    },
    "Metathesis": {
        "Cone Angle (°)": 0.20,
        "Electronic Parameter (cm⁻¹)": 0.25,
        "Bite Angle (°)": 0.10,
        "Steric Bulk (Å³)": 0.15,
        "Donor Strength (pKa)": 0.05,
        "Price Category": 0.10,
        "Coordination Mode": 0.15,
    },
    "C-H_Activation": {
        "Cone Angle (°)": 0.20,
        "Electronic Parameter (cm⁻¹)": 0.25,
        "Bite Angle (°)": 0.10,
        "Steric Bulk (Å³)": 0.10,
        "Donor Strength (pKa)": 0.15,
        "Price Category": 0.05,
        "Coordination Mode": 0.15,
    },
    "Carbonylation": {
        "Cone Angle (°)": 0.20,
        "Electronic Parameter (cm⁻¹)": 0.20,
        "Bite Angle (°)": 0.15,
        "Steric Bulk (Å³)": 0.15,
        "Donor Strength (pKa)": 0.15,
        "Price Category": 0.05,
        "Coordination Mode": 0.10,
    },
}


def parse_reaction_compatibility(compatibility_str, reaction_type):
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


def calculate_weighted_similarity(ligand1_features, ligand2_features, weights):
    """Calculate weighted similarity between two ligands based on reaction-specific weights"""
    similarity = 0.0

    feature_names = [
        "Cone Angle (°)",
        "Electronic Parameter (cm⁻¹)",
        "Bite Angle (°)",
        "Steric Bulk (Å³)",
        "Donor Strength (pKa)",
        "Price Category",
        "Coordination Mode",
    ]

    for i, feature in enumerate(feature_names):
        if feature in weights:
            # Normalized difference (0 = identical, 1 = maximum difference)
            max_val = max(abs(ligand1_features[i]), abs(ligand2_features[i]), 1)
            diff = abs(ligand1_features[i] - ligand2_features[i]) / max_val
            similarity += weights[feature] * (1 - diff)

    return similarity


def recommend_ligands_for_reaction(
    target_ligand=None, reaction_type="Cross-Coupling", top_n=5, min_compatibility=0.3
):
    """
    Recommend ligands for a specific reaction type.

    Parameters:
    - target_ligand: Name of reference ligand (optional)
    - reaction_type: Type of reaction ('Cross-Coupling', 'Hydrogenation', 'Metathesis',
                    'C-H_Activation', 'Carbonylation')
    - top_n: Number of recommendations to return
    - min_compatibility: Minimum compatibility score for the reaction type

    Returns:
    - List of dictionaries with ligand recommendations
    """
    df = create_ligand_dataframe()

    # Filter ligands by reaction compatibility
    compatible_ligands = []
    for idx, row in df.iterrows():
        compatibility = parse_reaction_compatibility(
            row["Reaction_Compatibility"], reaction_type
        )
        if compatibility >= min_compatibility:
            compatible_ligands.append(
                {
                    "index": idx,
                    "name": row["Ligand"],
                    "compatibility": compatibility,
                    "applications": row["Typical_Applications"],
                }
            )

    # Sort by compatibility score
    compatible_ligands.sort(key=lambda x: x["compatibility"], reverse=True)

    # If target ligand is specified, calculate similarity scores
    if target_ligand:
        target_idx = None
        for idx, name in enumerate(df["Ligand"]):
            if name.lower() == target_ligand.lower():
                target_idx = idx
                break

        if target_idx is not None:
            X = create_feature_matrix()
            target_features = X[target_idx]
            weights = REACTION_WEIGHTS.get(
                reaction_type, REACTION_WEIGHTS["Cross-Coupling"]
            )

            # Calculate similarity scores
            for ligand in compatible_ligands:
                if ligand["index"] != target_idx:  # Don't include the target itself
                    ligand_features = X[ligand["index"]]
                    similarity = calculate_weighted_similarity(
                        target_features, ligand_features, weights
                    )
                    ligand["similarity"] = similarity
                else:
                    ligand["similarity"] = 0  # Remove target from recommendations

            # Sort by similarity score (keeping only compatible ligands)
            compatible_ligands = [
                l for l in compatible_ligands if l["index"] != target_idx
            ]
            compatible_ligands.sort(
                key=lambda x: (x["compatibility"] * 0.6 + x["similarity"] * 0.4),
                reverse=True,
            )

    # Return top recommendations
    recommendations = []
    for i, ligand in enumerate(compatible_ligands[:top_n]):
        rec = {
            "rank": i + 1,
            "ligand": ligand["name"],
            "compatibility_score": round(ligand["compatibility"], 3),
            "applications": ligand["applications"],
            "reaction_suitability": reaction_type,
        }

        if "similarity" in ligand:
            rec["similarity_score"] = round(ligand["similarity"], 3)
            rec["combined_score"] = round(
                ligand["compatibility"] * 0.6 + ligand["similarity"] * 0.4, 3
            )

        recommendations.append(rec)

    return recommendations


def get_reaction_specific_ligands(reaction_type, property_preferences=None):
    """
    Get ligands optimized for a specific reaction type with optional property preferences.

    Parameters:
    - reaction_type: Type of reaction
    - property_preferences: Dict with preferred ranges for properties

    Returns:
    - Filtered ligand recommendations
    """
    recommendations = recommend_ligands_for_reaction(
        reaction_type=reaction_type, top_n=10, min_compatibility=0.4
    )

    if property_preferences:
        df = create_ligand_dataframe()
        filtered_recs = []

        for rec in recommendations:
            ligand_name = rec["ligand"]
            ligand_idx = df[df["Ligand"] == ligand_name].index[0]

            # Check property preferences
            meets_criteria = True
            if "cone_angle_max" in property_preferences:
                if (
                    df.loc[ligand_idx, "Cone Angle (°)"]
                    > property_preferences["cone_angle_max"]
                ):
                    meets_criteria = False

            if "price_category_max" in property_preferences:
                if (
                    df.loc[ligand_idx, "Price Category"]
                    > property_preferences["price_category_max"]
                ):
                    meets_criteria = False

            if "coordination_mode" in property_preferences:
                if (
                    df.loc[ligand_idx, "Coordination Mode"]
                    != property_preferences["coordination_mode"]
                ):
                    meets_criteria = False

            if meets_criteria:
                filtered_recs.append(rec)

        return filtered_recs[:5]  # Return top 5 that meet criteria

    return recommendations


# Define property weights
property_weights = {
    "Cone Angle (°)": 0.20,
    "Electronic Parameter (cm⁻¹)": 0.20,
    "Bite Angle (°)": 0.18,
    "Steric Bulk (Å³)": 0.15,
    "Donor Strength (pKa)": 0.15,
    "Price Category": 0.07,
    "Coordination Mode": 0.05,
}

# Export ligand database to Excel
excel_path = "ligand_database.xlsx"
with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
    # Main data sheet
    ligand_data.to_excel(writer, sheet_name="Ligand Properties", index=False)

    # Additional information sheet
    info_data = pd.DataFrame(
        {
            "Property": list(property_weights.keys()),
            "Weight": list(property_weights.values()),
            "Description": [
                "Tolman cone angle - measure of steric bulk",
                "Tolman Electronic Parameter - measure of electron density",
                "Natural bite angle for bidentate ligands",
                "Molecular volume/spatial requirement",
                "Basicity/electron-donating ability",
                "Relative cost category",
                "Binding mode to metal center",
            ],
        }
    )
    info_data.to_excel(writer, sheet_name="Property Information", index=False)

print(f"Ligand database has been exported to {excel_path}")


# Function to recommend similar ligands
def recommend_ligands(selected_ligand, num_recommendations=4):
    if selected_ligand not in ligand_data["Ligand"].values:
        print(f"Ligand {selected_ligand} not found in the database.")
        return []

    # Normalize the properties to 0-1 scale
    properties_normalized = ligand_data.copy()
    for column in property_weights.keys():
        min_val = properties_normalized[column].min()
        max_val = properties_normalized[column].max()
        properties_normalized[column] = (properties_normalized[column] - min_val) / (
            max_val - min_val
        )

    # Apply weights to properties
    properties_weighted = properties_normalized.copy()
    for column, weight in property_weights.items():
        properties_weighted[column] = properties_weighted[column] * weight

    # Extract weighted properties for similarity calculation
    ligand_properties = properties_weighted[list(property_weights.keys())].values
    ligand_names = ligand_data["Ligand"].values

    # Find the index of the selected ligand
    selected_idx = np.where(ligand_names == selected_ligand)[0][0]
    selected_vector = ligand_properties[selected_idx].reshape(1, -1)

    # Compute weighted Euclidean distances
    distances = cdist(selected_vector, ligand_properties, metric="euclidean")[0]

    # Sort and get the closest ligands (excluding itself)
    closest_indices = np.argsort(distances)[1 : num_recommendations + 1]
    recommended_ligands = ligand_names[closest_indices]

    # Add similarity scores
    similarity_scores = 1 / (1 + distances[closest_indices])
    recommendations = list(zip(recommended_ligands, similarity_scores))

    return recommendations


def analyze_ligand_clusters():
    # Normalize all properties for clustering
    properties_normalized = ligand_data.copy()
    for column in property_weights.keys():
        min_val = properties_normalized[column].min()
        max_val = properties_normalized[column].max()
        properties_normalized[column] = (properties_normalized[column] - min_val) / (
            max_val - min_val
        )

    # Apply weights to properties
    properties_weighted = properties_normalized.copy()
    for column, weight in property_weights.items():
        properties_weighted[column] = properties_weighted[column] * weight

    # Prepare data for clustering
    X = properties_weighted[list(property_weights.keys())].values

    # Perform hierarchical clustering
    Z = linkage(X, method="ward")

    # Create dendrogram
    plt.figure(figsize=(15, 10))
    plt.title("Hierarchical Clustering of Ligands")
    dendrogram(
        Z, labels=ligand_data["Ligand"].values, leaf_rotation=90, leaf_font_size=8
    )
    plt.xlabel("Ligands")
    plt.ylabel("Distance")

    # Adjust layout to prevent label cutoff
    plt.tight_layout()

    # Save the plot
    plt.savefig("ligand_clusters.png", dpi=300, bbox_inches="tight")
    plt.close()

    print(
        "Clustering analysis completed. The dendrogram has been saved as 'ligand_clusters.png'"
    )


def create_ligand_network(threshold=0.7):
    # Normalize all properties
    properties_normalized = ligand_data.copy()
    scaler = MinMaxScaler()
    properties = scaler.fit_transform(ligand_data[list(property_weights.keys())])

    # Calculate similarity matrix
    similarity_matrix = 1 / (1 + cdist(properties, properties, metric="euclidean"))

    # Create network graph
    G = nx.Graph()

    # Add nodes (ligands)
    for ligand in ligand_data["Ligand"]:
        G.add_node(ligand)

    # Add edges based on similarity threshold
    for i in range(len(ligand_data)):
        for j in range(i + 1, len(ligand_data)):
            if similarity_matrix[i, j] > threshold:
                G.add_edge(
                    ligand_data["Ligand"].iloc[i],
                    ligand_data["Ligand"].iloc[j],
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

    plt.title("Ligand Similarity Network\n(Edges show similarities above threshold)")
    plt.axis("off")

    # Save the network plot
    plt.savefig("ligand_network.png", dpi=300, bbox_inches="tight")
    plt.close()

    print(
        "Network analysis completed. The graph has been saved as 'ligand_network.png'"
    )
    return G


# Example Usage
if __name__ == "__main__":
    # Test ligand recommendation
    selected_ligand = "DPPF"
    recommended = recommend_ligands(selected_ligand, num_recommendations=4)
    print(f"\nRecommendations for {selected_ligand}:")
    for ligand, score in recommended:
        print(f"{ligand}: {score:.2f} similarity score")

    # Create visualizations
    analyze_ligand_clusters()
    create_ligand_network(threshold=0.7)
