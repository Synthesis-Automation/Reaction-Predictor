import os
import json
import numpy as np
import pandas as pd

# Optional dependencies (graceful degradation if missing)
try:  # scikit-learn for scaling
    from sklearn.preprocessing import MinMaxScaler  # type: ignore
except Exception:  # pragma: no cover
    MinMaxScaler = None  # type: ignore

try:  # SciPy for distances
    from scipy.spatial.distance import cdist  # type: ignore
except Exception:  # pragma: no cover
    cdist = None  # type: ignore

try:  # SciPy for clustering
    from scipy.cluster.hierarchy import dendrogram, linkage  # type: ignore
except Exception:  # pragma: no cover
    dendrogram = None  # type: ignore
    linkage = None  # type: ignore

try:  # Matplotlib for plots
    import matplotlib.pyplot as plt  # type: ignore
except Exception:  # pragma: no cover
    plt = None  # type: ignore

try:  # NetworkX for graphs
    import networkx as nx  # type: ignore
except Exception:  # pragma: no cover
    nx = None  # type: ignore

try:  # Excel writer backend
    import openpyxl  # noqa: F401
    OPENPYXL_AVAILABLE = True
except Exception:  # pragma: no cover
    OPENPYXL_AVAILABLE = False

EXPECTED_SOLVENT_COLUMNS = [
    "Solvent",
    "CAS Number",
    "Abbreviation",
    "Dielectric Constant",
    "Polarity Index",
    "Boiling Point (°C)",
    "Density (g/mL)",
    "Dipole Moment (D)",
    "Donor Number (DN)",
    "Hydrogen Bond Donor",
    "Reaction_Compatibility",
    "Typical_Applications",
]


def _default_data_paths():
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, 'data'))
    return (
        os.path.join(base_dir, 'solvents.json'),
        os.path.join(base_dir, 'solvents'),
    )


def create_solvent_dataframe(json_path: str | None = None, json_dir: str | None = None, prefer_builtin: bool = False):
    """
    Create the solvent DataFrame from JSON.

    Supports either:
      - data/solvents.json (array or {"solvents": [...]})
      - data/solvents/ (directory of *.json with one solvent per file)
    """
    default_file, default_dir = _default_data_paths()
    json_path = json_path or default_file
    json_dir = json_dir or default_dir

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
            "Boiling Point (°C)": entry.get('boiling_point_c') or entry.get('Boiling Point (°C)') or entry.get('Boiling Point (C)'),
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
    except Exception as e:  # pragma: no cover
        print(f"Warning: failed to load solvents JSON: {e}")
        records = []

    if records:
        df = pd.DataFrame.from_records(records)
        df = df[df['Solvent'].notna() & (df['Solvent'].astype(str).str.len() > 0)]
        return df.reset_index(drop=True)

    # Fallback: return empty DataFrame with expected columns
    return pd.DataFrame(columns=EXPECTED_SOLVENT_COLUMNS)


# Back-compat: module-level DataFrame used by some functions
solvent_data = create_solvent_dataframe()


def create_solvent_feature_matrix():
    """Create a normalized feature matrix for machine learning operations"""
    df = create_solvent_dataframe()

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
            features.append(row[col] if col in df.columns else 0)
        X.append(features)

    X = np.array(X)
    if MinMaxScaler is None:
        # No scaling available; return raw features
        return X
    scaler = MinMaxScaler()
    return scaler.fit_transform(X)


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
    reaction_types = [
        "Cross-Coupling",
        "Hydrogenation",
        "Metathesis",
        "C-H_Activation",
        "Carbonylation",
    ]

    # Map specialized types to base scoring buckets
    if isinstance(reaction_type, str) and reaction_type.lower() == 'ullmann':
        reaction_type = 'Cross-Coupling'

    if reaction_type not in reaction_types:
        return 0.5

    try:
        scores = [float(x) for x in str(compatibility_str).split(",")]
        reaction_idx = reaction_types.index(reaction_type)
        return scores[reaction_idx] if reaction_idx < len(scores) else 0.5
    except Exception:
        return 0.5


def calculate_solvent_weighted_similarity(solvent1_features, solvent2_features, weights):
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
            max_val = max(abs(solvent1_features[i]), abs(solvent2_features[i]), 1)
            diff = abs(solvent1_features[i] - solvent2_features[i]) / max_val
            similarity += weights[feature] * (1 - diff)

    return similarity


def recommend_solvents_for_reaction(target_solvent=None, reaction_type="Cross-Coupling", top_n=5, min_compatibility=0.3):
    df = create_solvent_dataframe()

    compatible_solvents = []
    for idx, row in df.iterrows():
        compatibility = parse_solvent_reaction_compatibility(row.get("Reaction_Compatibility", ""), reaction_type)
        if compatibility >= min_compatibility:
            compatible_solvents.append({
                "index": idx,
                "name": row.get("Solvent", row.get("name", "")),
                "compatibility": compatibility,
                "applications": row.get("Typical_Applications", ""),
                "abbreviation": row.get("Abbreviation", ""),
            })

    compatible_solvents.sort(key=lambda x: x["compatibility"], reverse=True)

    if target_solvent:
        target_idx = None
        for idx, name in enumerate(df["Solvent"]):
            if isinstance(name, str) and name.lower() == str(target_solvent).lower():
                target_idx = idx
                break

        if target_idx is not None and len(compatible_solvents) > 1 and MinMaxScaler is not None and cdist is not None:
            X = create_solvent_feature_matrix()
            if X.size and target_idx < len(X):
                target_features = X[target_idx]
                weights = SOLVENT_REACTION_WEIGHTS.get(reaction_type, SOLVENT_REACTION_WEIGHTS["Cross-Coupling"])

                for solvent in compatible_solvents:
                    if solvent["index"] != target_idx:
                        solvent_features = X[solvent["index"]]
                        similarity = calculate_solvent_weighted_similarity(target_features, solvent_features, weights)
                        solvent["similarity"] = similarity
                    else:
                        solvent["similarity"] = 0

                compatible_solvents = [s for s in compatible_solvents if s["index"] != target_idx]
                compatible_solvents.sort(key=lambda x: (x["compatibility"] * 0.6 + x.get("similarity", 0) * 0.4), reverse=True)

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
            rec["combined_score"] = round(solvent["compatibility"] * 0.6 + solvent.get("similarity", 0) * 0.4, 3)
        recommendations.append(rec)

    return recommendations


def get_reaction_specific_solvents(reaction_type, property_preferences=None):
    recommendations = recommend_solvents_for_reaction(reaction_type=reaction_type, top_n=10, min_compatibility=0.4)

    if property_preferences:
        df = create_solvent_dataframe()
        filtered_recs = []

        for rec in recommendations:
            solvent_name = rec["solvent"]
            matches = df[df["Solvent"] == solvent_name]
            if matches.empty:
                continue
            solvent_idx = matches.index[0]

            meets_criteria = True
            if "bp_max" in property_preferences:
                if df.loc[solvent_idx, "Boiling Point (°C)"] > property_preferences["bp_max"]:
                    meets_criteria = False
            if "bp_min" in property_preferences:
                if df.loc[solvent_idx, "Boiling Point (°C)"] < property_preferences["bp_min"]:
                    meets_criteria = False
            if "polarity_max" in property_preferences:
                if df.loc[solvent_idx, "Polarity Index"] > property_preferences["polarity_max"]:
                    meets_criteria = False
            if "polarity_min" in property_preferences:
                if df.loc[solvent_idx, "Polarity Index"] < property_preferences["polarity_min"]:
                    meets_criteria = False
            if "protic" in property_preferences:
                is_protic = float(df.loc[solvent_idx, "Hydrogen Bond Donor"] or 0) > 0.5
                if property_preferences["protic"] != is_protic:
                    meets_criteria = False

            if meets_criteria:
                filtered_recs.append(rec)

        return filtered_recs[:5]

    return recommendations


property_weights = {
    "CAS Number": 0.00,
    "Abbreviation": 0.00,
    "Dielectric Constant": 0.08,
    "Polarity Index": 0.27,
    "Boiling Point (°C)": 0.02,
    "Density (g/mL)": 0.02,
    "Dipole Moment (D)": 0.11,
    "Donor Number (DN)": 0.20,
    "Hydrogen Bond Donor": 0.30,
}


## Removed automatic Excel export (solvent_database.xlsx) on import as it's not relevant for GUI runtime.
## If needed, provide an explicit utility function elsewhere to generate the workbook on demand.


# Legacy helpers using optional libs; guarded for availability

def recommend_solvents(selected_solvent, num_recommendations=3):
    if 'Solvent' not in solvent_data.columns or selected_solvent not in solvent_data['Solvent'].values:
        print(f"Solvent {selected_solvent} not found in the database.")
        return []

    if cdist is None:
        return []

    numeric_columns = [col for col in property_weights.keys() if col not in ["CAS Number", "Abbreviation"]]
    properties_normalized = solvent_data.copy()

    for column in numeric_columns:
        min_val = properties_normalized[column].min()
        max_val = properties_normalized[column].max()
        if max_val == min_val:
            properties_normalized[column] = 0
        else:
            properties_normalized[column] = (properties_normalized[column] - min_val) / (max_val - min_val)

    properties_weighted = properties_normalized.copy()
    for column in numeric_columns:
        properties_weighted[column] = (properties_weighted[column] * property_weights[column])

    solvent_properties = properties_weighted[numeric_columns].values
    solvent_names = solvent_data["Solvent"].values

    selected_idx = np.where(solvent_names == selected_solvent)[0][0]
    selected_vector = solvent_properties[selected_idx].reshape(1, -1)

    distances = cdist(selected_vector, solvent_properties, metric="euclidean")[0]
    closest_indices = np.argsort(distances)[1 : num_recommendations + 1]
    recommended_solvents = solvent_names[closest_indices]
    similarity_scores = 1 / (1 + distances[closest_indices])
    return list(zip(recommended_solvents, similarity_scores))


def analyze_solvent_clusters():
    if linkage is None or plt is None:
        return

    numeric_columns = [col for col in property_weights.keys() if col not in ["CAS Number", "Abbreviation"]]
    properties_normalized = solvent_data.copy()
    for column in numeric_columns:
        min_val = properties_normalized[column].min()
        max_val = properties_normalized[column].max()
        if max_val == min_val:
            properties_normalized[column] = 0
        else:
            properties_normalized[column] = (properties_normalized[column] - min_val) / (max_val - min_val)

    properties_weighted = properties_normalized.copy()
    for column in numeric_columns:
        properties_weighted[column] = (properties_weighted[column] * property_weights[column])

    X = properties_weighted[numeric_columns].values
    Z = linkage(X, method="ward")
    plt.figure(figsize=(15, 10))
    plt.title("Hierarchical Clustering of Solvents")
    dendrogram(Z, labels=solvent_data["Solvent"].values, leaf_rotation=90, leaf_font_size=8)
    plt.xlabel("Solvents")
    plt.ylabel("Distance")
    plt.tight_layout()
    plt.savefig("solvent_clusters.png", dpi=300, bbox_inches="tight")
    plt.close()


def create_solvent_network(threshold=0.7):
    if cdist is None or nx is None or MinMaxScaler is None:
        return None

    numeric_columns = [col for col in property_weights.keys() if col not in ["CAS Number", "Abbreviation"]]
    properties_normalized = solvent_data[numeric_columns].copy()
    scaler = MinMaxScaler()
    properties = scaler.fit_transform(properties_normalized)

    similarity_matrix = 1 / (1 + cdist(properties, properties, metric="euclidean"))
    G = nx.Graph()

    for solvent in solvent_data["Solvent"]:
        G.add_node(solvent)

    for i in range(len(solvent_data)):
        for j in range(i + 1, len(solvent_data)):
            if similarity_matrix[i, j] > threshold:
                G.add_edge(solvent_data["Solvent"].iloc[i], solvent_data["Solvent"].iloc[j], weight=similarity_matrix[i, j])

    if plt is None:
        return G

    plt.figure(figsize=(20, 20))
    pos = nx.spring_layout(G, k=1, iterations=50)
    nx.draw_networkx_nodes(G, pos, node_size=1000, node_color="lightblue")
    nx.draw_networkx_edges(G, pos, alpha=0.2)
    nx.draw_networkx_labels(G, pos, font_size=8)
    plt.title("Solvent Similarity Network\n(Edges show similarities above threshold)")
    plt.axis("off")
    plt.savefig("solvent_network.png", dpi=300, bbox_inches="tight")
    plt.close()


 

def export_solvent_excel(excel_path: str = "solvent_database.xlsx") -> bool:
    """Export the solvent database to an Excel file on demand (no side effects at import).

    Returns True on success, False otherwise.
    """
    if not OPENPYXL_AVAILABLE:
        print("Note: openpyxl not installed; skipping solvent Excel export")
        return False
    try:
        with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
            # Main data sheet
            create_solvent_dataframe().to_excel(writer, sheet_name="Solvent Properties", index=False)

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
                    "Description": descriptions[: len(properties)],
                }
            )
            info_data.to_excel(writer, sheet_name="Property Information", index=False)
        print(f"Solvent database exported to {excel_path}")
        return True
    except Exception as e:
        print(f"Failed to export solvent Excel: {e}")
        return False


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
    selected_solvent = "Acetone"
    print(recommend_solvents_for_reaction(reaction_type='Cross-Coupling')[:3])
    if solvent_data is not None and not solvent_data.empty:
        print(recommend_solvents(selected_solvent, num_recommendations=3))
