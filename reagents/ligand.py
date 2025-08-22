import os
import json
import numpy as np
import pandas as pd

# Optional dependencies (graceful degradation)
try:
    from sklearn.preprocessing import MinMaxScaler  # type: ignore
except Exception:  # pragma: no cover
    MinMaxScaler = None  # type: ignore

try:
    from scipy.spatial.distance import cdist  # type: ignore
except Exception:  # pragma: no cover
    cdist = None  # type: ignore

try:
    from scipy.cluster.hierarchy import dendrogram, linkage  # type: ignore
except Exception:  # pragma: no cover
    dendrogram = None  # type: ignore
    linkage = None  # type: ignore

try:
    import matplotlib.pyplot as plt  # type: ignore
except Exception:  # pragma: no cover
    plt = None  # type: ignore

try:
    import networkx as nx  # type: ignore
except Exception:  # pragma: no cover
    nx = None  # type: ignore

try:
    import openpyxl  # noqa: F401
    OPENPYXL_AVAILABLE = True
except Exception:  # pragma: no cover
    OPENPYXL_AVAILABLE = False


EXPECTED_LIGAND_COLUMNS = [
    "Ligand",
    "Cone Angle (°)",
    "Electronic Parameter (cm⁻¹)",
    "Bite Angle (°)",
    "Steric Bulk (Å³)",
    "Donor Strength (pKa)",
    "Price Category",
    "Coordination Mode",
    "Reaction_Compatibility",
    "Typical_Applications",
]


def _default_paths():
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, 'data'))
    return (
        os.path.join(base_dir, 'ligands.json'),
        os.path.join(base_dir, 'ligands'),
    )


def create_ligand_dataframe(json_path: str | None = None, json_dir: str | None = None, prefer_builtin: bool = False):
    """
    Create the ligand DataFrame from JSON.

    Supports either:
      - data/ligands.json (array or {"ligands": [...]})
      - data/ligands/ (directory of *.json files)
    """
    default_file, default_dir = _default_paths()
    json_path = json_path or default_file
    json_dir = json_dir or default_dir

    def _normalize_entry(entry: dict) -> dict:
        name = entry.get('ligand') or entry.get('name') or entry.get('Ligand')
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
            "Ligand": name,
            "Cone Angle (°)": entry.get('cone_angle') or entry.get('Cone Angle (°)'),
            "Electronic Parameter (cm⁻¹)": entry.get('electronic_parameter') or entry.get('Electronic Parameter (cm⁻¹)'),
            "Bite Angle (°)": entry.get('bite_angle') or entry.get('Bite Angle (°)'),
            "Steric Bulk (Å³)": entry.get('steric_bulk') or entry.get('Steric Bulk (Å³)'),
            "Donor Strength (pKa)": entry.get('donor_pka') or entry.get('Donor Strength (pKa)'),
            "Price Category": entry.get('price_category') or entry.get('Price Category'),
            "Coordination Mode": entry.get('coordination_mode') or entry.get('Coordination Mode'),
            "Reaction_Compatibility": rc_str,
            "Typical_Applications": apps_str,
        }

    records: list[dict] = []
    try:
        if os.path.isfile(json_path):
            with open(json_path, 'r', encoding='utf-8') as f:
                data = json.load(f)
            if isinstance(data, dict) and 'ligands' in data:
                data = data['ligands']
            if isinstance(data, list):
                records = [_normalize_entry(d) for d in data]
        elif os.path.isdir(json_dir):
            for fn in os.listdir(json_dir):
                if fn.lower().endswith('.json'):
                    with open(os.path.join(json_dir, fn), 'r', encoding='utf-8') as f:
                        d = json.load(f)
                    if isinstance(d, dict) and 'ligand' in d and isinstance(d['ligand'], dict):
                        d = d['ligand']
                    records.append(_normalize_entry(d))
    except Exception as e:  # pragma: no cover
        print(f"Warning: failed to load ligands JSON: {e}")
        records = []

    if records:
        df = pd.DataFrame.from_records(records)
        df = df[df['Ligand'].notna() & (df['Ligand'].astype(str).str.len() > 0)]
        return df.reset_index(drop=True)

    # Fallback: empty DF with expected columns
    return pd.DataFrame(columns=EXPECTED_LIGAND_COLUMNS)


def create_feature_matrix():
    df = create_ligand_dataframe()
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
            features.append(row[col] if col in df.columns else 0)
        X.append(features)

    X = np.array(X)
    if MinMaxScaler is None:
        return X
    scaler = MinMaxScaler()
    return scaler.fit_transform(X)


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


def calculate_weighted_similarity(ligand1_features, ligand2_features, weights):
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
            max_val = max(abs(ligand1_features[i]), abs(ligand2_features[i]), 1)
            diff = abs(ligand1_features[i] - ligand2_features[i]) / max_val
            similarity += weights[feature] * (1 - diff)
    return similarity


def recommend_ligands_for_reaction(target_ligand=None, reaction_type="Cross-Coupling", top_n=5, min_compatibility=0.3):
    df = create_ligand_dataframe()

    compatible = []
    for idx, row in df.iterrows():
        compatibility = parse_reaction_compatibility(row.get("Reaction_Compatibility", ""), reaction_type)
        if compatibility >= min_compatibility:
            compatible.append({
                "index": idx,
                "name": row.get("Ligand", row.get("name", "")),
                "compatibility": compatibility,
                "applications": row.get("Typical_Applications", ""),
            })

    compatible.sort(key=lambda x: x["compatibility"], reverse=True)

    if target_ligand:
        target_idx = None
        for idx, name in enumerate(df.get("Ligand", [])):
            if isinstance(name, str) and name.lower() == str(target_ligand).lower():
                target_idx = idx
                break

        if target_idx is not None and len(compatible) > 1 and MinMaxScaler is not None:
            X = create_feature_matrix()
            if X.size and target_idx < len(X):
                target_features = X[target_idx]
                weights = REACTION_WEIGHTS.get(reaction_type, REACTION_WEIGHTS["Cross-Coupling"])
                for lig in compatible:
                    if lig["index"] != target_idx:
                        lig_features = X[lig["index"]]
                        sim = calculate_weighted_similarity(target_features, lig_features, weights)
                        lig["similarity"] = sim
                    else:
                        lig["similarity"] = 0

                compatible = [l for l in compatible if l["index"] != target_idx]
                compatible.sort(key=lambda x: (x["compatibility"] * 0.6 + x.get("similarity", 0) * 0.4), reverse=True)

    # Domain-specific adjustment: Ullmann reactions favor N-based ligands; penalize phosphines
    def _is_phosphine(name: str) -> bool:
        n = (name or "").lower()
        tokens = [
            'phos', 'phosphine', 'pph3', 'binap', 'dppf', 'dppp', 'dppe', 'xantphos', 'johnphos', 'davephos', 'pc y3', 'p(cy)3', 'ptbu', 'p(tbu)3'
        ]
        return any(t in n for t in tokens)

    def _is_n_ligand(name: str) -> bool:
        n = (name or "").lower()
        tokens = [
            'phen', 'phenanthroline', 'bipy', "bipyridine", 'proline', 'en', 'ethylenediamine', 'dmeda', 'diamine', 'pyridine'
        ]
        return any(t in n for t in tokens)

    if isinstance(reaction_type, str) and 'ullmann' in reaction_type.lower():
        for lig in compatible:
            nm = lig.get('name', '')
            boost = 0.0
            if _is_n_ligand(nm):
                boost += 0.2
            if _is_phosphine(nm):
                boost -= 0.2
            lig['adjusted'] = max(0.0, min(1.0, lig['compatibility'] + boost))
        # Sort by adjusted score if present
        compatible.sort(key=lambda x: x.get('adjusted', x['compatibility']), reverse=True)
    else:
        compatible.sort(key=lambda x: x["compatibility"], reverse=True)

    recs = []
    for i, lig in enumerate(compatible[:top_n]):
        rec = {
            "rank": i + 1,
            "ligand": lig["name"],
            "compatibility_score": round(lig.get('adjusted', lig["compatibility"]), 3),
            "applications": lig["applications"],
            "reaction_suitability": reaction_type,
        }
        if "similarity" in lig:
            rec["similarity_score"] = round(lig["similarity"], 3)
            base_comp = lig.get('adjusted', lig["compatibility"])
            rec["combined_score"] = round(base_comp * 0.6 + lig.get("similarity", 0) * 0.4, 3)
        recs.append(rec)

    # Ensure presence of at least some N-ligands for Ullmann
    if isinstance(reaction_type, str) and 'ullmann' in reaction_type.lower():
        n_count = sum(1 for r in recs if _is_n_ligand(r.get('ligand', '')))
        if n_count < max(1, top_n // 2):
            preferred = ["1,10-Phenanthroline", "2,2'-Bipyridine", "L-Proline", "Ethylenediamine", "DMEDA"]
            existing_names = {r['ligand'] for r in recs}
            for name in preferred:
                if name in existing_names:
                    continue
                # Only add if present in DB (to keep consistency), otherwise add curated placeholder
                row = df[df['Ligand'] == name]
                if not row.empty:
                    recs.append({
                        "rank": len(recs) + 1,
                        "ligand": name,
                        "compatibility_score": 0.8,
                        "applications": row.iloc[0].get("Typical_Applications", "Ullmann Cu-catalyzed"),
                        "reaction_suitability": reaction_type,
                    })
                else:
                    recs.append({
                        "rank": len(recs) + 1,
                        "ligand": name,
                        "compatibility_score": 0.75,
                        "applications": "Ullmann Cu-catalyzed",
                        "reaction_suitability": reaction_type,
                        "source": "curated"
                    })
                n_count += 1
                if n_count >= max(1, top_n // 2):
                    break
    return recs


def get_reaction_specific_ligands(reaction_type, property_preferences=None):
    recs = recommend_ligands_for_reaction(reaction_type=reaction_type, top_n=10, min_compatibility=0.4)

    if property_preferences:
        df = create_ligand_dataframe()
        filtered = []
        for rec in recs:
            ligand_name = rec["ligand"]
            matches = df[df["Ligand"] == ligand_name]
            if matches.empty:
                continue
            idx = matches.index[0]

            keep = True
            if "cone_angle_max" in property_preferences and pd.notna(df.loc[idx, "Cone Angle (°)"]):
                if df.loc[idx, "Cone Angle (°)"] > property_preferences["cone_angle_max"]:
                    keep = False
            if "price_category_max" in property_preferences and pd.notna(df.loc[idx, "Price Category"]):
                if df.loc[idx, "Price Category"] > property_preferences["price_category_max"]:
                    keep = False
            if "coordination_mode" in property_preferences and pd.notna(df.loc[idx, "Coordination Mode"]):
                if df.loc[idx, "Coordination Mode"] != property_preferences["coordination_mode"]:
                    keep = False

            if keep:
                filtered.append(rec)
        return filtered[:5]

    return recs


# Lightweight Excel export (optional)
def export_ligand_database_excel(path: str = "ligand_database.xlsx") -> None:
    if not OPENPYXL_AVAILABLE:
        print("Note: openpyxl not installed; skipping ligand Excel export")
        return
    df = create_ligand_dataframe()
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        df.to_excel(writer, sheet_name="Ligand Properties", index=False)

        info_data = pd.DataFrame(
            {
                "Property": [
                    "Cone Angle (°)",
                    "Electronic Parameter (cm⁻¹)",
                    "Bite Angle (°)",
                    "Steric Bulk (Å³)",
                    "Donor Strength (pKa)",
                    "Price Category",
                    "Coordination Mode",
                ],
                "Weight": [0.20, 0.20, 0.18, 0.15, 0.15, 0.07, 0.05],
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


# Legacy similarity-based helper using optional libs
def recommend_ligands(selected_ligand, num_recommendations=4):
    df = create_ligand_dataframe()
    if 'Ligand' not in df.columns or selected_ligand not in df['Ligand'].values:
        print(f"Ligand {selected_ligand} not found in the database.")
        return []

    props = df[[
        "Cone Angle (°)",
        "Electronic Parameter (cm⁻¹)",
        "Bite Angle (°)",
        "Steric Bulk (Å³)",
        "Donor Strength (pKa)",
        "Price Category",
        "Coordination Mode",
    ]].copy()
    # Normalize 0-1 per column
    for col in props.columns:
        min_val = props[col].min()
        max_val = props[col].max()
        if pd.isna(min_val) or pd.isna(max_val) or max_val == min_val:
            props[col] = 0
        else:
            props[col] = (props[col] - min_val) / (max_val - min_val)

    if cdist is None:
        return []

    names = df['Ligand'].values
    idx = np.where(names == selected_ligand)[0][0]
    X = props.values
    distances = cdist(X[idx:idx+1], X, metric="euclidean")[0]
    closest = np.argsort(distances)[1:num_recommendations+1]
    recs = names[closest]
    scores = 1 / (1 + distances[closest])
    return list(zip(recs, scores))


if __name__ == "__main__":
    # Simple smoke test
    print(recommend_ligands_for_reaction(reaction_type='Cross-Coupling')[:3])
