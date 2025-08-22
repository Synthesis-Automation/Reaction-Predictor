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
    import openpyxl  # noqa: F401
    OPENPYXL_AVAILABLE = True
except Exception:  # pragma: no cover
    OPENPYXL_AVAILABLE = False


EXPECTED_BASE_COLUMNS = [
    "Base",
    "Formula",
    "Type",  # Inorganic / Organic / Superbases
    "Basicity (pKaH)",
    "Nucleophilicity Index",
    "Solubility Class",
    "Hygroscopicity",
    "Price Category",
    "Reaction_Compatibility",
    "Typical_Applications",
]


def _default_paths():
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, 'data'))
    return (
        os.path.join(base_dir, 'bases.json'),
        os.path.join(base_dir, 'bases'),
    )


def create_base_dataframe(json_path: str | None = None, json_dir: str | None = None):
    """Create the base DataFrame from JSON.

    Supports either:
      - data/bases.json (array or {"bases": [...]})
      - data/bases/ (directory of *.json files)
    """
    default_file, default_dir = _default_paths()
    json_path = json_path or default_file
    json_dir = json_dir or default_dir

    def _normalize_entry(entry: dict) -> dict:
        name = entry.get('base') or entry.get('name') or entry.get('Base')
        formula = entry.get('formula') or entry.get('Formula')
        btype = entry.get('type') or entry.get('Type')
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
            "Base": name,
            "Formula": formula,
            "Type": btype,
            "Basicity (pKaH)": entry.get('basicity_pkah') or entry.get('Basicity (pKaH)'),
            "Nucleophilicity Index": entry.get('nucleophilicity_index') or entry.get('Nucleophilicity Index'),
            "Solubility Class": entry.get('solubility_class') or entry.get('Solubility Class'),
            "Hygroscopicity": entry.get('hygroscopicity') or entry.get('Hygroscopicity'),
            "Price Category": entry.get('price_category') or entry.get('Price Category'),
            "Reaction_Compatibility": rc_str,
            "Typical_Applications": apps_str,
        }

    records: list[dict] = []
    try:
        if os.path.isfile(json_path):
            with open(json_path, 'r', encoding='utf-8') as f:
                data = json.load(f)
            if isinstance(data, dict) and 'bases' in data:
                data = data['bases']
            if isinstance(data, list):
                records = [_normalize_entry(d) for d in data]
        elif os.path.isdir(json_dir):
            for fn in os.listdir(json_dir):
                if fn.lower().endswith('.json'):
                    with open(os.path.join(json_dir, fn), 'r', encoding='utf-8') as f:
                        d = json.load(f)
                    if isinstance(d, dict) and 'base' in d and isinstance(d['base'], dict):
                        d = d['base']
                    records.append(_normalize_entry(d))
    except Exception as e:  # pragma: no cover
        print(f"Warning: failed to load bases JSON: {e}")
        records = []

    if records:
        df = pd.DataFrame.from_records(records)
        df = df[df['Base'].notna() & (df['Base'].astype(str).str.len() > 0)]
        return df.reset_index(drop=True)

    return pd.DataFrame(columns=EXPECTED_BASE_COLUMNS)


def create_base_feature_matrix():
    df = create_base_dataframe()
    feature_columns = [
        "Basicity (pKaH)",
        "Nucleophilicity Index",
        # Map categorical Solubility Class to ordinal via simple mapping if present
    ]

    # Build numerical features with fallbacks
    X = []
    for _, row in df.iterrows():
        features = []
        # Basicity
        features.append(row.get("Basicity (pKaH)", 0) if "Basicity (pKaH)" in df.columns else 0)
        # Nucleophilicity
        features.append(row.get("Nucleophilicity Index", 0) if "Nucleophilicity Index" in df.columns else 0)
        X.append(features)

    X = np.array(X)
    if MinMaxScaler is None:
        return X
    scaler = MinMaxScaler()
    return scaler.fit_transform(X)


BASE_REACTION_WEIGHTS = {
    "Cross-Coupling": {
        "Basicity (pKaH)": 0.35,
        "Nucleophilicity Index": 0.25,
    },
    "Hydrogenation": {
        "Basicity (pKaH)": 0.10,
        "Nucleophilicity Index": 0.10,
    },
    "Metathesis": {
        "Basicity (pKaH)": 0.15,
        "Nucleophilicity Index": 0.10,
    },
    "C-H_Activation": {
        "Basicity (pKaH)": 0.25,
        "Nucleophilicity Index": 0.20,
    },
    "Carbonylation": {
        "Basicity (pKaH)": 0.20,
        "Nucleophilicity Index": 0.15,
    },
}


def parse_base_reaction_compatibility(compatibility_str, reaction_type):
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


def _weighted_similarity(v1, v2, weights: dict[str, float]):
    sim = 0.0
    names = ["Basicity (pKaH)", "Nucleophilicity Index"]
    for i, nm in enumerate(names):
        if nm in weights:
            max_val = max(abs(v1[i]), abs(v2[i]), 1)
            diff = abs(v1[i] - v2[i]) / max_val
            sim += weights[nm] * (1 - diff)
    return sim


def recommend_bases_for_reaction(target_base=None, reaction_type="Cross-Coupling", top_n=5, min_compatibility=0.3):
    df = create_base_dataframe()

    candidates = []
    for idx, row in df.iterrows():
        compatibility = parse_base_reaction_compatibility(row.get("Reaction_Compatibility", ""), reaction_type)
        if compatibility >= min_compatibility:
            candidates.append({
                "index": idx,
                "name": row.get("Base", row.get("name", "")),
                "compatibility": compatibility,
                "applications": row.get("Typical_Applications", ""),
                "type": row.get("Type", ""),
            })

    candidates.sort(key=lambda x: x["compatibility"], reverse=True)

    # Optional similarity re-ranking if a target base is provided
    if target_base:
        target_idx = None
        if 'Base' in df.columns:
            for i, nm in enumerate(df['Base']):
                if isinstance(nm, str) and nm.lower() == str(target_base).lower():
                    target_idx = i
                    break
        if target_idx is not None and len(candidates) > 1 and MinMaxScaler is not None and cdist is not None:
            X = create_base_feature_matrix()
            if X.size and target_idx < len(X):
                tgt = X[target_idx]
                weights = BASE_REACTION_WEIGHTS.get(reaction_type, BASE_REACTION_WEIGHTS["Cross-Coupling"])
                for base in candidates:
                    if base["index"] != target_idx:
                        v = X[base["index"]]
                        base["similarity"] = _weighted_similarity(tgt, v, weights)
                    else:
                        base["similarity"] = 0
                candidates = [b for b in candidates if b["index"] != target_idx]
                candidates.sort(key=lambda x: (x["compatibility"] * 0.6 + x.get("similarity", 0) * 0.4), reverse=True)

    # Domain-specific adjustment: Ullmann boosts for K2CO3, Cs2CO3, K3PO4, KOtBu
    def _is_ullmann_favored(name: str) -> bool:
        n = (name or "").lower()
        tokens = ["k2co3", "cs2co3", "k3po4", "kotbu", "naotbu", "potassium carbonate", "cesium carbonate", "potassium phosphate"]
        return any(t in n for t in tokens)

    if isinstance(reaction_type, str) and 'ullmann' in reaction_type.lower():
        for b in candidates:
            nm = b.get('name', '')
            boost = 0.0
            if _is_ullmann_favored(nm):
                boost += 0.2
            b['adjusted'] = max(0.0, min(1.0, b['compatibility'] + boost))
        candidates.sort(key=lambda x: x.get('adjusted', x['compatibility']), reverse=True)

    recs = []
    for i, b in enumerate(candidates[:top_n]):
        rec = {
            "rank": i + 1,
            "base": b["name"],
            "compatibility_score": round(b.get('adjusted', b["compatibility"]), 3),
            "applications": b["applications"],
            "reaction_suitability": reaction_type,
            "type": b.get("type", "")
        }
        if "similarity" in b:
            rec["similarity_score"] = round(b["similarity"], 3)
            base_comp = b.get('adjusted', b["compatibility"])
            rec["combined_score"] = round(base_comp * 0.6 + b.get("similarity", 0) * 0.4, 3)
        recs.append(rec)

    # If no data present, provide a minimal curated fallback for Ullmann
    if not recs and isinstance(reaction_type, str) and 'ullmann' in reaction_type.lower():
        recs = [
            {"rank": 1, "base": "K2CO3", "compatibility_score": 0.8, "applications": "Ullmann Cu-catalyzed", "reaction_suitability": reaction_type, "source": "curated"},
            {"rank": 2, "base": "Cs2CO3", "compatibility_score": 0.78, "applications": "Ullmann Cu-catalyzed", "reaction_suitability": reaction_type, "source": "curated"},
            {"rank": 3, "base": "K3PO4", "compatibility_score": 0.75, "applications": "Ullmann Cu-catalyzed", "reaction_suitability": reaction_type, "source": "curated"},
            {"rank": 4, "base": "KOtBu", "compatibility_score": 0.72, "applications": "Ullmann Cu-catalyzed", "reaction_suitability": reaction_type, "source": "curated"},
        ][:top_n]

    return recs


def get_reaction_specific_bases(reaction_type, property_preferences=None):
    recs = recommend_bases_for_reaction(reaction_type=reaction_type, top_n=10, min_compatibility=0.4)

    if property_preferences:
        df = create_base_dataframe()
        filtered = []
        for rec in recs:
            base_name = rec["base"]
            matches = df[df["Base"] == base_name]
            if matches.empty:
                continue
            idx = matches.index[0]

            keep = True
            # Example filters
            if "pkah_min" in property_preferences and pd.notna(df.loc[idx, "Basicity (pKaH)"]):
                if df.loc[idx, "Basicity (pKaH)"] < property_preferences["pkah_min"]:
                    keep = False
            if "type_in" in property_preferences and pd.notna(df.loc[idx, "Type"]):
                if str(df.loc[idx, "Type"]) not in set(property_preferences["type_in"]):
                    keep = False

            if keep:
                filtered.append(rec)
        return filtered[:5]

    return recs


def export_base_database_excel(path: str = "base_database.xlsx") -> None:
    if not OPENPYXL_AVAILABLE:
        print("Note: openpyxl not installed; skipping base Excel export")
        return
    df = create_base_dataframe()
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        df.to_excel(writer, sheet_name="Base Properties", index=False)

        info_data = pd.DataFrame(
            {
                "Property": [
                    "Basicity (pKaH)",
                    "Nucleophilicity Index",
                    "Type",
                ],
                "Weight": [0.3, 0.2, 0.0],
                "Description": [
                    "Conjugate acid pKa (higher = stronger base)",
                    "Relative nucleophilicity index",
                    "Inorganic/Organic category",
                ],
            }
        )
        info_data.to_excel(writer, sheet_name="Property Information", index=False)


if __name__ == "__main__":
    # Simple smoke
    print(recommend_bases_for_reaction(reaction_type='Cross-Coupling')[:3])
