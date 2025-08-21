"""
Migrate in-code ligand and solvent datasets to JSON files for easier editing.
- Writes combined files: data/ligands.json, data/solvents.json
- Writes per-item files: data/ligands/*.json, data/solvents/*.json
"""
from __future__ import annotations
import os
import json
from typing import Dict, Any, List

# Import DataFrame builders from existing modules
from reagents.ligand import create_ligand_dataframe
from reagents.solvent import create_solvent_dataframe

ROOT = os.path.abspath(os.path.dirname(__file__))
DATA_DIR = os.path.join(ROOT, "data")
LIGANDS_DIR = os.path.join(DATA_DIR, "ligands")
SOLVENTS_DIR = os.path.join(DATA_DIR, "solvents")

REACTION_ORDER = ["Cross-Coupling", "Hydrogenation", "Metathesis", "C-H_Activation", "Carbonylation"]


def ensure_dirs() -> None:
    os.makedirs(DATA_DIR, exist_ok=True)
    os.makedirs(LIGANDS_DIR, exist_ok=True)
    os.makedirs(SOLVENTS_DIR, exist_ok=True)


def parse_rc_str(rc: str) -> Dict[str, float]:
    try:
        parts = [float(x) for x in str(rc).split(",")]
        out = {}
        for i, key in enumerate(REACTION_ORDER):
            out[key] = float(parts[i]) if i < len(parts) else 0.5
        return out
    except Exception:
        return {k: 0.5 for k in REACTION_ORDER}


def split_apps(apps: Any) -> List[str]:
    if isinstance(apps, list):
        return [str(a) for a in apps]
    if isinstance(apps, str):
        # split on comma if it looks like a list
        return [s.strip() for s in apps.split(",") if s.strip()] if "," in apps else ([apps.strip()] if apps.strip() else [])
    return []


def migrate_ligands() -> Dict[str, Any]:
    # Force reading from built-in tables to ensure full migration
    df = create_ligand_dataframe(prefer_builtin=True)
    items: List[Dict[str, Any]] = []
    for _, row in df.iterrows():
        item = {
            "ligand": row.get("Ligand"),
            "cone_angle": row.get("Cone Angle (°)"),
            "electronic_parameter": row.get("Electronic Parameter (cm⁻¹)"),
            "bite_angle": row.get("Bite Angle (°)"),
            "steric_bulk": row.get("Steric Bulk (Å³)"),
            "donor_pka": row.get("Donor Strength (pKa)"),
            "price_category": row.get("Price Category"),
            "coordination_mode": row.get("Coordination Mode"),
            "reaction_compatibility": parse_rc_str(row.get("Reaction_Compatibility", "")),
            "typical_applications": split_apps(row.get("Typical_Applications", "")),
        }
        if item["ligand"]:
            items.append(item)
    # write combined
    with open(os.path.join(DATA_DIR, "ligands.json"), "w", encoding="utf-8") as f:
        json.dump({"ligands": items}, f, ensure_ascii=False, indent=2)
    # write per-item
    for it in items:
        name = str(it["ligand"]).replace("/", "-")
        with open(os.path.join(LIGANDS_DIR, f"{name}.json"), "w", encoding="utf-8") as f:
            json.dump({"ligand": it}, f, ensure_ascii=False, indent=2)
    return {"count": len(items)}


def migrate_solvents() -> Dict[str, Any]:
    # Force reading from built-in tables to ensure full migration
    df = create_solvent_dataframe(prefer_builtin=True)
    items: List[Dict[str, Any]] = []
    for _, row in df.iterrows():
        item = {
            "solvent": row.get("Solvent"),
            "abbreviation": row.get("Abbreviation"),
            "cas": row.get("CAS Number"),
            "dielectric_constant": row.get("Dielectric Constant"),
            "polarity_index": row.get("Polarity Index"),
            "boiling_point_c": row.get("Boiling Point (°C)"),
            "density_g_ml": row.get("Density (g/mL)"),
            "dipole_moment_d": row.get("Dipole Moment (D)"),
            "donor_number_dn": row.get("Donor Number (DN)"),
            "hydrogen_bond_donor": row.get("Hydrogen Bond Donor"),
            "reaction_compatibility": parse_rc_str(row.get("Reaction_Compatibility", "")),
            "typical_applications": split_apps(row.get("Typical_Applications", "")),
        }
        if item["solvent"]:
            items.append(item)
    # combined
    with open(os.path.join(DATA_DIR, "solvents.json"), "w", encoding="utf-8") as f:
        json.dump({"solvents": items}, f, ensure_ascii=False, indent=2)
    # per-item
    for it in items:
        name = str(it["solvent"]).replace("/", "-")
        with open(os.path.join(SOLVENTS_DIR, f"{name}.json"), "w", encoding="utf-8") as f:
            json.dump({"solvent": it}, f, ensure_ascii=False, indent=2)
    return {"count": len(items)}


if __name__ == "__main__":
    ensure_dirs()
    # Backup current JSON if exists
    for fname in ("ligands.json", "solvents.json"):
        p = os.path.join(DATA_DIR, fname)
        if os.path.exists(p):
            bak = p + ".bak"
            try:
                with open(p, "rb") as src, open(bak, "wb") as dst:
                    dst.write(src.read())
                print(f"Backed up {fname} -> {fname}.bak")
            except Exception as e:
                print(f"Warning: could not backup {fname}: {e}")
    lig = migrate_ligands()
    solv = migrate_solvents()
    print(f"Migrated ligands: {lig['count']}")
    print(f"Migrated solvents: {solv['count']}")
