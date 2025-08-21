# JSON data formats for ligands and solvents

You can provide editable JSON instead of modifying Python tables.

Preferred locations:
- data/ligands.json or data/ligands/ (one JSON per ligand)
- data/solvents.json or data/solvents/ (one JSON per solvent)

Examples:

Ligand (object fields):
{
  "ligand": "DPPF",
  "cone_angle": 125,
  "electronic_parameter": 2072.1,
  "bite_angle": 96,
  "steric_bulk": 425,
  "donor_pka": 2.80,
  "price_category": 2,
  "coordination_mode": 2,
  "reaction_compatibility": {
    "Cross-Coupling": 0.8,
    "Hydrogenation": 0.8,
    "Metathesis": 0.3,
    "C-H_Activation": 0.6,
    "Carbonylation": 0.8
  },
  "typical_applications": ["Cross-coupling", "carbonylation"]
}

Solvent (object fields):
{
  "solvent": "DMF",
  "abbreviation": "DMF",
  "cas": "68-12-2",
  "dielectric_constant": 36.7,
  "polarity_index": 6.4,
  "boiling_point_c": 153,
  "density_g_ml": 0.944,
  "dipole_moment_d": 3.82,
  "donor_number_dn": 26.6,
  "hydrogen_bond_donor": 0,
  "reaction_compatibility": [0.9, 0.6, 0.2, 0.6, 0.8],
  "typical_applications": ["Polar aprotic solvent", "cross-coupling"]
}

Notes:
- Arrays or dicts for reaction_compatibility are accepted; they’ll be normalized internally.
- If both JSON and built-in tables exist, JSON takes precedence.
- Per-file JSON: put files under data/ligands/ or data/solvents/ with any name ending in .json.
