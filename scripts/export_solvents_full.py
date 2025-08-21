import os
import sys
import json

# Ensure project root is on sys.path so we can import the top-level 'reagents' package
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from reagents.solvent import solvent_data

ORDER = ["Cross-Coupling","Hydrogenation","Metathesis","C-H_Activation","Carbonylation"]

data_dir = os.path.join(ROOT, 'data')
os.makedirs(data_dir, exist_ok=True)
out_path = os.path.join(data_dir, 'solvents.json')

items = []
for _, row in solvent_data.iterrows():
    rc = row.get('Reaction_Compatibility')
    try:
        parts = [float(x) for x in str(rc).split(',')]
    except Exception:
        parts = [0.5]*5
    rc_dict = {k:(parts[i] if i < len(parts) else 0.5) for i,k in enumerate(ORDER)}
    apps = row.get('Typical_Applications','')
    if isinstance(apps, str):
        apps = [a.strip() for a in apps.split(',') if a.strip()]

    items.append({
        'solvent': row.get('Solvent'),
        'abbreviation': row.get('Abbreviation'),
        'cas': row.get('CAS Number'),
        'dielectric_constant': row.get('Dielectric Constant'),
        'polarity_index': row.get('Polarity Index'),
        'boiling_point_c': row.get('Boiling Point (°C)') or row.get('Boiling Point (C)'),
        'density_g_ml': row.get('Density (g/mL)'),
        'dipole_moment_d': row.get('Dipole Moment (D)'),
        'donor_number_dn': row.get('Donor Number (DN)'),
        'hydrogen_bond_donor': row.get('Hydrogen Bond Donor'),
        'reaction_compatibility': rc_dict,
        'typical_applications': apps,
    })

with open(out_path, 'w', encoding='utf-8') as f:
    json.dump({'solvents': items}, f, ensure_ascii=False, indent=2)

print(f"Wrote {len(items)} solvents to {out_path}")
