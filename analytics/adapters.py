from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Any
import csv
import os

from .normalization import canonicalize


@dataclass
class UnifiedRow:
    reaction_type: str
    metal: str | None
    catalyst: str | None
    ligand: str | None
    base: str | None
    solvent: str | None
    additives: List[str]
    temperature_c: float | None
    time_h: float | None
    cat_loading_molpct: float | None
    base_equiv: float | None
    yield_pct: float | None
    source: str | None


ULLMANN_COLUMN_MAP = {
    # source columns -> unified fields
    "ReactionType": "reaction_type",
    "CoreGeneric": "metal",  # often Cu/Cu(I)/Cu(II)
    "Ligand": "ligand",
    "ReagentRaw": "base",  # may include BASE etc.; we will postprocess
    "Solvent": "solvent",
    "Temperature_C": "temperature_c",
    "Time_h": "time_h",
    "Yield_%": "yield_pct",
    "Reference": "source",
}


def _read_csv_rows(path: str) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    with open(path, newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for r in reader:
            rows.append(r)
    return rows


def adapt_ullmann(path: str) -> List[UnifiedRow]:
    raw = _read_csv_rows(path)
    out: List[UnifiedRow] = []

    for r in raw:
        rt = (r.get("ReactionType") or "").strip() or "Ullmann"
        metal = (r.get("CoreGeneric") or r.get("CoreDetail") or "").strip()
        ligand = (r.get("Ligand") or "").strip()
        base = (r.get("ReagentRaw") or "").strip()
        solvent = (r.get("Solvent") or "").strip()
        temp = r.get("Temperature_C")
        time = r.get("Time_h")
        yld = r.get("Yield_%")
        source = (r.get("Reference") or "").strip()

        def _as_float(x):
            try:
                return float(x) if x not in (None, "", [], {}) else None
            except Exception:
                return None

        row = UnifiedRow(
            reaction_type=rt,
            metal=metal or None,
            catalyst=None,
            ligand=ligand or None,
            base=base or None,
            solvent=solvent or None,
            additives=[],
            temperature_c=_as_float(temp),
            time_h=_as_float(time),
            cat_loading_molpct=None,
            base_equiv=None,
            yield_pct=_as_float(yld),
            source=source or None,
        )
        out.append(row)

    return out


def adapt_dataset_for_type(reaction_type: str, path: str) -> List[UnifiedRow]:
    # Only Ullmann in Milestone 1
    if canonicalize(reaction_type).startswith("ullmann"):
        return adapt_ullmann(path)
    raise ValueError(f"No adapter for reaction type: {reaction_type}")
