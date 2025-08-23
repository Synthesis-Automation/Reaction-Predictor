from __future__ import annotations

import json
import os
from collections import Counter
import csv
from dataclasses import asdict
from datetime import datetime
from typing import Dict, List, Tuple

from .adapters import adapt_dataset_for_type, UnifiedRow
from .normalization import canonicalize, map_base, map_solvent, map_ligand, map_metal
from dataset_registry import resolve_dataset_path


def _explode(values: List[str]) -> List[str]:
    out: List[str] = []
    for v in values:
        if not v:
            continue
        # values in Ullmann CSV often look like JSON-like lists; just keep raw token
        out.append(v)
    return out


def _normalize_field(name: str | None, kind: str) -> List[str]:
    if not name:
        return []
    # The dataset has JSON-like arrays in strings; strip brackets/quotes lightly
    s = str(name).strip()
    s = s.strip('[]')
    parts = [p.strip().strip('"').strip("'") for p in s.split(',') if p.strip()]
    normed: List[str] = []
    for p in parts if parts else [s]:
        if not p:
            continue
        if kind == 'base':
            normed.append(map_base(p))
        elif kind == 'solvent':
            normed.append(map_solvent(p))
        elif kind == 'ligand':
            normed.append(map_ligand(p))
        elif kind == 'metal':
            normed.append(map_metal(p))
        else:
            normed.append(canonicalize(p))
    # dedupe
    seen = set()
    uniq: List[str] = []
    for t in normed:
        if t and t not in seen:
            seen.add(t)
            uniq.append(t)
    return uniq


def aggregate_ullmann(rows: List[UnifiedRow]) -> Dict:
    # Filter rows with Ullmann-like type
    kept: List[UnifiedRow] = [r for r in rows if 'ullmann' in canonicalize(r.reaction_type)]

    cnt_metal = Counter()
    cnt_ligand = Counter()
    cnt_base = Counter()
    cnt_solvent = Counter()

    temps: List[float] = []
    times: List[float] = []
    yields: List[float] = []

    # Co-occurrence counters
    co_lig_sol = Counter()
    co_base_sol = Counter()
    co_cat_lig = Counter()

    for r in kept:
        metals = _normalize_field(r.metal, 'metal')
        ligands = _normalize_field(r.ligand, 'ligand')
        bases = _normalize_field(r.base, 'base')
        solvents = _normalize_field(r.solvent, 'solvent')

        for t in metals:
            cnt_metal[t] += 1
        for t in ligands:
            cnt_ligand[t] += 1
        for t in bases:
            cnt_base[t] += 1
        for t in solvents:
            cnt_solvent[t] += 1

        # Co-occurrence updates (cartesian pairs per row)
        for lg in ligands:
            for sv in solvents:
                co_lig_sol[(lg, sv)] += 1
        for bs in bases:
            for sv in solvents:
                co_base_sol[(bs, sv)] += 1
        for mt in metals:
            for lg in ligands:
                co_cat_lig[(mt, lg)] += 1

        if r.temperature_c is not None:
            temps.append(r.temperature_c)
        if r.time_h is not None:
            times.append(r.time_h)
        if r.yield_pct is not None:
            yields.append(r.yield_pct)

    def _top(counter: Counter, total: int) -> List[Dict]:
        items = []
        if total <= 0:
            return items
        for name, count in counter.most_common():
            pct = count / total
            items.append({"name": name, "count": count, "pct": round(pct, 4)})
        return items

    total_rows = len(kept)

    def _pctile(arr: List[float], q: float) -> float | None:
        if not arr:
            return None
        arr_sorted = sorted(arr)
        idx = int((len(arr_sorted) - 1) * q)
        return arr_sorted[idx]

    def _winsorize(arr: List[float], low: float = 0.05, high: float = 0.95) -> List[float]:
        if not arr:
            return arr
        s = sorted(arr)
        lo_v = _pctile(s, low)
        hi_v = _pctile(s, high)
        if lo_v is None or hi_v is None:
            return s
        out = []
        for v in s:
            if v < lo_v:
                out.append(lo_v)
            elif v > hi_v:
                out.append(hi_v)
            else:
                out.append(v)
        return out

    # Winsorize numeric arrays to reduce outlier impact
    temps_w = _winsorize(temps)
    times_w = _winsorize(times)
    yields_w = _winsorize(yields)

    numeric_stats = {
        "temperature_c": {
            "median": _pctile(temps_w, 0.5),
            "p25": _pctile(temps_w, 0.25),
            "p75": _pctile(temps_w, 0.75),
            "n": len(temps),
        },
        "time_h": {
            "median": _pctile(times_w, 0.5),
            "p25": _pctile(times_w, 0.25),
            "p75": _pctile(times_w, 0.75),
            "n": len(times),
        },
        "yield_pct": {
            "median": _pctile(yields_w, 0.5),
            "p25": _pctile(yields_w, 0.25),
            "p75": _pctile(yields_w, 0.75),
            "n": len(yields),
        },
    }

    summary = {
        "summary": {
            "total_rows": total_rows,
            "analyzed_rows": total_rows,
            "generated_at": datetime.utcnow().isoformat() + "Z",
            "dataset_versions": {},  # filled by CLI
        },
        "top": {
            "metals": _top(cnt_metal, total_rows),
            "ligands": _top(cnt_ligand, total_rows),
            "bases": _top(cnt_base, total_rows),
            "solvents": _top(cnt_solvent, total_rows),
            "additives": [],
        },
        "cooccurrence": {
            "ligand_solvent": [
                {"a": a, "b": b, "count": c, "pct": round(c / total_rows, 4) if total_rows else 0.0}
                for (a, b), c in sorted(co_lig_sol.items(), key=lambda kv: kv[1], reverse=True)
            ],
            "base_solvent": [
                {"a": a, "b": b, "count": c, "pct": round(c / total_rows, 4) if total_rows else 0.0}
                for (a, b), c in sorted(co_base_sol.items(), key=lambda kv: kv[1], reverse=True)
            ],
            "catalyst_ligand": [
                {"a": a, "b": b, "count": c, "pct": round(c / total_rows, 4) if total_rows else 0.0}
                for (a, b), c in sorted(co_cat_lig.items(), key=lambda kv: kv[1], reverse=True)
            ],
        },
        "numeric_stats": numeric_stats,
        "notes": {
            "normalizations_applied": [
                "canonicalize tokens", "synonym maps (in-module)", "split CSV list-like fields"
            ],
            "missing_columns": [],
        },
    }

    return summary


def _write_csvs(summary: Dict, version_dir: str) -> None:
    os.makedirs(version_dir, exist_ok=True)
    # Tops
    for kind in ("metals", "ligands", "bases", "solvents"):
        rows = (summary.get("top", {}) or {}).get(kind, [])
        path = os.path.join(version_dir, f"top_{kind}.csv")
        with open(path, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["name", "count", "pct"])
            for item in rows:
                w.writerow([item.get("name"), item.get("count"), item.get("pct")])

    # Co-occurrence
    co = summary.get("cooccurrence", {}) or {}
    for cname in ("ligand_solvent", "base_solvent", "catalyst_ligand"):
        rows = co.get(cname, [])
        path = os.path.join(version_dir, f"co_{cname}.csv")
        with open(path, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["a", "b", "count", "pct"])
            for item in rows:
                w.writerow([item.get("a"), item.get("b"), item.get("count"), item.get("pct")])

    # Numeric stats
    nums = summary.get("numeric_stats", {}) or {}
    path = os.path.join(version_dir, "numeric_stats.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["metric", "median", "p25", "p75", "n"])
        for metric, stats in nums.items():
            w.writerow([
                metric,
                (stats or {}).get("median"),
                (stats or {}).get("p25"),
                (stats or {}).get("p75"),
                (stats or {}).get("n"),
            ])


def run_and_export_ullmann(out_dir: str) -> str:
    path = resolve_dataset_path("Ullmann")
    if not path:
        raise FileNotFoundError("Ullmann dataset not found")

    rows = adapt_dataset_for_type("Ullmann", path)
    summary = aggregate_ullmann(rows)

    # Versioned folder and latest.json
    base_out = os.path.join(out_dir, "data", "analytics", "Ullmann")
    ts = datetime.utcnow().strftime("%Y%m%d-%H%M%S")
    version_dir = os.path.join(base_out, ts)
    os.makedirs(version_dir, exist_ok=True)

    # Write summary.json
    with open(os.path.join(version_dir, "summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    # Also write CSVs into the version folder
    _write_csvs(summary, version_dir)

    # Update latest.json
    latest_path = os.path.join(base_out, "latest.json")
    with open(latest_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    return latest_path


if __name__ == "__main__":
    here = os.path.dirname(__file__)
    root = os.path.abspath(os.path.join(here, os.pardir))
    out = run_and_export_ullmann(root)
    print(f"Wrote analytics to: {out}")
