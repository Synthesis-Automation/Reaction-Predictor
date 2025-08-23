# Dataset Analytics and Evidence Priors — Implementation Plan

Date: 2025-08-23

## Goals

- Build a robust, repeatable analysis over reaction datasets (e.g., Ullmann) to surface frequently used metals/catalysts, ligands, bases, solvents, and key numeric conditions.
- Persist analytics to versioned artifacts and a `latest.json` summary per reaction type.
- Use these summaries as conservative evidence priors to gently boost recommendations in the enhanced engine.

## Scope

- Phase 1 focuses on Ullmann; design is generic to extend to other reaction types.
- Read-only analysis of CSV datasets in `data/reaction_dataset/`; no changes to raw data.
- Non-invasive boosts (bounded) that augment existing compatibility scoring.

## Deliverables

- Normalization utilities for reagents and numeric condition parsing.
- Dataset adapters to unify columns into a shared schema.
- Aggregation pipeline producing top lists, co-occurrences, and numeric stats.
- Versioned exports under `data/analytics/<reaction_type>/...` plus `latest.json`.
- Engine hook to load `latest.json` and apply frequency-based priors.
- CLI script `scripts/analyze_dataset.py` and a VS Code task to run it.
- Minimal unit tests and a smoke validation.

## Architecture Overview

### 1. Normalization (analytics/normalization.py)

- Token canonicalization: casefold, strip punctuation/whitespace, unify unicode.
- Synonym maps (data/normalization/*.yaml):
  - Bases: KOtBu/t-BuOK, K2CO3, Cs2CO3, K3PO4, Na2CO3, Et3N/TEA, DIPEA, DBU.
  - Solvents: DMF, NMP (N-Methyl-2-pyrrolidone), MeTHF/2-MeTHF, i-PrOH/IPA, DMSO, MeCN/ACN, Toluene/PhMe.
  - Ligands: BINAP variants, XPhos/SPhos/BrettPhos, PPh3/P(p-tol)3/P(o-tol)3, "none".
  - Metals/catalysts: Cu/CuI/CuBr, Pd(PPh3)4→Pd+PPh3, Ni(cod)2→Ni.
- Parsers: split mixtures (e.g., "DMF/THF", "toluene:MeCN"), strip descriptors ("dry"), extract numbers and units (°C, h, mol%, equiv).

### 2. Adapters (analytics/adapters.py)

- Map dataset-specific columns to a unified schema:
  - reaction_type, metal, catalyst, ligand, base, solvent, additives,
    temperature_c, time_h, cat_loading_molpct, base_equiv, yield_pct, source.
- Integrate with `dataset_registry.py` to enumerate CSVs and provide column maps.
- Heuristic fallback if columns missing; record in `notes.missing_columns`.

### 3. Aggregation (analytics/aggregate.py)

- Filter by normalized reaction_type.
- Normalize fields and explode multi-values.
- Counts: top metals, ligands, bases, solvents, additives (with count, proportion).
- Co-occurrence matrices: ligand–solvent, base–solvent, catalyst–ligand.
- Numeric stats with robust parsing and light winsorization (1st–99th pct):
  - temperature_c, time_h, cat_loading_molpct, base_equiv, yield_pct.
- Support thresholds (e.g., ≥1–2%) for surfacing; keep raw counts in JSON.

### 4. Export & Versioning

- Write to `data/analytics/<reaction_type>/<YYYYMMDD-HHMMSS>/`:
  - summary.json (master), top_*.csv, co_*.csv, numeric_stats.csv.
- Maintain `data/analytics/<reaction_type>/latest.json` for fast engine loading.
- Include dataset file mtimes and a content hash; optional cache skip when unchanged.

### 5. Engine Integration (enhanced_recommendation_engine.py)

- `_load_analytics_summary(reaction_type)` loads `latest.json` if present.
- Convert top lists to frequency-based priors and apply conservative boosts:
  - final_score = base_score × (1 + w_freq × g(freq_norm)), cap to [0,1].
  - Defaults: w_freq=0.15, g(x)=sqrt(x); ignore items below min support.
- Apply to solvents, bases, ligands, and optional catalyst/metal.
- Re-rank post-boost; fall back to current behavior when no summary exists.
- Config toggles per reagent type for priors.

## Data Contracts

Inputs

- CSVs under `data/reaction_dataset/` as registered via `dataset_registry.py`.

Outputs

- JSON: `data/analytics/<reaction_type>/latest.json` with:
  - summary: { total_rows, analyzed_rows, generated_at, dataset_versions }
  - top: { metals[], ligands[], bases[], solvents[], additives[] } items { name, count, pct }
  - cooccurrence: { ligand_solvent[], base_solvent[], catalyst_ligand[] } items { a, b, count, pct }
  - numeric_stats: { temperature_c, time_h, cat_loading_molpct, base_equiv, yield_pct } each with { median, p25, p75, n }
  - notes: { normalizations_applied[], missing_columns[] }
- CSVs: top_*.csv, co_*.csv, numeric_stats.csv in the same folder.

Error Modes

- Missing columns: skip affected metrics, record in notes, continue.
- Unparseable numeric values: drop and count as missing.
- Empty or low-support datasets: produce empty top lists; engine will fallback.

Success Criteria

- `latest.json` exists with non-empty top lists for Ullmann.
- Engine loads priors without errors and gently re-ranks results.

## Implementation Steps (Milestones)

Milestone 1: Foundations (Day 1)

- Create `analytics/normalization.py` with core token normalization and synonym maps (YAML).
- Add `analytics/adapters.py` for Ullmann dataset -> unified schema.
- Add `analytics/aggregate.py` to compute top lists and numeric stats.
- Export versioned artifacts and `latest.json` for Ullmann.

Milestone 2: Engine Hook (Day 2)

- Implement `_load_analytics_summary()` and apply priors to solvents and bases.
- Add config weights/gates; cap scores; re-rank.
- Smoke validation showing expected gentle boosts for common Ullmann items (DMF, K2CO3, etc.).

Milestone 3: Coverage & Co-occurrence (Day 3)

- Add co-occurrence tables; extend priors to ligands and catalyst/metal.
- Improve numeric parsers and winsorization.
- Unit tests for normalization and parsers; small synthetic aggregator test.

Milestone 4: CLI & Tasks (Day 3–4)

- New script `scripts/analyze_dataset.py`:
  - `--reaction-type Ullmann` or `--all`, `--min-support`, `--no-cache`, `--out-dir`.
- Add VS Code task: "Analyze: Ullmann (PYTHONPATH)".
- README section: dataset analytics usage and troubleshooting.

## Validation & Tests

- Unit tests (pytest or simple asserts) for:
  - Synonym mapping (bases/solvents/ligands) and multi-solvent splitting.
  - Numeric parsing and unit conversion.
  - Aggregator counts on a tiny fixture.
- Smokes:
  - Run analyzer for Ullmann, confirm `latest.json` created and non-empty top lists.
  - Run recommendation pre/post analytics; confirm boosted ordering for frequent items.

## Risks & Edge Cases

- Mixed/ambiguous naming: mitigate with evolving synonym maps (YAML) and logging unknowns.
- Multi-component solvents and ratios: treat each as present; optionally weight by first/ratio if parsed.
- Dataset drift: adapters isolate column-map changes; record missing columns in notes.
- Over-boosting: keep weights conservative and cap final scores to [0,1].

## Runbook (once implemented)

- Analyze Ullmann only:

```powershell
$env:PYTHONPATH='.'; python scripts/analyze_dataset.py --reaction-type "Ullmann"
```

- Analyze all known types:

```powershell
$env:PYTHONPATH='.'; python scripts/analyze_dataset.py --all
```

- Use in recommendations: the engine will auto-load `data/analytics/<type>/latest.json` when present.

## Ownership & Next Steps

- Owner: Reaction-Predictor maintainers.
- Next: Implement Milestone 1; iterate synonyms as new datasets arrive; extend to additional reaction types.
