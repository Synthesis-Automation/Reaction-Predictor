
# QUARC‑OSS Integration Plan (Reagents/Catalysts First)

This document is a step‑by‑step plan to integrate QUARC‑OSS into our recommendation system, focused on reagents/catalysts first. Temperature and equivalents will be integrated later. It also embeds an iterative chemist feedback loop to continuously improve quality.

## Goals
- Add data‑driven agent (reagents/catalysts) recommendations via QUARC‑OSS (FFN pipeline first).
- Keep current heuristics/analytics as fallback and allow hybrid re‑ranking.
- Preserve current JSON export schema; no breaking changes.
- Make the system testable (automated metrics) and reviewable (chemist grading pack).

## Scope for Phase 1 (Now)
- Agents only: catalyst, ligand, base (solvent optional; keep current solvent behavior initially).
- Use QUARC‑OSS pretrained pipeline; no NameRxn/Pistachio dependency.
- No retraining yet; inference‑only with analytics re‑ranking; caching; CLI/engine flags.

---

## Phase 0 — Preparation (no behavior change)

Deliverables:
- Env toggles and flags (declared; no-op initially).
- Controlled vocabulary for agents (roles + synonyms + CAS).

Steps:
1) Choose QUARC flavor
   - Default: QUARC‑OSS (no NameRxn). Defer GNN until atom‑mapping available.
2) Config and env
   - Env vars (document):
     - `QUARC_OSS_HOME` (path to cloned quarc-oss)
     - `QUARC_OSS_CONFIG` (pipeline YAML, e.g., ffn_pipeline.yaml)
     - `QUARC_OSS_PYTHON` (python in its conda env)
     - `QUARC_OSS_CHECKPOINTS` (optional override)
   - CLI/Engine flags (declare only now):
     - `--use-quarc [auto|always|off]` (default: auto)
     - `--quarc-config <path>`
     - `--quarc-topk <N>`
3) Vocabulary seed files
   - Build `data/vocab/agents_vocab.json`: list of { name, role, synonyms[], cas? }.
   - Use `data/cas_dictionary.csv`, `data/ligands.json`, `data/solvents.json` for enrichment.
   - Roles: catalyst, ligand, base, solvent, additive.

---

## Phase 1 — Adapter + Engine Wiring (agents only)

Deliverables:
- `integration/quarc_oss_adapter.py` (inference only; agents)
- Engine hook to call adapter and merge agents
- Analytics re‑ranking + soft penalty
- Caching and graceful fallback
- CLI flags operational
- README section “Using QUARC‑OSS (agents only)”

Steps:
1) Implement adapter
   - Contract:
     - Input: `{ rxn_smiles: str, top_k: int=5 }`
     - Output: `[{ name: str, role?: str, score?: float }]`
     - Errors: missing models/env/timeout → return `[]` and structured error object
   - Invocation:
     - Prefer Python API if quarc‑oss exposes one; else subprocess `scripts/inference.py` with `--config-path`, `--input` temp JSON, `--output`, `--top-k`.
     - Timeout (e.g., 30s). Log stderr (truncate) and classify errors.
   - Post‑processing:
     - Decode agent tokens → human names (use quarc‑oss encoder artifacts if provided).
     - Assign roles using vocabulary heuristics (regex/synonyms).
     - Keep catalyst/ligand/base; treat others as additive (ignored for now).
2) Engine integration (`enhanced_recommendation_engine.py`)
   - Add `_get_quarc_agents(smiles, k)`; guarded by `--use-quarc` and env presence.
   - Merge:
     - Promote returned catalyst/ligand/base into our lists; de‑duplicate by canonical name.
     - Re‑rank with analytics priors for the detected reaction family (e.g., Ullmann) and apply soft penalties for unsupported items.
     - Keep solvent from current engine/analytics initially.
   - Update `combined_conditions` to prioritize QUARC‑backed combos.
   - Add `meta.providers: ["quarc_oss","analytics"]` in exports when used.
3) Caching
   - Cache by normalized SMILES: `.cache/quarc_oss/<hash>.json` (input, predictions, created_at, adapter_version).
   - Evict LRU or cap size (e.g., last 500 entries).
4) CLI flags
   - `predict_cli.py`: parse `--use-quarc`, `--quarc-config`, `--quarc-topk`, forward to engine.
   - Keep output UTF‑8 safe and support `-o` writing.
5) Docs
   - README: how to set envs, create conda env for quarc‑oss, and run with flags on Windows (PowerShell examples).

Acceptance:
- If quarc‑oss not configured, behavior unchanged.
- If configured, recommendations include quarc‑backed agents; export remains schema‑compatible.

---

## Phase 2 — Evaluation + 20‑Reaction Chemist Grading Pack

Deliverables:
- `scripts/evaluate_recommendations.py` (automated metrics for agents)
- `scripts/grading/make_panel.py` (sample 20 diverse reactions; export XLSX/CSV/HTML pack)
- VS Code tasks to run both

Steps:
1) Automated eval (agents only)
   - Build gold from `data/reaction_dataset/*`; normalize roles using the vocabulary.
   - Metrics per role: Top‑K accuracy, MRR@K, Recall@K, Jaccard@K; combined triple match@K (optional).
   - Stratify by reaction family; bootstrap 95% CIs.
2) Human grading (n=20)
   - Sampling: stratify (e.g., 5 Ullmann, 5 amide, 5 Buchwald, 5 other), ensure diversity via ECFP4 Tanimoto (< 0.7), fixed seed.
   - Emit sheet with top‑3 per role and rubric (1–10); include notes column.
   - Collector script to merge grades and compute mean/median and acceptable@K (≥7/10).

Acceptance:
- Reports produced under `reports/`, with summary JSON and per‑reaction CSV.

---

## Phase 3 — Feedback Loop (Weight Tuning + Calibration)

Deliverables:
- `data/feedback/grades.csv` schema: reaction_id, smiles, type, candidates JSON, grade, reviewer, notes, timestamp, engine_version.
- `scripts/feedback/ingest_grades.py` (validate & merge)
- `scripts/feedback/optimize_weights.py` (Optuna)
- `weights/engine_weights.json` (w_quarc vs w_analytics per role, penalty factors)
- Optional: Platt/Isotonic calibration artifacts `calibration/*.pkl`

Steps:
1) Ingest & label
   - Convert grades to targets: accept (grade ≥ 7), preference pairs (higher > lower), or normalized grade.
2) Tune weights
   - Objective: maximize accept@3 or mean grade on validation.
   - Params: `w_freq_*`, `penalty_factor_*`, synergy bonuses, `w_quarc` per role.
3) Calibrate scores
   - Fit calibration to map scores → P(accept≥7); show confidence bands in UI (later).
4) Deploy weights
   - Engine loads `weights/engine_weights.json` if present; versioned in export meta.

Acceptance:
- +5–10% accept@3 vs baseline in validation (target; iterate if unmet).

---

## Phase 4 — Learning‑to‑Rank Re‑ranker (Optional but High‑Impact)

Deliverables:
- `scripts/feedback/train_ranker.py` → `models/ranker_{role}.pkl`
- Engine re‑ranking hook per role

Steps:
1) Features per candidate: engine score, analytics pct, quarc score, role confidence, synonym/CAS presence; simple chemistry flags.
2) Train LambdaMART/XGBoost pairwise using preference pairs (from grades) or graded relevance.
3) Inference: re‑rank top‑K per role before building `combined_conditions`.

Acceptance:
- +10–15% accept@3 uplift on validation set.

---

## Phase 5 — Domain Adaptation (Train QUARC‑OSS on Our Data) [Deferred]

Deliverables:
- `scripts/datasets/convert_to_quarc_oss.py` (create training records and agent encoder)
- Fine‑tuned checkpoints under `data/quarc_oss_checkpoints/`

Steps:
1) Convert `data/reaction_dataset/*` into QUARC‑OSS training samples (agents only).
2) Train/fine‑tune FFN stages; evaluate, calibrate; point adapter to tuned pipeline.

---

## Phase 6 — Temperature & Equivalents [Deferred]

- Integrate QUARC predictions for temperature and equivalents; extend export (`conditions`) with numeric values/ranges.
- Add MAE/RMSE metrics and tolerance bands (±10 °C, ±0.2 equiv) to the eval harness.

---

## Risks & Mitigations
- Env friction on Windows → document conda setup; support subprocess via `QUARC_OSS_PYTHON`.
- Token/name mismatches → expand synonyms vocabulary; keep “other” bucket and avoid hard failures.
- Performance variance → enable caching and timeouts; graceful fallback to current engine.
- Data sparsity long‑tail → cap vocab per role; apply soft penalties using analytics support.

---

## Milestones & Exit Criteria
- M1 (Phase 1): QUARC‑OSS agents integrated; export unchanged; CLI flags live; caching; README updated.
- M2 (Phase 2): Automated eval + 20‑reaction grading pack; baseline metrics with CIs.
- M3 (Phase 3): Weight tuning & calibration; measurable uplift on validation.
- M4 (Phase 4, optional): LTR re‑ranker delivers additional uplift.

---

## Try It (once Phase 1 lands)

PowerShell examples (UTF‑8 safe):

```powershell
# Set env (example)
$env:QUARC_OSS_HOME='C:\\path\\to\\quarc-oss'
$env:QUARC_OSS_PYTHON='C:\\Miniconda3\\envs\\quarc\\python.exe'
$env:QUARC_OSS_CONFIG='C:\\path\\to\\quarc-oss\\configs\\ffn_pipeline.yaml'

# Inline JSON and write to file
$payload = '{"reaction_smiles":"Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1","selected_reaction_type":"C-N Coupling - Ullmann"}'
python predict_cli.py --use-quarc auto --quarc-topk 5 $payload -o .\\out.json
Get-Content .\\out.json | Out-String -Width 4096 | Write-Output
```

---

## File/Module Plan (to be added)
- `integration/quarc_oss_adapter.py` — inference wrapper + role mapping + caching
- `scripts/evaluate_recommendations.py` — automated metrics
- `scripts/grading/make_panel.py` — 20‑reaction grading pack
- `scripts/feedback/ingest_grades.py` — merge and validate grades
- `scripts/feedback/optimize_weights.py` — Optuna tuning
- `scripts/feedback/train_ranker.py` — train LTR re‑Ranker
- `weights/engine_weights.json` — tuned weights (loaded by engine)

---

## Ownership & Review
- Owner: Recommendation team
- Reviewer: Lead chemist (grading rubric, panel selection)
- Cadence: bi‑weekly review of metrics/grades; monthly weight tuning release

---

Last updated: 2025‑08‑25
