#!/usr/bin/env python3
"""
QUARC-OSS adapter for agent (reagents/catalysts) inference.
- Supports subprocess invocation of quarc-oss pipeline
- Returns normalized list of { name, role?, score? }
- Caches per SMILES hash under .cache/quarc_oss
- Gracefully falls back on errors (returns [], plus error dict)

Phase 1 scope: agents only (catalyst, ligand, base). Solvent stays engine-side.
"""
from __future__ import annotations

import os
import json
import time
import hashlib
import tempfile
import subprocess
from typing import Any, Dict, List, Optional, Tuple

_ROOT = os.path.abspath(os.path.dirname(__file__) + os.sep + "..")

CACHE_DIR = os.path.join(_ROOT, ".cache", "quarc_oss")
ADAPTER_VERSION = "1.0.0"

# Load vocabulary for basic role mapping
_VOCAB_PATH = os.path.join(_ROOT, "data", "vocab", "agents_vocab.json")
try:
    with open(_VOCAB_PATH, "r", encoding="utf-8") as f:
        _VOCAB = json.load(f)
except Exception:
    _VOCAB = {"items": []}

# Build quick synonym lookup
_syn_to_role: Dict[str, str] = {}
_name_to_role: Dict[str, str] = {}
for it in _VOCAB.get("items", []):
    role = (it.get("role") or "").strip().lower()
    nm = (it.get("name") or "").strip().lower()
    if nm:
        _name_to_role[nm] = role
    for syn in it.get("synonyms", []) or []:
        s = (str(syn) or "").strip().lower()
        if s:
            _syn_to_role[s] = role


def _hash_smiles(smi: str) -> str:
    return hashlib.sha256((smi or "").encode("utf-8")).hexdigest()[:16]


def _ensure_cache_dir():
    os.makedirs(CACHE_DIR, exist_ok=True)


def _cache_path(key: str) -> str:
    return os.path.join(CACHE_DIR, f"{key}.json")


def _load_cache(key: str) -> Optional[dict]:
    try:
        with open(_cache_path(key), "r", encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return None


def _save_cache(key: str, data: dict) -> None:
    try:
        _ensure_cache_dir()
        with open(_cache_path(key), "w", encoding="utf-8") as f:
            json.dump(data, f, ensure_ascii=False)
    except Exception:
        pass


def _prune_cache(max_entries: int = 500) -> None:
    try:
        if not os.path.isdir(CACHE_DIR):
            return
        files = [os.path.join(CACHE_DIR, fn) for fn in os.listdir(CACHE_DIR) if fn.endswith('.json')]
        if len(files) <= max_entries:
            return
        files.sort(key=lambda p: os.path.getmtime(p), reverse=True)
        for old in files[max_entries:]:
            try:
                os.remove(old)
            except Exception:
                continue
    except Exception:
        pass


def _env() -> dict:
    return {
        "home": os.environ.get("QUARC_OSS_HOME"),
        "python": os.environ.get("QUARC_OSS_PYTHON"),
        "config": os.environ.get("QUARC_OSS_CONFIG"),
        "checkpoints": os.environ.get("QUARC_OSS_CHECKPOINTS"),
    }


def _classify_role(name: str) -> Optional[str]:
    low = (name or "").strip().lower()
    if low in _name_to_role:
        return _name_to_role[low]
    if low in _syn_to_role:
        return _syn_to_role[low]
    # simple heuristics
    if any(tok in low for tok in ["phos", "binap", "bipy", "phen", "dp"]):
        return "ligand"
    if any(tok in low for tok in ["k2co3", "cs2co3", "kotbu", "k3po4", "na2co3", "et3n", "dipea", "dbu", "koh"]):
        return "base"
    if any(tok in low for tok in ["cu", "palladium", "pdoac", "pd("]):
        return "catalyst"
    return None


def run_inference(rxn_smiles: str, top_k: int = 5, config_path: Optional[str] = None, timeout_s: int = 30) -> Tuple[List[Dict], Optional[Dict]]:
    """Invoke quarc-oss to get agent suggestions.

    Returns (agents, error). On success, error is None. On failure, agents=[], error has details.
    """
    key = _hash_smiles(rxn_smiles)
    cache = _load_cache(key)
    if cache and isinstance(cache, dict) and cache.get("adapter_version") == ADAPTER_VERSION:
        return cache.get("agents", []), None

    env = _env()
    py = env.get("python") or "python"
    cfg = config_path or env.get("config")
    home = env.get("home")

    if not home or not os.path.isdir(home) or not cfg or not os.path.exists(cfg):
        # Not configured; graceful no-op
        return [], {"type": "not_configured", "message": "QUARC-OSS not configured (set QUARC_OSS_HOME and QUARC_OSS_CONFIG)"}

    with tempfile.TemporaryDirectory() as td:
        inp = os.path.join(td, "input.json")
        outp = os.path.join(td, "output.json")
        payload = {"rxn_smiles": rxn_smiles, "top_k": int(top_k)}
        with open(inp, "w", encoding="utf-8") as f:
            json.dump(payload, f, ensure_ascii=False)
        # Assume a generic quarc-oss inference script; adapt path if needed later
        script = os.path.join(home, "scripts", "inference.py")
        if not os.path.exists(script):
            return [], {"type": "missing_script", "message": f"inference.py not found under {home}"}
        cmd = [py, script, "--config-path", cfg, "--input", inp, "--output", outp, "--top-k", str(int(top_k))]
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout_s)
        except subprocess.TimeoutExpired:
            return [], {"type": "timeout", "message": f"QUARC inference timed out after {timeout_s}s"}
        if proc.returncode != 0:
            msg = (proc.stderr or proc.stdout or "error").splitlines()[:10]
            return [], {"type": "process_error", "message": "\n".join(msg)}
        try:
            with open(outp, "r", encoding="utf-8") as f:
                raw = json.load(f)
        except Exception as e:
            return [], {"type": "bad_output", "message": f"Failed to read output: {e}"}

    # Expected raw format: { agents: [ { token/name, score?, role? }, ... ] }
    raw_agents = (raw.get("agents") if isinstance(raw, dict) else None) or []
    agents: List[Dict[str, Any]] = []
    for item in raw_agents:
        try:
            name = item.get("name") or item.get("token") or ""
            name = str(name).strip()
            if not name:
                continue
            role = item.get("role") or _classify_role(name)
            score = item.get("score")
            agents.append({"name": name, "role": role, "score": float(score) if score is not None else None})
        except Exception:
            continue

    record = {
        "adapter_version": ADAPTER_VERSION,
        "created_at": time.strftime('%Y-%m-%dT%H:%M:%S'),
        "input": {"rxn_smiles": rxn_smiles, "top_k": int(top_k)},
        "agents": agents,
    }
    _save_cache(key, record)
    _prune_cache(500)
    return agents, None


def get_agents_for_engine(rxn_smiles: str, top_k: int = 5, config_path: Optional[str] = None) -> Tuple[List[Dict], Optional[Dict]]:
    """Engine-facing wrapper; hides timeout and returns normalized list.
    """
    return run_inference(rxn_smiles, top_k=top_k, config_path=config_path, timeout_s=30)


if __name__ == "__main__":
    # Simple smoke
    smi = "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1"
    ag, err = get_agents_for_engine(smi, top_k=3)
    print(json.dumps({"agents": ag, "error": err}, ensure_ascii=False))
