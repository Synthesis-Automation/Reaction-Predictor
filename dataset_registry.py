"""
Dataset Registry
================

Centralized mapping from reaction types to dataset filenames and helpers to
resolve the correct on-disk path. Extend `DATASET_MAP` as new reaction types
and datasets are added.

Usage:
    from dataset_registry import resolve_dataset_path
    path = resolve_dataset_path(reaction_type, base_dir=__file__dir__)
"""

from __future__ import annotations

import os
from typing import Optional, Dict, List

# Map GUI reaction type strings to dataset basenames
# Add new entries as you create more datasets.
DATASET_MAP: Dict[str, str] = {
    # Buchwald (aliases)
    "C-N Coupling - Buchwald-Hartwig": "Buchwald-2021-2014.tsv",
    "Buchwald-Hartwig Amination": "Buchwald-2021-2014.tsv",

    # Ullmann (C-N and C-O variants share the same dataset file for now)
    "C-N Coupling - Ullmann": "Ullman-2020-2024.tsv",
    # Old label kept for back-compat in samples/docs
    "C-O Coupling - Ullmann Ether": "Ullman-2020-2024.tsv",
    # New reorganized label (metal in parentheses removed via normalization below)
    "C-O Coupling - Ullmann": "Ullman-2020-2024.tsv",
    "Ullmann Reaction": "Ullman-2020-2024.tsv",
    "Ullmann Ether Synthesis": "Ullman-2020-2024.tsv",

    # Amide formation / amidation
    "Amidation": "Amide-formation.tsv",
    "Amide formation": "Amide-formation.tsv",
}

# Fallback keyword routing for unexpected labels
KEYWORD_FALLBACKS: List[tuple[str, str]] = [
    ("ullmann", "Ullman-2020-2024.tsv"),
    ("buchwald", "Buchwald-2021-2014.tsv"),
    ("amid", "Amide-formation.tsv"),
    # Treat generic cross-coupling and Chan-Lam as falling back to the
    # Buchwald dataset for similarity browsing when nothing else is available
    ("cross-coupling", "Buchwald-2021-2014.tsv"),
    ("chan-lam", "Buchwald-2021-2014.tsv"),
]


def _candidate_paths(base_dir: str, basename: str) -> List[str]:
    return [
        os.path.join(base_dir, 'data', 'reaction_dataset', basename),
        os.path.join(base_dir, 'data', basename),
    ]


def resolve_dataset_path(reaction_type: Optional[str], base_dir: Optional[str] = None) -> Optional[str]:
    """Resolve dataset path for a given reaction type.

    - Checks explicit mapping in DATASET_MAP.
    - Falls back to keyword routing in KEYWORD_FALLBACKS.
    - Searches common locations: data/reaction_dataset/, data/.
    Returns absolute path if found; otherwise None.
    """
    if not base_dir:
        base_dir = os.path.dirname(__file__)

    label = (reaction_type or "").strip()

    # Normalize labels like "C-N Coupling - Buchwald-Hartwig (Pd)" ->
    # "C-N Coupling - Buchwald-Hartwig" to tolerate metal tags in UI labels
    try:
        if label.endswith(')') and ' (' in label:
            label_base = label[: label.rfind(' (')]
        else:
            label_base = label
    except Exception:
        label_base = label

    # Direct match
    if label_base in DATASET_MAP:
        basename = DATASET_MAP[label_base]
        for p in _candidate_paths(base_dir, basename):
            if os.path.exists(p):
                return p

    # Keyword fallback (case-insensitive contains)
    low = label.lower()
    for key, basename in KEYWORD_FALLBACKS:
        if key in low:
            for p in _candidate_paths(base_dir, basename):
                if os.path.exists(p):
                    return p

    return None


def list_available_datasets(base_dir: Optional[str] = None) -> Dict[str, str]:
    """Return a mapping of reaction type -> resolved dataset path for those present on disk."""
    if not base_dir:
        base_dir = os.path.dirname(__file__)
    found: Dict[str, str] = {}
    for k, basename in DATASET_MAP.items():
        for p in _candidate_paths(base_dir, basename):
            if os.path.exists(p):
                found[k] = p
                break
    return found
