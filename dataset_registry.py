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
    "C-N Coupling - Buchwald-Hartwig": "buchwald_reactions.csv",
    "Buchwald-Hartwig Amination": "buchwald_reactions.csv",

    # Ullmann (C-N and C-O variants share the same dataset file for now)
    "C-N Coupling - Ullmann": "Ullman_2024_full.csv",
    "C-O Coupling - Ullmann Ether": "Ullman_2024_full.csv",
    "Ullmann Reaction": "Ullman_2024_full.csv",
    "Ullmann Ether Synthesis": "Ullman_2024_full.csv",
}

# Fallback keyword routing for unexpected labels
KEYWORD_FALLBACKS: List[tuple[str, str]] = [
    ("ullmann", "Ullman_2024_full.csv"),
    ("buchwald", "buchwald_reactions.csv"),
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

    # Direct match
    if label in DATASET_MAP:
        basename = DATASET_MAP[label]
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
