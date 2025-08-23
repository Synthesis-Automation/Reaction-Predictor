from __future__ import annotations

import re
import unicodedata
from typing import Dict, List, Tuple

# Simple in-module synonym maps; can be upgraded to YAML files later
BASE_SYNONYMS = {
    "k2co3": ["potassium carbonate", "k_2co_3", "potassium carbonate (k2co3)", "potassiumcarbonate"],
    "cs2co3": ["cesium carbonate", "cesium carbonate (cs2co3)", "caesium carbonate"],
    "na2co3": ["sodium carbonate", "sodium carbonate (na2co3)"],
    "kotbu": ["potassium tert-butoxide", "potassium t-butoxide", "t-buok", "buok", "tbuok"],
    "naotbu": ["sodium tert-butoxide", "sodium t-butoxide"],
    "k3po4": ["tripotassium phosphate", "potassium phosphate", "kpo4", "tripotassiumphosphate"],
    "koac": ["potassium acetate", "k acetate"],
    "naome": ["sodium methoxide", "na ome"],
    "koh": ["potassium hydroxide"],
    "naoh": ["sodium hydroxide"],
    "tea": ["triethylamine", "et3n"],
    "dipea": ["diisopropylethylamine", "hunig's base"],
}

SOLVENT_SYNONYMS = {
    "dmf": ["n,n-dimethylformamide", "dimethylformamide"],
    "dmso": ["dimethyl sulfoxide", "dms o", "dm s o", "dimethylsulfoxide", "dimethyl sulphoxide"],
    "toluene": ["phme"],
    "meoh": ["methanol"],
    "etoh": ["ethanol"],
    "acn": ["acetonitrile", "me cn", "mecn", "me-c n"],
    "nmp": ["n-methyl-2-pyrrolidone", "n-methylpyrrolidone"],
}

LIGAND_SYNONYMS = {
    "l-proline": ["l proline", "proline"],
    "phen": ["1,10-phenanthroline", "o-phenanthroline", "phenanthroline"],
    "bipy": ["2,2'-bipyridine", "bipyridine"],
    "pPh3": ["triphenylphosphine", "p(ph)3", "p ph3"],
    "xphos": ["x-phos"],
}

METAL_SYNONYMS = {
    "cu": ["copper", "cui", "cu(i)", "cu(ii)", "cuprous iodide", "cupric acetate"],
    "pd": ["palladium"],
    "ni": ["nickel"],
}

_RE_NON_ALNUM = re.compile(r"[^a-z0-9]+")


def canonicalize(token: str) -> str:
    """Lowercase, strip punctuation/whitespace, normalize unicode.
    Also collapse runs of non-alphanum to a single space and then remove spaces.
    """
    if token is None:
        return ""
    s = unicodedata.normalize("NFKC", str(token)).casefold().strip()
    s = s.replace("µ", "u")
    s = _RE_NON_ALNUM.sub(" ", s)
    s = "".join(s.split())
    return s


def synonym_lookup(tok: str, table: Dict[str, List[str]]) -> str:
    ctok = canonicalize(tok)
    for k, vals in table.items():
        if ctok == k:
            return k
        for v in vals:
            if canonicalize(v) == ctok:
                return k
    return ctok


def normalize_mixture(text: str) -> List[str]:
    """Split mixtures like "DMF/THF", "toluene:MeCN", or comma-separated lists.
    Returns a list of canonicalized tokens.
    """
    if not text:
        return []
    parts = re.split(r"[\/:,+]|;|\\|\s+and\s+", text)
    return [canonicalize(p) for p in parts if p.strip()]


def parse_numeric(text: str) -> Tuple[float | None, str | None]:
    """Extract a float and a unit suffix if present. Returns (value, unit).
    Examples: "110 C" -> (110.0, "c"), "12 h" -> (12.0, "h"), "5 mol%" -> (5.0, "mol%")
    """
    if text is None:
        return None, None
    s = unicodedata.normalize("NFKC", str(text)).strip()
    # tolerate commas and localized decimals
    s = s.replace(",", ".")
    # quick lower for unit detection
    slow = s.casefold()
    m = re.search(r"([\-+]?[0-9]*\.?[0-9]+)", slow)
    if not m:
        return None, None
    val = float(m.group(1))
    # crude unit capture
    unit = None
    if "mol" in slow:
        unit = "mol%"
    elif "°c" in slow or slow.endswith(" c") or " deg c" in slow or "c" == slow:
        unit = "c"
    elif "h" in slow or "hr" in slow or slow.endswith(" hours"):
        unit = "h"
    elif "equiv" in slow or re.search(r"\beq\b", slow):
        unit = "equiv"
    return val, unit


def map_base(name: str) -> str:
    return synonym_lookup(name, BASE_SYNONYMS)


def map_solvent(name: str) -> str:
    # Fix common spacing issues first
    name = name.replace("DMS O", "DMSO")
    return synonym_lookup(name, SOLVENT_SYNONYMS)


def map_ligand(name: str) -> str:
    return synonym_lookup(name, LIGAND_SYNONYMS)


def map_metal(name: str) -> str:
    return synonym_lookup(name, METAL_SYNONYMS)
