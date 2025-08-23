from analytics.normalization import canonicalize, parse_numeric, map_solvent, map_base


def test_canonicalize_basic():
    assert canonicalize(" DMS O ") == canonicalize("dmso")
    assert canonicalize("K 2 CO 3") == canonicalize("k2co3")


def test_parse_numeric_units():
    assert parse_numeric("110 C")[1] == "c"
    assert parse_numeric("12 h")[1] == "h"
    assert parse_numeric("5 mol%")[1] == "mol%"
    assert parse_numeric("2.0 equiv")[1] == "equiv"


def test_solvent_mapping_dmso():
    assert map_solvent("DMS O") == map_solvent("dmso")


def test_base_mapping_aliases():
    assert map_base("potassium carbonate") == "k2co3"
    assert map_base("K3PO4") == "k3po4"
