import os
from analytics.adapters import adapt_dataset_for_type


def test_adapt_ullmann_loads(tmp_path):
    # Create a tiny CSV with required columns
    p = tmp_path / "mini_ullmann.csv"
    p.write_text(
        "ReactionType,CoreGeneric,Ligand,ReagentRaw,Solvent,Temperature_C,Time_h,Yield_%,Reference\n"
        "Ullmann,CuI,1,10-Phenanthroline,K2CO3,DMSO,120,12,78,Ref1\n",
        encoding="utf-8",
    )
    rows = adapt_dataset_for_type("Ullmann", str(p))
    assert rows and rows[0].reaction_type
