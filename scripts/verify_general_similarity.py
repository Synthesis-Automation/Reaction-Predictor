#!/usr/bin/env python3
import json

from enhanced_recommendation_engine import create_recommendation_engine


def main():
    engine = create_recommendation_engine()
    rxn = 'Brc1ccccc1.NH2c1ccccc1>>Nc1ccccc1'

    # Targeted Ullmann path
    recs = engine.get_recommendations(rxn, 'C-N Coupling - Ullmann')
    print('Ullmann TYPE:', recs.get('reaction_type'))
    ligs = recs.get('ligand_recommendations', [])
    print('Ullmann TOP LIG:', [(r.get('ligand'), r.get('compatibility_score')) for r in ligs[:5]])

    # Auto-detect + general similarity fallback
    recs2 = engine.get_recommendations(rxn, 'Auto-detect')
    print('Auto TYPE:', recs2.get('reaction_type'))
    gen = recs2.get('general_recommendations')
    if gen:
        print('General Ligands:', gen.get('ligand_recommendations'))
        print('General Solvents:', gen.get('solvent_recommendations'))
        print('General Bases:', gen.get('base_recommendations'))
        print('Similarity Hits (top):', gen.get('top_hits'))
    else:
        print('General recommendation not available (RDKit missing or no similar reactions found).')


if __name__ == '__main__':
    main()
