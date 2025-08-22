import os, sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))
from enhanced_recommendation_engine import create_recommendation_engine

if __name__ == "__main__":
    e = create_recommendation_engine()
    rxn = 'Brc1ccccc1.NH2c1ccccc1>>Nc1ccccc1'
    recs = e.get_recommendations(rxn, 'C-N Coupling - Ullmann')
    ligs = recs.get('ligand_recommendations', [])[:5]
    print('TYPE', recs.get('reaction_type'))
    print('TOP_LIGANDS', [(r['ligand'], r.get('compatibility_score')) for r in ligs])
