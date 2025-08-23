from enhanced_recommendation_engine import create_recommendation_engine

e = create_recommendation_engine()
recs = e.get_recommendations('Brc1ccccc1.NH2c1ccccc1>>Nc1ccccc1','C-N Coupling - Ullmann')

ligs = recs.get('ligand_recommendations', [])[:5]
sols = recs.get('solvent_recommendations', [])[:5]
bases = recs.get('base_recommendations', [])[:5]

print('TYPE', recs.get('reaction_type'))
print('LIG', [(r.get('ligand'), r.get('compatibility_score')) for r in ligs])
print('SOL', [(r.get('solvent'), r.get('compatibility_score')) for r in sols])
print('BAS', [(r.get('base'), r.get('compatibility_score')) for r in bases])
