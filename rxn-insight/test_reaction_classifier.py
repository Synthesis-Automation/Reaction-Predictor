from rxn_insight.reaction import Reaction
for i in range(10):
    smile = input("Enter your reaction smiles: ")
    rxn = Reaction(smile)

    info = rxn.get_reaction_info()  # dict with class + name + features
    # print(info["CLASS"], info["NAME"])
    print(info)
