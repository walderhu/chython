# from chython import smiles

# m = smiles('CC=CC')
# m.clean2d()
# m



# import pikachu
# import pikachu.general
# import pikachu.parsers
# import pikachu.parsers.coconut

# # pikachu.general.draw_smiles('CC=O')
# pikachu.parsers.


import pikachu 

def clean2d(smiles):
    m = pikachu.Molecule
    mol = pk.Molecule.from_smiles(smiles)
    mol.clean2d()
    xy = [[atom.x, atom.y] for atom in mol.atoms]
    return xy