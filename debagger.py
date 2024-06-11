from .chython import smiles
lst1 = 'CCN1C[C@H](OC1=O)[C@H](O)[C@H](CC1CCCCC1)NC(=O)[C@H](CC1=CC=CS1)N=C(O)[C@@H](CC(=O)N1CCOCC1)CC1=CC=CC=C1'
lst2 = 'C([H])([H])([H])([H])'
smile = lst2

m = smiles(smile)
m.clean2d()
