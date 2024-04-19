from pikachu.drawing.drawing import Drawer
from pikachu.smiles.smiles import read_smiles

def clean2d(smiles):
    drawer = Drawer(
       structure=read_smiles(smiles)
    )
    drawer.process_structure()
  
    drawn_atoms = drawer.structure.get_drawn_atoms()
    xy = [[atom.draw.position.x, atom.draw.position.y] for atom in drawn_atoms]
    return xy
    
smiles = 'CC=CC'
print(*clean2d(smiles), sep='\n')