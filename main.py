from pikachu.drawing.drawing import Drawer
from pikachu.smiles.smiles import read_smiles

def clean2d(smiles):
    drawer = Drawer(
       structure=read_smiles(smiles)
    )
    drawer.process_structure()
  
    drawn_atoms = drawer.structure.get_drawn_atoms()
    coordinates = [(atom.draw.position.x, atom.draw.position.y) for atom in drawn_atoms]
    return coordinates
    
smiles = 'CC=CC'
clean2d(smiles)