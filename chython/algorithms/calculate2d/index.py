from .pikachu.drawing.drawing import Drawer
from .pikachu.smiles.smiles import read_smiles
def clean2d(smiles):
    structure=read_smiles(smiles)
    drawer = Drawer(structure)
    drawer._process_structure()
    drawn_atoms = drawer.structure.get_drawn_atoms()

    xy = []
    for atom in drawn_atoms:
        x, y = atom.draw.position.x, atom.draw.position.y
        xy.append([x, y])
    return xy

__all__ = ['clean2d']