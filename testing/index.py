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

__all__ = ['clean2d']

"""создаем экземпляр класса Drawer (аналог DrawerBase), 
в который передаем параметр smiles через функцию read_smiles, 
которая  парсит строку на соответствие стандарту
и в случае успеха возвращает объект класса Structure.
второй аргумент конструктора это класс Options, который по умолчанию 
устанавливает отрисовку цветом white.
process_structure является аналогом processGraph в js
позже получаем все атомы нашей молекулы и создаем массив координат их вершин
"""


