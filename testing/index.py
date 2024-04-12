from pikachu.drawing.drawing import Drawer
from pikachu.smiles.smiles import read_smiles

def clean2d(smiles):
    drawer = Drawer(
       structure=read_smiles(smiles_string=smiles)
    )
    drawer.process_structure()
  
    vertices = drawer.structure.graph
    
    # drawer = DrawerBase({})
    # parsed = Parser.parse(smiles)
    # drawer.initDraw(parsed, 'light', False)
    # drawer.processGraph()

    vertices = drawer.graph.vertices
    xy = list()
    for i in range(len(vertices)):
      position = vertices[i].position
      xy.append([position.x, position.y])
# OR
    # xy = [[point.position.x, point.position.y] for point in vertices]


    return xy

__all__ = ['clean2d']

"""создаем экземпляр класса Drawer (аналог DrawerBase), 
в который передаем параметр smiles через функцию read_smiles, 
которая  парсит строку на соответствие стандарту
и в случае успеха возвращает объект класса Structure.
второй аргумент конструктора это класс Options, который по умолчанию 
устанавливает отрисовку цветом white.
process_structure является аналогом processGraph в js"""


