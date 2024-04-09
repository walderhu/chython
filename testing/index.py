import pikachu

def clean2d(smiles):
    # drawer = pikachu.drawing
    # drawer = DrawerBase({})
    # parsed = Parser.parse(smiles)

    # drawer.initDraw(parsed, 'light', False)
    # drawer.processGraph()
    # vertices = drawer.graph.vertices


    xy = list()
    for i in range(len(vertices)):
      position = vertices[i].position
      xy.append([position.x, position.y])

    return xy

__all__ = ['clean2d']