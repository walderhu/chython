import unittest
from chython import smiles


class TestAddFunction(unittest.TestCase):

    def test(self):
        lst = [
            'CCCCCCCCCC(=O)N[C@@H](CC1=CNC2=CC=CC=C21)C(=O)N[C@@H](CC(=O)N)C(=O)N[C@@H](CC(=O)O)C(=O)N[C@H]3[C@H](OC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@H](NC(=O)CNC(=O)[C@@H](NC(=O)[C@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)CNC3=O)CCCN)CC(=O)O)C)CC(=O)O)CO)[C@H](C)CC(=O)O)CC(=O)C4=CC=CC=C4N)C',
            'CC=O',
            '[Na+].CN(C)c2ccc(/N=N/c1ccc(cc1)S([O-])(=O)=O)cc2',
            'CC=CC'#,
            # 'C1CCCCC1', # Циклогексан
            # 'C1CCCCC1C1CCCCC1', # Циклогексадиен
            # 'C1CCCCC1C1CCCCC1C1CCCCC1', # Циклогексатриен
            # 'c1ccccc1', # Бензол
            # 'c1ccccn1', # Пиридин
            # 'c1cccc2c1cccc2', # Фуран
            # 'C1=CC=CC=C1C=O', # Бензол с формальдегидной группой
                  ]
        for smile in lst:
            m = smiles(smile)
            m.clean2d()


if __name__ == '__main__':
    unittest.main()
