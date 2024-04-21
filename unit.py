import unittest
from chython import smiles


class TestAddFunction(unittest.TestCase):

    def test(self):
        smile = 'CCCCCCCCCC(=O)N[C@@H](CC1=CNC2=CC=CC=C21)C(=O)N[C@@H](CC(=O)N)C(=O)N[C@@H](CC(=O)O)C(=O)N[C@H]3[C@H](OC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@H](NC(=O)CNC(=O)[C@@H](NC(=O)[C@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)CNC3=O)CCCN)CC(=O)O)C)CC(=O)O)CO)[C@H](C)CC(=O)O)CC(=O)C4=CC=CC=C4N)C'
        m = smiles(smile)
        m.clean2d()

    def test2(self):
        smile = 'CC=O'
        m = smiles(smile)
        m.clean2d()

if __name__ == '__main__':
    unittest.main()