from .index import clean2d as cl2d
from random import random
from typing import TYPE_CHECKING, Union

try:
    from importlib.resources import files
except ImportError:  # python3.8
    from importlib_resources import files


if TYPE_CHECKING:
    from chython import ReactionContainer, MoleculeContainer


class Calculate2DMolecule:
    __slots__ = ()

    def clean2d(self: Union['MoleculeContainer', 'Calculate2DMolecule']):
        plane = {}
        entry = iter(sorted(self, key=lambda n: len(self._bonds[n])))
        for _ in range(min(5, len(self))):
            smiles, order = self.__clean2d_prepare(next(entry))
            xy = cl2d(smiles)
            if xy is not None:
                break

        shift_x, shift_y = xy[0]
        for n, (x, y) in zip(order, xy):
            plane[n] = (x - shift_x, shift_y - y)

        bonds = []
        for n, m, _ in self.bonds():
            xn, yn = plane[n]
            xm, ym = plane[m]
            bonds.append(((xm - xn) ** 2 + (ym - yn) ** 2) ** 0.5)
        if bonds:
            bond_reduce = sum(bonds) / len(bonds) / 0.825
        else:
            bond_reduce = 1.0

        self_plane = self._plane
        for n, (x, y) in plane.items():
            self_plane[n] = (x / bond_reduce, y / bond_reduce)

        if self.connected_components_count > 1:
            shift_x = 0.0
            for c in self.connected_components:
                shift_x = self._fix_plane_mean(shift_x, component=c) + 0.9
        self.__dict__.pop('__cached_method__repr_svg_', None)


    def _fix_plane_mean(self: 'MoleculeContainer', shift_x: float, shift_y=0., component=None) -> float:
        plane = self._plane
        if component is None:
            component = plane

        left_atom = min(component, key=lambda x: plane[x][0])
        right_atom = max(component, key=lambda x: plane[x][0])

        min_x = plane[left_atom][0] - shift_x
        if len(self._atoms[left_atom].atomic_symbol) == 2:
            min_x -= 0.2

        max_x = plane[right_atom][0] - min_x
        min_y = min(plane[x][1] for x in component)
        max_y = max(plane[x][1] for x in component)
        mean_y = (max_y + min_y) / 2 - shift_y
        for n in component:
            x, y = plane[n]
            plane[n] = (x - min_x, y - mean_y)

        if -0.18 <= plane[right_atom][1] <= 0.18:
            factor = self._hydrogens[right_atom]
            if factor == 1:
                max_x += 0.15
            elif factor:
                max_x += 0.25
        return max_x

    def _fix_plane_min(self: 'MoleculeContainer', shift_x: float, shift_y=0., component=None) -> float:
        plane = self._plane
        if component is None:
            component = plane

        right_atom = max(component, key=lambda x: plane[x][0])
        min_x = min(plane[x][0] for x in component) - shift_x
        max_x = plane[right_atom][0] - min_x
        min_y = min(plane[x][1] for x in component) - shift_y

        for n in component:
            x, y = plane[n]
            plane[n] = (x - min_x, y - min_y)

        if shift_y - 0.18 <= plane[right_atom][1] <= shift_y + 0.18:
            factor = self._hydrogens[right_atom]
            if factor == 1:
                max_x += 0.15
            elif factor:
                max_x += 0.25
        return max_x

    def __clean2d_prepare(self: 'MoleculeContainer', entry):
        hydrogens = self._hydrogens
        charges = self._charges
        allenes_stereo = self._allenes_stereo
        atoms_stereo = self._atoms_stereo
        self._charges = self._hydrogens = {n: 0 for n in hydrogens}
        self._atoms_stereo = self._allenes_stereo = {}
        w = {n: random() for n in hydrogens}
        w[entry] = -1
        try:
            smiles, order = self._smiles(w.__getitem__, random=True, _return_order=True)
        finally:
            self._hydrogens = hydrogens
            self._charges = charges
            self._allenes_stereo = allenes_stereo
            self._atoms_stereo = atoms_stereo
        return ''.join(smiles).replace('~', '-'), order

class Calculate2DReaction:
    __slots__ = ()

    def clean2d(self: 'ReactionContainer'):
        for m in self.molecules():
            m.clean2d()
        self.fix_positions()

    def fix_positions(self: 'ReactionContainer'):
        shift_x = 0
        reactants = self.reactants
        amount = len(reactants) - 1
        signs = []
        for m in reactants:
            max_x = m._fix_plane_mean(shift_x)
            if amount:
                max_x += .2
                signs.append(max_x)
                amount -= 1
            shift_x = max_x + 1
        arrow_min = shift_x

        if self.reagents:
            shift_x += .4
            for m in self.reagents:
                max_x = m._fix_plane_min(shift_x, .5)
                shift_x = max_x + 1
            shift_x += .4
            if shift_x - arrow_min < 3:
                shift_x = arrow_min + 3
        else:
            shift_x += 3
        arrow_max = shift_x - 1

        products = self.products
        amount = len(products) - 1
        for m in products:
            max_x = m._fix_plane_mean(shift_x)
            if amount:
                max_x += .2
                signs.append(max_x)
                amount -= 1
            shift_x = max_x + 1
        self._arrow = (arrow_min, arrow_max)
        self._signs = tuple(signs)
        self.flush_cache()


__all__ = ['Calculate2DMolecule', 'Calculate2DReaction']