import sys

from ..chem.bond_properties import BOND_PROPERTIES
from ..errors import StructureError, KekulisationError
from ..chem.atom import Atom
from ..chem.bond import Bond
from ..chem.kekulisation import Match
from ..chem.rings.ring_identification import is_aromatic
from ..chem.rings.find_cycles import Cycles
from ..chem.aromatic_system import AromaticSystem

sys.setrecursionlimit(1000000)


class Structure:
    """
    Class to store a molecular structure

    Attributes
    ----------
    graph: dict of {atom: [atom, ->], ->}
        with each atom (Atom object) pointing to a list of atoms
        covalently bonded to that atom
    bonds: dict of {bond nr: bond, ->}
        with bond nr (int) an arbitrarily chosen and unique number for a bond,
        and bond (Bond object) the bond associated with that number
    bond_lookup: dict of {atom: {atom: bond, ->}, ->}
        with each pair of atoms (Atom object) pointing to the bond (Bond)
        between those atoms
    """

    def __init__(self, graph=None, bonds=None, bond_lookup=None):

        if graph:
            self.graph = graph
            self.set_atom_neighbours()
            self.set_atoms()
        else:
            self.graph = {}
            self.atoms = {}
        if bonds:
            self.bonds = bonds
            self.make_bond_lookup()
        else:
            self.bonds = {}

        if bond_lookup:
            self.bond_lookup = bond_lookup
        else:
            self.make_bond_lookup()

        self.annotations = set()

        self.aromatic_cycles = []
        self.aromatic_systems = []

    def get_atoms(self):
        atoms = []
        for atom in self.graph:
            atoms.append(atom)

        return atoms

    def bond_exists(self, atom_1, atom_2):
        if atom_1 in self.bond_lookup:
            if atom_2 in self.bond_lookup[atom_1]:
                return True
            return False
        return False

    def get_subtree(self, atom, parent_atom):
        pass

    def set_atoms(self):
        self.atoms = {}
        for atom in self.graph:
            self.atoms[atom.nr] = atom
            
    def deepcopy(self):

        new_graph = {}
        new_bonds = {}
        new_atoms = {}

        for atom_nr, atom in self.atoms.items():
            new_atoms[atom_nr] = atom.copy()
        
        for atom_1, atoms in self.graph.items():
            new_atom_1 = new_atoms[atom_1.nr]
            new_graph[new_atom_1] = []
            for atom_2 in atoms:
                new_atom_2 = new_atoms[atom_2.nr]
                new_graph[new_atom_1].append(new_atom_2)
                
        for bond_nr, bond in self.bonds.items():
            new_atom_1 = new_atoms[bond.atom_1.nr]
            new_atom_2 = new_atoms[bond.atom_2.nr]
            new_bond = Bond(new_atom_1, new_atom_2, bond.type, bond.nr)
            new_bond.chiral = bond.chiral
            new_bond.chiral_symbol = bond.chiral_symbol

            for atom_1, atoms_and_chirality in bond.chiral_dict.items():
                new_1 = None

                if isinstance(atom_1, Atom):
                    new_1 = new_atoms[atom_1.nr]
                else:
                    pass

                assert new_1

                new_bond.chiral_dict[new_1] = {}

                for atom_2, chirality in atoms_and_chirality.items():
                    if isinstance(atom_2, Atom):
                        new_2 = new_atoms[atom_2.nr]
                        new_bond.chiral_dict[new_1][new_2] = chirality
                    else:
                        pass

            new_bonds[bond_nr] = new_bond
            if new_bond not in new_atom_1.bonds:
                new_atom_1.bonds.append(new_bond)
            if new_bond not in new_atom_2.bonds:
                new_atom_2.bonds.append(new_bond)
                
        new_structure = Structure(new_graph, new_bonds)

        new_structure.add_shells_non_hydrogens()
        new_structure.add_shells()

        new_structure.form_pi_bonds()
        new_structure.hybridise_atoms()
        new_structure.promote_pi_bonds()

        new_structure.find_cycles()
        new_structure.aromatic_cycles = new_structure.find_aromatic_cycles()
        new_structure.aromatic_systems = new_structure.find_aromatic_systems()

        new_structure.form_sigma_bonds()
        new_structure.drop_electrons()
        new_structure.make_lone_pairs()

        new_structure.set_atom_neighbours()
        new_structure.make_bond_lookup()
        new_structure.set_connectivities()

        for annotation, default in self.annotations:
            new_structure.annotations.add((annotation, default))

        return new_structure

    def copy(self):
        new_graph = {}
        new_bonds = {}

        for atom_1, atoms in self.graph.items():
            new_graph[atom_1] = []
            for atom_2 in atoms:
                new_graph[atom_1].append(self.atoms[atom_2.nr])

        for bond_nr, bond in self.bonds.items():
            new_bonds[bond_nr] = bond

        new_structure = Structure(new_graph, new_bonds)
        new_structure.cycles = self.cycles

        return new_structure

    def refresh_structure(self, find_cycles=False):
        new_graph = {}
        self.set_atoms()

        for atom_1, atoms in self.graph.items():
            new_graph[atom_1] = []
            for atom_2 in atoms:
                new_graph[atom_1].append(self.atoms[atom_2.nr])

        for bond_nr, bond in self.bonds.items():
            bond.atom_1 = self.atoms[bond.atom_1.nr]
            bond.atom_2 = self.atoms[bond.atom_2.nr]
            if bond not in bond.atom_1.bonds:
                bond.atom_1.bonds.append(bond)
            if bond not in bond.atom_2.bonds:
                bond.atom_2.bonds.append(bond)

        self.graph = new_graph

        self.make_bond_lookup()
        self.set_atom_neighbours()
        self.set_connectivities()
        for bond in self.bonds.values():
            bond.set_bond_summary()

        if find_cycles:
            self.find_cycles()

    def get_next_in_ring(self, ring, current_atom, previous_atom):
        neighbours = self.graph[current_atom]

        for neighbour in neighbours:
            for member in ring.members:
                if neighbour == member:
                    if previous_atom != neighbour:
                        return neighbour

        return None

    def get_drawn_atoms(self):
        atoms = []
        for atom in self.graph:
            if atom.draw.is_drawn:
                atoms.append(atom)
        return atoms

    def set_double_bond_chirality(self):
        for bond_nr, bond in self.bonds.items():
            # iterate over all double bonds
            if bond.type == 'double' or bond.type == 'triple':
                if len(bond.atom_1.bonds) + len(bond.atom_1.lone_pairs) == 3 and \
                        len(bond.atom_2.bonds) + len(bond.atom_2.lone_pairs) == 3 and \
                        len(bond.atom_1.lone_pairs) < 2 and len(bond.atom_2.lone_pairs) < 2:
                    # define atoms adjacent to the atoms involved in the double bond
                    # also keep track of the chiral symbol that defines these bonds

                    atom_1_1 = None
                    atom_1_2 = None
                    chiral_1_1 = None
                    chiral_1_2 = None

                    atom_2_1 = None
                    atom_2_2 = None
                    chiral_2_1 = None
                    chiral_2_2 = None

                    if bond.atom_1.lone_pairs:
                        atom_1_1 = bond.atom_1.lone_pairs[0]
                    if bond.atom_2.lone_pairs:
                        atom_2_1 = bond.atom_2.lone_pairs[0]
                    # Check bonds adjacent to the first atom
                    for bond_1 in bond.atom_1.bonds:

                        if bond_1.type == 'single':
                            # Looks at the bonds between the atom adjacent to the stereobond and its neighbours
                            if bond.atom_1 == bond_1.atom_1:
                                if bond_1.chiral_symbol == '/':
                                    direction = 'up'
                                elif bond_1.chiral_symbol == '\\':
                                    direction = 'down'
                                else:
                                    direction = None
                                # First time it runs through this it will define atom_1_1
                                if not atom_1_1:
                                    atom_1_1 = bond_1.atom_2
                                    chiral_1_1 = direction
                                # Second time it runs through this, it will define atom_1_2
                                else:
                                    atom_1_2 = bond_1.atom_2
                                    chiral_1_2 = direction

                            elif bond.atom_1 == bond_1.atom_2:

                                if bond_1.chiral_symbol == '/':
                                    direction = 'down'
                                elif bond_1.chiral_symbol == '\\':
                                    direction = 'up'
                                else:
                                    direction = None

                                if not atom_1_1:
                                    atom_1_1 = bond_1.atom_1
                                    chiral_1_1 = direction
                                else:
                                    atom_1_2 = bond_1.atom_1
                                    chiral_1_2 = direction

                    for bond_2 in bond.atom_2.bonds:

                        if bond_2.type == 'single':
                            if bond.atom_2 == bond_2.atom_1:
                                if bond_2.chiral_symbol == '/':
                                    direction = 'up'
                                elif bond_2.chiral_symbol == '\\':
                                    direction = 'down'
                                else:
                                    direction = None

                                if not atom_2_1:
                                    atom_2_1 = bond_2.atom_2
                                    chiral_2_1 = direction
                                else:
                                    atom_2_2 = bond_2.atom_2
                                    chiral_2_2 = direction

                            elif bond.atom_2 == bond_2.atom_2:
                                if bond_2.chiral_symbol == '/':
                                    direction = 'down'
                                elif bond_2.chiral_symbol == '\\':
                                    direction = 'up'
                                else:
                                    direction = None

                                if not atom_2_1:
                                    atom_2_1 = bond_2.atom_1
                                    chiral_2_1 = direction
                                else:
                                    atom_2_2 = bond_2.atom_1
                                    chiral_2_2 = direction

                    chiral_1 = False
                    chiral_2 = False

                    if chiral_1_1 or chiral_1_2:
                        chiral_1 = True

                    if chiral_2_1 or chiral_2_2:
                        chiral_2 = True

                    if chiral_1 and chiral_2:
                        if chiral_1_1 == chiral_1_2:
                            raise StructureError('chiral double bond')
                        if chiral_2_2 == chiral_2_1:
                            raise StructureError('chiral double bond')

                        if chiral_1_1:
                            first_atom = atom_1_1
                            first_other_atom = atom_1_2
                            first_chiral_symbol = chiral_1_1

                            if not chiral_1_2 and isinstance(atom_1_2, Atom):
                                # Make sure where chiral symbols are not defined, they are added

                                if (atom_1_1.nr > bond.atom_1.nr and atom_1_2.nr > bond.atom_1.nr) or \
                                        (atom_1_1.nr < bond.atom_1.nr and atom_1_2.nr < bond.atom_1.nr):
                                    if self.bond_lookup[bond.atom_1][atom_1_1].chiral_symbol == '/':
                                        self.bond_lookup[bond.atom_1][atom_1_2].chiral_symbol = '\\'
                                    else:
                                        self.bond_lookup[bond.atom_1][atom_1_2].chiral_symbol = '/'

                                else:
                                    if self.bond_lookup[bond.atom_1][atom_1_1].chiral_symbol == '/':
                                        self.bond_lookup[bond.atom_1][atom_1_2].chiral_symbol = '/'
                                    else:
                                        self.bond_lookup[bond.atom_1][atom_1_2].chiral_symbol = '\\'

                        else:
                            first_atom = atom_1_2
                            first_other_atom = atom_1_1
                            first_chiral_symbol = chiral_1_2

                            # Make sure where chiral symbols are not defined, they are added

                            if isinstance(atom_1_1, Atom):

                                if (atom_1_1.nr > bond.atom_1.nr and atom_1_2.nr > bond.atom_1.nr) or \
                                   (atom_1_1.nr < bond.atom_1.nr and atom_1_2.nr < bond.atom_1.nr):
                                    if self.bond_lookup[bond.atom_1][atom_1_2].chiral_symbol == '/':
                                        self.bond_lookup[bond.atom_1][atom_1_1].chiral_symbol = '\\'
                                    else:
                                        self.bond_lookup[bond.atom_1][atom_1_1].chiral_symbol = '/'

                                else:
                                    if self.bond_lookup[bond.atom_1][atom_1_2].chiral_symbol == '/':
                                        self.bond_lookup[bond.atom_1][atom_1_1].chiral_symbol = '/'
                                    else:
                                        self.bond_lookup[bond.atom_1][atom_1_1].chiral_symbol = '\\'

                        if chiral_2_1:
                            second_atom = atom_2_1
                            second_other_atom = atom_2_2
                            second_chiral_symbol = chiral_2_1

                            if not chiral_2_2 and isinstance(atom_2_2, Atom):
                                # Make sure where chiral symbols are not defined, they are added
                                if (atom_2_1.nr > bond.atom_2.nr and atom_2_2.nr > bond.atom_2.nr) or \
                                        (atom_2_1.nr < bond.atom_2.nr and atom_2_2.nr < bond.atom_2.nr):
                                    if self.bond_lookup[bond.atom_2][atom_2_1].chiral_symbol == '/':
                                        self.bond_lookup[bond.atom_2][atom_2_2].chiral_symbol = '\\'
                                    else:
                                        self.bond_lookup[bond.atom_2][atom_2_2].chiral_symbol = '/'

                                else:
                                    if self.bond_lookup[bond.atom_2][atom_2_1].chiral_symbol == '/':
                                        self.bond_lookup[bond.atom_2][atom_2_2].chiral_symbol = '/'
                                    else:
                                        self.bond_lookup[bond.atom_2][atom_2_2].chiral_symbol = '\\'

                        else:
                            second_atom = atom_2_2
                            second_other_atom = atom_2_1
                            second_chiral_symbol = chiral_2_2
                            # Make sure where chiral symbols are not defined, they are added
                            if isinstance(atom_2_1, Atom):
                                if (atom_2_1.nr > bond.atom_2.nr and atom_2_2.nr > bond.atom_2.nr) or \
                                        (atom_2_1.nr < bond.atom_2.nr and atom_2_2.nr < bond.atom_2.nr):
                                    if self.bond_lookup[bond.atom_2][atom_2_2].chiral_symbol == '/':
                                        self.bond_lookup[bond.atom_2][atom_2_1].chiral_symbol = '\\'
                                    else:
                                        self.bond_lookup[bond.atom_2][atom_2_1].chiral_symbol = '/'

                                else:
                                    if self.bond_lookup[bond.atom_2][atom_2_2].chiral_symbol == '/':
                                        self.bond_lookup[bond.atom_2][atom_2_1].chiral_symbol = '/'
                                    else:
                                        self.bond_lookup[bond.atom_2][atom_2_1].chiral_symbol = '\\'

                        if isinstance(first_atom, Atom):
                            bond.chiral_dict[first_atom] = {}
                        if isinstance(first_other_atom, Atom):
                            bond.chiral_dict[first_other_atom] = {}
                        if isinstance(second_atom, Atom):
                            bond.chiral_dict[second_atom] = {}
                        if isinstance(second_other_atom, Atom):
                            bond.chiral_dict[second_other_atom] = {}

                        if first_chiral_symbol == second_chiral_symbol:
                            if isinstance(first_atom, Atom) and isinstance(second_atom, Atom):
                                bond.chiral_dict[first_atom][second_atom] = 'cis'
                                bond.chiral_dict[second_atom][first_atom] = 'cis'

                            if isinstance(first_other_atom, Atom) and isinstance(second_other_atom, Atom):
                                bond.chiral_dict[first_other_atom][second_other_atom] = 'cis'
                                bond.chiral_dict[second_other_atom][first_other_atom] = 'cis'

                            if isinstance(first_atom, Atom) and isinstance(second_other_atom, Atom):
                                bond.chiral_dict[first_atom][second_other_atom] = 'trans'
                                bond.chiral_dict[second_other_atom][first_atom] = 'trans'

                            if isinstance(first_other_atom, Atom) and isinstance(second_atom, Atom):
                                bond.chiral_dict[first_other_atom][second_atom] = 'trans'
                                bond.chiral_dict[second_atom][first_other_atom] = 'trans'

                        else:
                            if isinstance(first_atom, Atom) and isinstance(second_atom, Atom):
                                bond.chiral_dict[first_atom][second_atom] = 'trans'
                                bond.chiral_dict[second_atom][first_atom] = 'trans'

                            if isinstance(first_other_atom, Atom) and isinstance(second_other_atom, Atom):
                                bond.chiral_dict[first_other_atom][second_other_atom] = 'trans'
                                bond.chiral_dict[second_other_atom][first_other_atom] = 'trans'

                            if isinstance(first_atom, Atom) and isinstance(second_other_atom, Atom):
                                bond.chiral_dict[first_atom][second_other_atom] = 'cis'
                                bond.chiral_dict[second_other_atom][first_atom] = 'cis'

                            if isinstance(first_other_atom, Atom) and isinstance(second_atom, Atom):
                                bond.chiral_dict[first_other_atom][second_atom] = 'cis'
                                bond.chiral_dict[second_atom][first_other_atom] = 'cis'

                        bond.chiral = True

    def make_bond_lookup(self):
        """
        Update bond_lookup from current collection of bonds
        """
        self.bond_lookup = {}
        for bond in self.bonds:
            atom_1 = self.bonds[bond].atom_1
            atom_2 = self.bonds[bond].atom_2
            if atom_1 not in self.bond_lookup:
                self.bond_lookup[atom_1] = {}
            if atom_2 not in self.bond_lookup:
                self.bond_lookup[atom_2] = {}
            self.bond_lookup[atom_1][atom_2] = self.bonds[bond]
            self.bond_lookup[atom_2][atom_1] = self.bonds[bond]

    def find_cycles(self):
        """
        Find all cycles in a structure and store them
        """
        self.cycles = Cycles(self)
        self.sssr = self.cycles.find_sssr()

    @staticmethod
    def promote_lone_pairs_in_aromatic_cycles(cycles):
        for cycle in cycles:
            for atom in cycle:
                if atom.hybridisation == 'sp3':
                    atom._promote_lone_pair_to_p_orbital()
                    if atom.type == 'N':
                        atom.pyrrole = True
                    elif atom.type == 'S':
                        atom.thiophene = True
                    elif atom.type == 'O':
                        atom.furan = True


    def find_aromatic_cycles(self):
        """
        Returns cycles that are aromatic

        Output
        ------
        aromatic_cycles: list of [[Atom, ->], ->], with each list of
            atoms the atoms that comprise an aromatic cycle

        """
        aromatic_cycles = set()
        previous_nr_aromatic_cycles = -1
        current_nr_aromatic_cycles = 0

        while previous_nr_aromatic_cycles != current_nr_aromatic_cycles:
            previous_nr_aromatic_cycles = current_nr_aromatic_cycles
            for cycle in self.sssr:
                if tuple(cycle) not in aromatic_cycles and is_aromatic(cycle):
                    self.make_cycle_aromatic(cycle)
                    aromatic_cycles.add(tuple(cycle))

            current_nr_aromatic_cycles = len(aromatic_cycles)

        self.promote_lone_pairs_in_aromatic_cycles(aromatic_cycles)

        return aromatic_cycles

    def find_aromatic_systems(self):
        """
        Return ring systems that are aromatic

        Output
        ------
        aromatic_ring_systems: list of [[Atom, ->], ->], with each list of
            atoms the atoms that comprise an aromatic ring system
        """
        previous_system_nr = -1
        current_system_nr = 1

        aromatic_systems = [list(aromatic_cycle) for aromatic_cycle in self.aromatic_cycles]

        while current_system_nr != previous_system_nr:
            previous_system_nr = current_system_nr
            indices_to_remove = None
            new_system = None
            for i, system_1 in enumerate(aromatic_systems):
                system_found = False
                for j, system_2 in enumerate(aromatic_systems):
                    if i != j:
                        if len(set(system_1).intersection(set(system_2))) >= 2:
                            indices_to_remove = [i, j]
                            new_system = list(set(system_1[:] + system_2[:]))
                            system_found = True
                            break
                if system_found:
                    break

            if new_system:
                indices_to_remove.sort(reverse=True)
                for index in indices_to_remove:
                    aromatic_systems.pop(index)

                aromatic_systems.append(new_system)
                
            current_system_nr = len(aromatic_systems)

        aromatic_ring_systems = []

        for i, aromatic_system in enumerate(aromatic_systems):
            system = AromaticSystem(i, aromatic_system)
            aromatic_ring_systems.append(system)

        return aromatic_ring_systems

    def find_double_bond_sequences(self):
        double_bond_fragments = []
        for bond in self.bonds.values():
            if bond.type == 'single':
                stereobond_1 = None
                stereobond_2 = None
                for bond_1 in bond.atom_1.bonds:
                    if bond_1.chiral:
                        stereobond_1 = bond_1
                        break

                for bond_2 in bond.atom_2.bonds:
                    if bond_2.chiral:
                        stereobond_2 = bond_2
                        break

                if stereobond_1 and stereobond_2:
                    double_bond_fragments.append([stereobond_1, stereobond_2])

        previous_fragment_nr = -1
        fragment_nr = len(double_bond_fragments)

        while fragment_nr != previous_fragment_nr:
            previous_fragment_nr = fragment_nr

            indices_to_remove = None
            new_fragment = None

            for i, fragment_1 in enumerate(double_bond_fragments):
                found_fragment = False
                for j, fragment_2 in enumerate(double_bond_fragments):
                    if i != j:
                        if fragment_1[-1] == fragment_2[0]:
                            new_fragment = fragment_1[:] + fragment_2[1:]
                        elif fragment_1[-1] == fragment_2[-1]:
                            new_fragment = fragment_1[:] + list(reversed(fragment_2[:-1]))
                        elif fragment_1[0] == fragment_2[0]:
                            new_fragment = list(reversed(fragment_2[1:])) + fragment_1[:]
                        elif fragment_1[0] == fragment_2[-1]:
                            new_fragment = fragment_2[:-1] + fragment_1[:]

                        if new_fragment:
                            indices_to_remove = [i, j]
                            found_fragment = True
                            break
                if found_fragment:
                    break

            if indices_to_remove:
                assert new_fragment
                # sort indices from large to small to make sure you first remove the last index
                indices_to_remove.sort(reverse=True)
                double_bond_fragments.pop(indices_to_remove[0])
                double_bond_fragments.pop(indices_to_remove[1])
                double_bond_fragments.append(new_fragment)
                fragment_nr = len(double_bond_fragments)

        return double_bond_fragments

    @staticmethod
    def make_cycle_aromatic(cycle):
        for atom_1 in cycle:
            atom_1.aromatic = True
            for atom_2 in cycle:
                if atom_1 != atom_2:
                    bond = atom_1.get_bond(atom_2)
                    if bond and bond.type != 'aromatic':
                        bond.make_aromatic()

    def refine_structure(self):
        """
        """

        self.add_shells()
        self.add_hydrogens()
        self.add_shells()

        self.sort_by_nr()

        self.form_pi_bonds()
        self.hybridise_atoms()
        self.promote_pi_bonds()

        self.find_cycles()
        self.aromatic_cycles = self.find_aromatic_cycles()
        self.aromatic_systems = self.find_aromatic_systems()

        self.form_sigma_bonds()
        self.drop_electrons()
        self.make_lone_pairs()

        self.set_atom_neighbours()

        self.make_bond_lookup()

        self.set_connectivities()
        self.set_atoms()

    def make_lone_pairs(self):
        for atom in self.graph:
            atom._set_lone_pairs()

    def promote_pi_bonds(self):
        for atom in self.graph:
            atom._promote_pi_bonds_to_d_orbitals()

    def set_connectivities(self):
        for atom in self.graph:
            atom.set_connectivity()

    def set_atom_neighbours(self):
        for atom in self.graph:
            atom._set_neighbours(self)

    def get_atom(self, atom):
        return self.atoms[atom.nr]

    def get_bond(self, bond):
        return self.bonds[bond.nr]

    def drop_electrons(self):
        for atom in self.graph:
            atom._drop_electrons()

    def add_shells_non_hydrogens(self):
        for atom in self.graph:
            if atom.type != 'H':
                atom._add_electron_shells()

    def add_shells(self):
        for atom in self.graph:
            if not atom.shells:
                atom._add_electron_shells()

    def hybridise_atoms(self, atoms=None):

        if not atoms:
            for atom in self.graph:
                atom.hybridise()

        else:
            for atom in atoms:
                atom.hybridise()

    def add_hydrogens(self):
        atom_numbers = [atom.nr for atom in self.graph]

        if len(atom_numbers) > 0:
            max_atom_nr = max(atom_numbers)
        else:
            max_atom_nr = -1

        bond_numbers = list(self.bonds.keys())

        if len(bond_numbers) > 0:
            max_bond_nr = max(bond_numbers)

        else:
            max_bond_nr = -1

        for atom in list(self.graph.keys()):
            hydrogens_to_add = atom._get_nr_implicit_hydrogens()

            for i in range(hydrogens_to_add):
                max_atom_nr += 1
                max_bond_nr += 1
                hydrogen = Atom('H', max_atom_nr, None, 0, False)
                self.add_bond(atom, hydrogen, 'single', max_bond_nr)

    def form_pi_bonds(self):
        for bond_nr in self.bonds:
            bond = self.bonds[bond_nr]
            if bond.type != 'single':
                bond.combine_p_orbitals()

    def form_sigma_bonds(self):
        for bond_nr, bond in self.bonds.items():
            bond.combine_hybrid_orbitals()

    def add_disconnected_atom(self, atom):
        self.graph[atom] = []
        self.atoms[atom.nr] = atom

    def sort_by_nr(self):
        for atom in self.graph:
            self.graph[atom].sort(key=lambda x: x.nr)

    def add_bond(self, atom_1, atom_2, bond_type, bond_nr, chiral_symbol=None):
        if atom_1 in self.graph:
            self.graph[atom_1].append(atom_2)
        else:
            self.graph[atom_1] = [atom_2]

        if atom_2 in self.graph:
            self.graph[atom_2].append(atom_1)
        else:
            self.graph[atom_2] = [atom_1]

        bond = Bond(atom_1, atom_2, bond_type, bond_nr)

        bond.chiral_symbol = chiral_symbol

        atom_1._add_bond(bond)
        atom_2._add_bond(bond)

        self.bonds[bond_nr] = bond

        if atom_1 not in self.bond_lookup:
            self.bond_lookup[atom_1] = {}
        if atom_2 not in self.bond_lookup:
            self.bond_lookup[atom_2] = {}

        self.bond_lookup[atom_1][atom_2] = bond
        self.bond_lookup[atom_2][atom_1] = bond

    def find_pi_subgraph(self, prune=True):
        pi_subgraph = {}

        for bond in self.bonds.values():
            if bond.type == 'aromatic':

                # prune the subgraph as kekulisation can only occur in atoms
                # that have an unpaired electron

                unpaired_electrons_1 = 0
                unpaired_electrons_2 = 0

                if len(bond.aromatic_system.get_contributed_electrons(bond.atom_1)) == 1:
                    unpaired_electrons_1 += 1

                if len(bond.aromatic_system.get_contributed_electrons(bond.atom_2)) == 1:
                    unpaired_electrons_2 += 1

                if unpaired_electrons_1 and unpaired_electrons_2:

                    if bond.atom_1 not in pi_subgraph:
                        pi_subgraph[bond.atom_1] = []
                    if bond.atom_2 not in pi_subgraph:
                        pi_subgraph[bond.atom_2] = []

                    pi_subgraph[bond.atom_1].append(bond.atom_2)
                    pi_subgraph[bond.atom_2].append(bond.atom_1)

                elif not prune:

                    if bond.atom_1 not in pi_subgraph:
                        pi_subgraph[bond.atom_1] = []
                    if bond.atom_2 not in pi_subgraph:
                        pi_subgraph[bond.atom_2] = []

                    pi_subgraph[bond.atom_1].append(bond.atom_2)
                    pi_subgraph[bond.atom_2].append(bond.atom_1)

        return pi_subgraph

    def kekulise(self):

        kekule_structure = self.deepcopy()

        pruned = kekule_structure.find_pi_subgraph(prune=True)
        unpruned = kekule_structure.find_pi_subgraph(prune=False)

        aromatic_unmatched = set(unpruned.keys()) - set(pruned.keys())

        matching = Match.from_structure(kekule_structure)
        unmatched_nodes = matching.unmatched_nodes()
        if unmatched_nodes != 0:
            raise Exception("This structure cannot be kekulised!")
        else:
            double_bond_pairs = set()
            single_bond_pairs = set()

            for node in matching.nodes:
                double_bond_pair = tuple(sorted([node.atom, node.mate.atom], key=lambda x: x.nr))
                if double_bond_pair not in double_bond_pairs:
                    double_bond_pairs.add(double_bond_pair)

                for neighbour in node.neighbors:
                    if neighbour.index != node.mate.index:
                        single_bond_pair = tuple(sorted([node.atom, neighbour.atom], key=lambda x: x.nr))
                        if single_bond_pair not in single_bond_pairs:
                            single_bond_pairs.add(single_bond_pair)

            # heteroatoms containing lone pairs, sp2-hybridised carbons

            for atom in aromatic_unmatched:
                for neighbour in atom.neighbours:
                    if neighbour in atom.aromatic_system.atoms:
                        single_bond_pair = tuple(sorted([atom, neighbour], key=lambda x: x.nr))
                        if single_bond_pair not in single_bond_pairs:
                            single_bond_pairs.add(single_bond_pair)

        for aromatic_system in kekule_structure.aromatic_systems:
            aromatic_system.relocalise_electrons()

        for pair in double_bond_pairs:

            new_atom_1 = kekule_structure.atoms[pair[0].nr]
            new_atom_2 = kekule_structure.atoms[pair[1].nr]
            
            bond = kekule_structure.bond_lookup[new_atom_1][new_atom_2]
            bond.type = 'double'
            bond.aromatic = False
            
            bond.atom_1.aromatic = False
            bond.atom_2.aromatic = False

            bond.set_bond_summary()
            
            orbitals_1 = new_atom_1._get_orbitals('p')
            orbitals_2 = new_atom_2._get_orbitals('p')
            
            if orbitals_1 and orbitals_2:
                orbital_1 = orbitals_1[0]
                orbital_2 = orbitals_2[0]
                
                if not len(orbital_1.electrons) == 1 or not len(orbital_2.electrons) == 1:
                    raise KekulisationError(bond.aromatic_system.__repr__())

                orbital_1.add_electron(orbital_2.electrons[0])
                orbital_2.add_electron(orbital_1.electrons[0])

                orbital_1.set_bond(bond, 'pi')
                orbital_2.set_bond(bond, 'pi')

                bond.electrons.append(orbital_1.electrons[0])
                bond.electrons.append(orbital_2.electrons[0])

                bond.set_bond_summary()
                bond.aromatic_system = None

        for pair in single_bond_pairs:
            new_atom_1 = kekule_structure.atoms[pair[0].nr]
            new_atom_2 = kekule_structure.atoms[pair[1].nr]

            bond = kekule_structure.bond_lookup[new_atom_1][new_atom_2]
            bond.type = 'single'

            bond.aromatic = False
            bond.atom_1.aromatic = False
            bond.atom_2.aromatic = False
            bond.atom_1.pyrrole = False
            bond.atom_2.pyrrole = False
            bond.atom_1.furan = False
            bond.atom_2.furan = False
            bond.atom_1.thiophene = False
            bond.atom_2.thiophene = False

            bond.aromatic_system = None

            bond.set_bond_summary()

        kekule_structure.aromatic_systems = []
        kekule_structure.aromatic_cycles = []

        return kekule_structure

    def find_start_nodes(self, paths):
        """Return atoms that still have outgoing bonds within an existing path

        Input:
        paths: list of [atom, ->], with each atom a tuple of (str, int), with
            str atom type and int atom number
        bond_dict: dict of {atom: remaining_bonds, ->}, with atom tuple of
            (str, int), with str atom type and int atom number, and remaining
            bonds int

        Output:
        start_atoms: list of [atom, ->], with each atom a tuple of (str, int),
            with str atom type and int atom number


        """

        start_atoms = []
        for path in paths:
            for atom in path:
                if self.bond_nr_dict[atom] > 0:
                    start_atoms.append(atom)

        return start_atoms

    def remove_connectors(self):
        """Remove nodes that only have incoming edges from graph

        Input:
        working_graph: dict of {node: [node, ->], ->}, representing a graph

        """

        for node in self.graph:
            for next_node in self.graph[node]:
                if next_node not in self.graph:
                    self.graph[node].remove(next_node)

    def find_new_start_node(self):
        """Return list of nodes that still have outgoing edges

        Input:
        edges_dict: dict of {node: int, ->}, with int representing the number of
            outgoing edges of that node

        Output:
        start_nodes: list of [node, ->], with each node an immutable object
        """
        start_nodes = []
        for atom in self.graph:
            if self.bond_nr_dict[atom] > 0:
                start_nodes.append(atom)

        return start_nodes