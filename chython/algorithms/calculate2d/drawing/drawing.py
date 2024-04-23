#!/usr/bin/env python
from typing import Union, List, Tuple, Dict, Set, Generator, Optional, TYPE_CHECKING
import math
from ..drawing.rings import Ring, RingOverlap, find_neighbouring_rings, rings_connected_by_bridge, find_bridged_systems
from ..math_functions import Vector, Polygon
from ..chem.atom import Atom
from ..chem.bond import Bond
from ..chem.structure import Structure
if TYPE_CHECKING:
    pass
class Options:

    def __init__(self):
        self.width: int = 500
        self.height: int = 500
        self.bond_thickness: int = 2
        self.bond_length: int = 15
        self.chiral_bond_width: float = self.bond_length * 0.1
        self.bond_length_squared: int = self.bond_length ** 2
        self.short_bond_length: float = 0.50
        self.double_bond_length: float = 0.80
        self.bond_spacing: float = 0.18 * self.bond_length
        self.isomeric: bool = True
        self.padding: int = 30
        self.font_size_large: int = 5
        self.font_size_small: int = 3
        self.kk_threshold: float = 0.1
        self.kk_inner_threshold: float = 0.1
        self.kk_max_iteration: int = 2000
        self.kk_max_inner_iteration: int = 50
        self.kk_max_energy: float = 1e9
        self.overlap_sensitivity: float = 0.10
        self.overlap_resolution_iterations: int = 5
        self.background_color: str = 'white'
        self.draw_hydrogens: bool = False
        self.finetune: bool = True
        self.strict_mode: bool = False
        self.svg_font: str = "verdana"
        self.svg_font_size: int = 8
        self.svg_font_size_small: int = 6
        self.svg_letter_spacing: int = -2


class Drawer:
    def __init__(self, structure: Structure, options: Union[Options, None] = None,
                 coords_only: bool = False, multiple: bool = False, kekulise: bool = True) -> None:
        """
        Parameters
        ----------
        structure: Structure instance
        options: Options instance, use this to define custom drawing parameters
        coords_only: bool, if True, only assign positions to atoms; draw the structure otherwise
        multiple: bool, if True, multiple structures are drawn, False otherwise.
             Do not toggle this parameter. If you want to draw multiple structures, use the function
             draw_multiple instead
        kekulise: bool, if True, kekulise the structure before drawing, False otherwise
        """
        if options is None:
            self.options: Options = Options()
        else:
            self.options: Options = options
        if kekulise:
            self.structure = structure.kekulise()
        else:
            self.structure = structure
        self.multiple: bool = multiple
        # Used for tracking rings in the structure
        self.rings: List[Ring] = []
        self.ring_overlaps: List[RingOverlap] = []
        self.id_to_ring: Dict[int, Ring] = {}
        self.has_bridged_ring: bool = False
        self.ring_id_tracker: int = 0
        self.ring_overlap_id_tracker: int = 0
        # Used for tracking rings in the structure once bridged ring systems have been established
        self.original_rings: List[Ring] = []
        self.original_ring_overlaps: List[RingOverlap] = []
        self.drawn_atoms: List[Atom] = []
        self.drawn_bonds: List[Bond] = []
        self.total_overlap_score: float = 0.0
        # TODO: Check if we cannot use the atom dictionary in Structure instead
        self.atom_nr_to_atom: Dict[int, Atom] = {}
        # Used to keep track of bonds that are under steric restraint
        self.chiral_bonds: List[Bond] = []
        self.chiral_bond_to_orientation: Dict[Bond, Tuple[str, Atom]] = {}
        self.fixed_chiral_bonds: Set[Bond] = set()
        # Used for grouping atoms and bonds within an SVG file
        self.svg_groups: Dict[str, Dict[str, List[str]]] = {}
        self.annotation_to_colour: Optional[Dict[str, str]] = None
        self.annotation: Optional[str] = None
        self.structure_id: Optional[str] = None
        # self.svg_style: str = """<style> line {stroke: black; stroke_width: 1px;} polygon {fill: black;} </style>"""
        self.svg_style: str = ""
        self.draw(coords_only)

    def set_annotation_for_grouping(self, annotation: str) -> None:
        self.svg_groups = {}
        self.annotation = annotation
    
    def find_shortest_path(self, atom_1: Atom, atom_2: Atom, path_type: str = 'bond') -> List[Union[Bond, Atom]]:
        """
        Return the shortest path between two atoms, either as a list of bonds or as a list of atoms
        Parameters
        ----------
        atom_1: Atom instance, must be in structure
        atom_2: Atom instance, must be in structure
        path_type: str, 'bond' or 'atom'
        Returns
        list of atoms or bonds describing the shortest path between atom_1 and atom_2
        """
        distances: Dict[Atom, float] = {}
        previous_hop: Dict[Atom, Optional[Atom]] = {}
        unvisited: Set[Atom] = set()

        for atom in self.structure.graph:
            distances[atom] = float('inf')
            previous_hop[atom] = None
            unvisited.add(atom)
        distances[atom_1] = 0.0
        while unvisited:
            current_atom: Optional[Atom] = None
            minimum = float('inf')
            # Find the atom with the smallest distance value that has not yet been visited
            for atom in unvisited:
                dist = distances[atom]
                if dist < minimum:
                    current_atom = atom
                    minimum = dist
            if current_atom is None:
                break            
            if current_atom == atom_2:
                break
            unvisited.remove(current_atom)
            # If there exists a shorter path between the source atom and the neighbours, update distance
            for neighbour in self.structure.graph[current_atom]:
                if neighbour in unvisited:
                    alternative_distance: float = distances[current_atom] + 1.0

                    if alternative_distance < distances[neighbour]:
                        distances[neighbour] = alternative_distance
                        previous_hop[neighbour] = current_atom
        # Construct the path of atoms
        path_atoms: List[Atom] = []
        current_atom: Optional[Atom] = atom_2
        if previous_hop[current_atom] or current_atom == atom_1:
            while current_atom:
                path_atoms.insert(0, current_atom)
                current_atom = previous_hop[current_atom]
        if path_type == 'bond':
            path: List[Union[Bond, Atom]] = []
            for i in range(1, len(path_atoms)):
                atom_1 = path_atoms[i - 1]
                atom_2 = path_atoms[i]
                bond = self.structure.bond_lookup[atom_1][atom_2]
                path.append(bond)

            return path
        elif path_type == 'atom':
            return path_atoms
        else:
            raise ValueError("Path type must be 'bond' or 'atom'.")
        
    def _finetune_overlap_resolution(self) -> None:
        if self.total_overlap_score > self.options.overlap_sensitivity:
            clashing_atoms = self._find_clashing_atoms()
            best_bonds: List[Bond] = []
            for atom_1, atom_2 in clashing_atoms:
                # shortest_path = self.find_shortest_path(atom_1, atom_2)
                if self.structure.is_connected(atom_1, atom_2):
                    shortest_path = self.find_shortest_path(atom_1, atom_2)
                rotatable_bonds: List[Bond] = []
                distances: List[float] = []
                for i, bond in enumerate(shortest_path):
                    distance_1: int = i
                    distance_2: int = len(shortest_path) - i
                    average_distance = len(shortest_path) / 2
                    distance_metric = abs(average_distance - distance_1) + abs(average_distance - distance_2)
                    if self.bond_is_rotatable(bond):
                        rotatable_bonds.append(bond)
                        distances.append(distance_metric)
                best_bond: Optional[Bond] = None
                optimal_distance: float = float('inf')
                for i, distance in enumerate(distances):
                    if distance < optimal_distance:
                        best_bond = rotatable_bonds[i]
                        optimal_distance = distance
                if best_bond is not None:
                    best_bonds.append(best_bond)
            best_bonds = list(set(best_bonds))
            for best_bond in best_bonds:
                if self.total_overlap_score > self.options.overlap_sensitivity:
                    subtree_size_1 = self.get_subgraph_size(best_bond.atom_1, {best_bond.atom_2})
                    subtree_size_2 = self.get_subgraph_size(best_bond.atom_2, {best_bond.atom_1})
                    if subtree_size_1 < subtree_size_2:
                        rotating_atom = best_bond.atom_1
                        parent_atom = best_bond.atom_2
                    else:
                        rotating_atom = best_bond.atom_2
                        parent_atom = best_bond.atom_1
                    overlap_score, _, _ = self.get_overlap_score()
                    scores: List[float] = [overlap_score]
                    # Attempt 12 rotations
                    for i in range(12):
                        self.rotate_subtree(rotating_atom, parent_atom, math.radians(30), parent_atom.draw.position)
                        new_overlap_score, _, _ = self.get_overlap_score()
                        scores.append(new_overlap_score)
                    assert len(scores) == 13
                    scores = scores[:12]
                    best_i = 0
                    best_score = scores[0]
                    for i, score in enumerate(scores):
                        if score < best_score:
                            best_score = score
                            best_i = i
                    self.total_overlap_score = best_score
                    self.rotate_subtree(rotating_atom, parent_atom, math.radians(30 * best_i + 1),
                                        parent_atom.draw.position)
   
    def _find_clashing_atoms(self) -> List[Tuple[Atom, Atom]]:
        clashing_atoms: List[Tuple[Atom, Atom]] = []
        for i, atom_1 in enumerate(self.drawn_atoms):
            for j in range(i + 1, len(self.drawn_atoms)):
                atom_2 = self.drawn_atoms[j]
                if not self.structure.bond_exists(atom_1, atom_2):
                    distance = Vector.subtract_vectors(atom_1.draw.position, atom_2.draw.position).get_squared_length()
                    if distance < 0.8 * self.options.bond_length_squared:
                        clashing_atoms.append((atom_1, atom_2))
        return clashing_atoms
    
    def set_chiral_bonds(self) -> None:
        for atom in self.structure.graph:
            if atom.chiral:
                bond, wedge = self._determine_chirality(atom)
                self.chiral_bonds.append(bond)
                self.chiral_bond_to_orientation[bond] = (wedge, atom)

    def draw(self, coords_only: bool = False) -> None:
        if not self.options.draw_hydrogens:
            self.hide_hydrogens()
        self.get_atom_nr_to_atom()
        self.define_rings()
        if not self.multiple:
            self._process_structure()
            self.set_chiral_bonds()
        else:
            self.restore_ring_information()
    
    @staticmethod
    def _chiral_bond_drawn_correctly(bond: Bond) -> bool:
        assert bond.chiral
        must_be_fixed = False
        for neighbour_1 in bond.atom_1.neighbours:
            if neighbour_1 != bond.atom_2:
                for neighbour_2 in bond.atom_2.neighbours:
                    if neighbour_2 != bond.atom_1:
                        if neighbour_1.draw.is_drawn and neighbour_2.draw.is_drawn:
                            placement_1 = Vector.get_position_relative_to_line(bond.atom_1.draw.position,
                                                                               bond.atom_2.draw.position,
                                                                               neighbour_1.draw.position)
                            placement_2 = Vector.get_position_relative_to_line(bond.atom_1.draw.position,
                                                                               bond.atom_2.draw.position,
                                                                               neighbour_2.draw.position)
                            orientation = bond.chiral_dict[neighbour_1][neighbour_2]
                            if orientation == 'cis':
                                if placement_1 != placement_2:
                                    must_be_fixed = True
                            else:
                                if placement_1 == placement_2:
                                    must_be_fixed = True
        if must_be_fixed:
            return False
        else:
            return True
    
    def _fix_chiral_bond(self, double_bond: Bond) -> None:
        if len(double_bond.atom_1.draw.rings) and len(double_bond.atom_2.draw.rings) and \
                len(set(double_bond.atom_1.draw.rings).intersection(set(double_bond.atom_2.draw.rings))) >= 1:
            self._flip_stereobond_in_ring(double_bond)
        else:
            if len(double_bond.atom_1.draw.rings) > 0:
                parent_atom = double_bond.atom_2
                root_atom = double_bond.atom_1
            else:
                parent_atom = double_bond.atom_1
                root_atom = double_bond.atom_2
            neighbours = parent_atom.drawn_neighbours[:]
            neighbours.remove(root_atom)
            if len(neighbours) == 1:
                neighbour = neighbours[0]
                self._flip_subtree(neighbour, root_atom, parent_atom)
            elif len(neighbours) == 2 and len(
                    set(neighbours[0].draw.rings).intersection(set(neighbours[1].draw.rings))) >= 1:
                self._flip_subtree(neighbours[0], root_atom, parent_atom)
            elif len(neighbours) == 2:
                neighbour_1 = neighbours[0]
                neighbour_2 = neighbours[1]
                self._flip_subtree(neighbour_1, root_atom, parent_atom)
                self._flip_subtree(neighbour_2, root_atom, parent_atom)
            self.fixed_chiral_bonds.add(double_bond)
    
    def _fix_chiral_bonds_in_rings(self) -> None:
        double_bond_sequences = self.structure.find_double_bond_sequences()
        for double_bond_sequence in double_bond_sequences:
            for double_bond in double_bond_sequence:
                chirality_correct = self._chiral_bond_drawn_correctly(double_bond)
                if chirality_correct:
                    self.fixed_chiral_bonds.add(double_bond)
                else:
                    self._fix_chiral_bond(double_bond)
        for bond in self.structure.bonds.values():
            if bond.type == 'double' and bond.chiral and bond not in self.fixed_chiral_bonds:
                chirality_correct = self._chiral_bond_drawn_correctly(bond)
                if chirality_correct:
                    self.fixed_chiral_bonds.add(bond)
                else:
                    self._fix_chiral_bond(bond)
  
    @staticmethod
    def _is_terminal(atom: Atom) -> bool:
        if len(atom.drawn_neighbours) <= 1:
            return True
        return False
    
    def _process_structure(self) -> None:
        self._position()
        self.structure.refresh_structure()
        self.restore_ring_information()
        self.restore_ring_information()
        self._fix_chiral_bonds_in_rings()
        self.resolve_primary_overlaps()
        self.total_overlap_score, sorted_overlap_scores, atom_to_scores = self.get_overlap_score()
        for i in range(self.options.overlap_resolution_iterations):
            for bond in self.drawn_bonds:
                if self.can_rotate_around_bond(bond):
                    tree_depth_1: int = self.get_subgraph_size(bond.atom_1, {bond.atom_2})
                    tree_depth_2: int = self.get_subgraph_size(bond.atom_2, {bond.atom_1})
                    atom_1_rotatable: bool = True
                    atom_2_rotatable: bool = True
                    for neighbouring_bond in bond.atom_1.bonds:
                        if neighbouring_bond.type == 'double' and neighbouring_bond.chiral:
                            atom_1_rotatable = False
                    for neighbouring_bond in bond.atom_2.bonds:
                        if neighbouring_bond.type == 'double' and neighbouring_bond.chiral:
                            atom_2_rotatable = False
                    if not atom_1_rotatable and not atom_2_rotatable:
                        continue
                    elif atom_1_rotatable and not atom_2_rotatable:
                        atom_2 = bond.atom_1
                        atom_1 = bond.atom_2
                    elif atom_2_rotatable and not atom_1_rotatable:
                        atom_1 = bond.atom_1
                        atom_2 = bond.atom_2
                    else:
                        atom_1 = bond.atom_2
                        atom_2 = bond.atom_1
                        if tree_depth_1 > tree_depth_2:
                            atom_1 = bond.atom_1
                            atom_2 = bond.atom_2
                    subtree_overlap_score, _ = self.get_subtree_overlap_score(atom_2, atom_1, atom_to_scores)
                    if subtree_overlap_score > self.options.overlap_sensitivity:
                        neighbours_2 = atom_2.drawn_neighbours[:]
                        neighbours_2.remove(atom_1)
                        if len(neighbours_2) == 1:
                            neighbour = neighbours_2[0]
                            angle = neighbour.draw.position.get_rotation_away_from_vector(atom_1.draw.position,
                                                                                          atom_2.draw.position,
                                                                                          math.radians(120))
                            self.rotate_subtree(neighbour, atom_2, angle, atom_2.draw.position)
                            new_overlap_score, _, _ = self.get_overlap_score()
                            if new_overlap_score > self.total_overlap_score:
                                self.rotate_subtree(neighbour, atom_2, -angle, atom_2.draw.position)
                            else:
                                self.total_overlap_score = new_overlap_score
                        elif len(neighbours_2) == 2:
                            if atom_2.draw.rings and atom_1.draw.rings:
                                continue
                            neighbour_1: Atom = neighbours_2[0]
                            neighbour_2: Atom = neighbours_2[1]
                            if len(neighbour_1.draw.rings) == 1 and len(neighbour_2.draw.rings) == 1:
                                if neighbour_1.draw.rings[0] != neighbour_2.draw.rings[0]:
                                    continue
                            elif neighbour_1.draw.rings or neighbour_2.draw.rings:
                                continue
                            else:
                                angle_1 = neighbour_1.draw.position.get_rotation_away_from_vector(atom_1.draw.position,
                                                                                                  atom_2.draw.position,
                                                                                                  math.radians(120))
                                angle_2 = neighbour_2.draw.position.get_rotation_away_from_vector(atom_1.draw.position,
                                                                                                  atom_2.draw.position,
                                                                                                  math.radians(120))
                                self.rotate_subtree(neighbour_1, atom_2, angle_1, atom_2.draw.position)
                                self.rotate_subtree(neighbour_2, atom_2, angle_2, atom_2.draw.position)
                                new_overlap_score, _, _ = self.get_overlap_score()
                                if new_overlap_score > self.total_overlap_score:
                                    self.rotate_subtree(neighbour_1, atom_2, -angle_1, atom_2.draw.position)
                                    self.rotate_subtree(neighbour_2, atom_2, -angle_2, atom_2.draw.position)
                                else:
                                    self.total_overlap_score = new_overlap_score
                        self.total_overlap_score, sorted_overlap_scores, atom_to_scores = self.get_overlap_score()
        for _ in range(self.options.overlap_resolution_iterations):
            if self.options.finetune:
                self._finetune_overlap_resolution()
                self.total_overlap_score, sorted_overlap_scores, atom_to_scores = self.get_overlap_score()
        for i in range(self.options.overlap_resolution_iterations):
            self.resolve_secondary_overlaps(sorted_overlap_scores)
    
    def _position(self) -> None:
        start_atom: Optional[Atom] = None
        for atom in self.structure.graph:
            if atom.draw.bridged_ring is not None:
                start_atom = atom
                break
        for ring in self.rings:
            if ring.bridged:
                start_atom = ring.members[0]
        if len(self.rings) > 0 and start_atom is None:
            start_atom = self.rings[0].members[0]
        if start_atom is None:
            for atom in self.drawn_atoms:
                if self._is_terminal(atom):
                    start_atom = atom
                    break
        if start_atom is None:
            if not self.drawn_atoms:
                raise ValueError("No atoms to draw")
            else:
                start_atom = list(self.drawn_atoms)[0]
        self.create_next_bond(start_atom, None, 0.0)
    
    def create_next_bond(self, atom: Atom, previous_atom: Optional[Atom] = None, angle: float = 0.0,
                         previous_branch_shortest: bool = False) -> None:
        if atom.draw.positioned:
            return
        if previous_atom is None:
            dummy: Vector = Vector(self.options.bond_length, 0)
            dummy.rotate(math.radians(-60.0))
            atom.draw.previous_position = dummy
            atom.draw.previous_atom = None
            atom.draw.set_position(Vector(self.options.bond_length, 0))
            atom.draw.angle = math.radians(-60.0)
            if atom.draw.bridged_ring is None:
                atom.draw.positioned = True
        elif len(previous_atom.draw.rings) > 0:
            neighbours = previous_atom.drawn_neighbours
            joined_vertex = None
            position: Vector = Vector(0, 0)
            if previous_atom.draw.bridged_ring is None and len(previous_atom.draw.rings) > 1:
                for neighbour in neighbours:
                    if len(set(neighbour.draw.rings) & set(previous_atom.draw.rings)) == len(previous_atom.draw.rings):
                        joined_vertex = neighbour
                        break
            if not joined_vertex:
                for neighbour in neighbours:
                    if neighbour.draw.positioned and self.atoms_are_in_same_ring(neighbour, previous_atom):
                        position.add(Vector.subtract_vectors(neighbour.draw.position, previous_atom.draw.position))
                position.invert()
                position.normalise()
                position.multiply_by_scalar(self.options.bond_length)
                position.add(previous_atom.draw.position)
            else:
                position = joined_vertex.draw.position.copy()
                position.rotate_around_vector(math.pi, previous_atom.draw.position)
            atom.draw.set_previous_position(previous_atom)
            atom.draw.set_position(position)
            atom.draw.positioned = True
        else:
            position: Vector = Vector(self.options.bond_length, 0)
            position.rotate(angle)
            position.add(previous_atom.draw.position)
            atom.draw.set_position(position)
            atom.draw.set_previous_position(previous_atom)
            atom.draw.positioned = True
        if len(atom.draw.rings) > 0:
            if atom.draw.bridged_ring:
                next_ring: Ring = self.id_to_ring[atom.draw.bridged_ring]
            else:
                next_ring: Ring = self.id_to_ring[atom.draw.rings[0]]
            if not next_ring.positioned:
                next_center = Vector.subtract_vectors(atom.draw.previous_position, atom.draw.position)
                next_center.invert()
                next_center.normalise()
                radius: float = Polygon.find_polygon_radius(self.options.bond_length, len(next_ring.members))
                next_center.multiply_by_scalar(radius)
                next_center.add(atom.draw.position)
                self.create_ring(next_ring, next_center, atom)
        else:
            neighbours: List[Atom] = atom.drawn_neighbours[:]
            if previous_atom:
                if previous_atom in neighbours:
                    neighbours.remove(previous_atom)
            previous_angle: float = atom.draw.get_angle()
            if len(neighbours) == 1:
                next_atom: Atom = neighbours[0]
                current_bond: Bond = self.structure.bond_lookup[atom][next_atom]
                previous_bond: Optional[Bond] = None
                if previous_atom:
                    previous_bond = self.structure.bond_lookup[previous_atom][atom]
                if current_bond.type == 'triple' or (previous_bond and previous_bond.type == 'triple') or \
                        (current_bond.type == 'double' and previous_bond and previous_bond.type == 'double' and
                         previous_atom and len(previous_atom.draw.rings) == 0 and
                         len(atom.neighbours) == 2):
                    if current_bond.type == 'double' and previous_bond.type == 'double':
                        atom.draw.draw_explicit = True
                    if current_bond.type == 'triple':
                        atom.draw.draw_explicit = True
                        next_atom.draw.draw_explicit = True
                    if previous_atom:
                        previous_bond.draw.center = True
                    current_bond.draw.center = True
                    if current_bond.type == 'double' or current_bond.type == 'triple' or \
                            (previous_atom and previous_bond.type == 'triple'):
                        next_atom.draw.angle = 0.0
                    self.create_next_bond(next_atom, atom, previous_angle + next_atom.draw.angle)
                elif previous_atom and len(previous_atom.draw.rings) > 0:
                    proposed_angle_1: float = math.radians(60.0)
                    proposed_angle_2: float = proposed_angle_1 * -1
                    proposed_vector_1 = Vector(self.options.bond_length, 0)
                    proposed_vector_2 = Vector(self.options.bond_length, 0)
                    proposed_vector_1.rotate(proposed_angle_1 + atom.draw.get_angle())
                    proposed_vector_2.rotate(proposed_angle_2 + atom.draw.get_angle())
                    proposed_vector_1.add(atom.draw.position)
                    proposed_vector_2.add(atom.draw.position)
                    centre_of_mass: Vector = self.get_current_centre_of_mass()
                    distance_1: float = proposed_vector_1.get_squared_distance(centre_of_mass)
                    distance_2: float = proposed_vector_2.get_squared_distance(centre_of_mass)
                    if distance_1 < distance_2:
                        next_atom.draw.angle = proposed_angle_2
                    else:
                        next_atom.draw.angle = proposed_angle_1
                    self.create_next_bond(next_atom, atom, previous_angle + next_atom.draw.angle)
                else:
                    proposed_angle: float = atom.draw.angle
                    if previous_atom and len(previous_atom.drawn_neighbours) > 3:
                        if round(proposed_angle, 2) > 0.00:
                            proposed_angle = min([math.radians(60), proposed_angle])
                        elif round(proposed_angle, 2) < 0.00:
                            proposed_angle = max([-math.radians(60), proposed_angle])
                        else:
                            proposed_angle = math.radians(60)
                    elif proposed_angle is None:
                        last_angled_atom = self.get_last_atom_with_angle(atom)
                        proposed_angle = last_angled_atom.draw.angle
                        if proposed_angle is None:
                            proposed_angle = math.radians(60)
                    rotatable = True
                    if previous_atom:
                        bond = self.structure.bond_lookup[previous_atom][atom]
                        if bond.type == 'double' and bond.chiral:
                            rotatable = False
                            previous_previous_atom = previous_atom.draw.previous_atom
                            if previous_previous_atom:
                                configuration = bond.chiral_dict[previous_previous_atom][next_atom]
                                if configuration == 'cis':
                                    proposed_angle = -proposed_angle
                    if rotatable:
                        if previous_branch_shortest:
                            next_atom.draw.angle = proposed_angle
                        else:
                            next_atom.draw.angle = -proposed_angle
                    else:
                        next_atom.draw.angle = -proposed_angle
                    self.create_next_bond(next_atom, atom, previous_angle + next_atom.draw.angle)
            elif len(neighbours) == 2:
                proposed_angle = atom.draw.angle
                if not proposed_angle:
                    proposed_angle = math.radians(60)
                neighbour_1, neighbour_2 = neighbours
                subgraph_1_size = self.get_subgraph_size(neighbour_1, {atom})
                subgraph_2_size = self.get_subgraph_size(neighbour_2, {atom})
                if previous_atom:
                    subgraph_3_size = self.get_subgraph_size(previous_atom, {atom})
                else:
                    subgraph_3_size = 0
                cis_atom_index = 0
                trans_atom_index = 1
                if neighbour_2.type == 'C' and neighbour_1.type != 'C' and subgraph_2_size > 1 and subgraph_1_size < 5:
                    cis_atom_index = 1
                    trans_atom_index = 0
                elif neighbour_2.type != 'C' and neighbour_1.type == 'C' and subgraph_1_size > 1 and subgraph_2_size < 5:
                    cis_atom_index = 0
                    trans_atom_index = 1
                elif subgraph_2_size > subgraph_1_size:
                    cis_atom_index = 1
                    trans_atom_index = 0
                cis_atom = neighbours[cis_atom_index]
                trans_atom = neighbours[trans_atom_index]
                previous_branch_shortest = False
                if subgraph_3_size < subgraph_2_size and subgraph_3_size < subgraph_1_size:
                    previous_branch_shortest = True
                trans_atom.draw.angle = proposed_angle
                cis_atom.draw.angle = -proposed_angle
                cis_bond = self.structure.bond_lookup[atom][cis_atom]
                trans_bond = self.structure.bond_lookup[atom][trans_atom]
                if cis_bond.type == 'single' and trans_bond.type == 'single':
                    if previous_atom:
                        previous_bond = self.structure.bond_lookup[atom][previous_atom]
                        if previous_bond.type == 'double' and previous_bond.chiral:
                            if previous_atom.draw.previous_atom:
                                configuration_cis_atom = previous_bond.chiral_dict[previous_atom.draw.previous_atom][cis_atom]
                                if configuration_cis_atom == 'cis':
                                    trans_atom.draw.angle = -proposed_angle
                                    cis_atom.draw.angle = proposed_angle
                self.create_next_bond(trans_atom, atom, previous_angle + trans_atom.draw.angle, previous_branch_shortest)
                self.create_next_bond(cis_atom, atom, previous_angle + cis_atom.draw.angle, previous_branch_shortest)
            elif len(neighbours) == 3:
                subgraph_1_size = self.get_subgraph_size(neighbours[0], {atom})
                subgraph_2_size = self.get_subgraph_size(neighbours[1], {atom})
                subgraph_3_size = self.get_subgraph_size(neighbours[2], {atom})
                straight_atom = neighbours[0]
                left_atom = neighbours[1]
                right_atom = neighbours[2]
                if subgraph_2_size > subgraph_1_size and subgraph_2_size > subgraph_3_size:
                    straight_atom = neighbours[1]
                    left_atom = neighbours[0]
                    right_atom = neighbours[2]
                elif subgraph_3_size > subgraph_1_size and subgraph_3_size > subgraph_2_size:
                    straight_atom = neighbours[2]
                    left_atom = neighbours[0]
                    right_atom = neighbours[1]
                if previous_atom and len(previous_atom.draw.rings) < 1\
                        and len(straight_atom.draw.rings) < 1\
                        and len(left_atom.draw.rings) < 1\
                        and len(right_atom.draw.rings) < 1\
                        and self.get_subgraph_size(left_atom, {atom}) == 1\
                        and self.get_subgraph_size(right_atom, {atom}) == 1\
                        and self.get_subgraph_size(straight_atom, {atom}) > 1:
                    straight_atom.draw.angle = atom.draw.angle * -1
                    if atom.draw.angle >= 0:
                        left_atom.draw.angle = math.radians(30)
                        right_atom.draw.angle = math.radians(90)
                    else:
                        left_atom.draw.angle = math.radians(-30)
                        right_atom.draw.angle = math.radians(-90)
                else:
                    straight_atom.draw.angle = 0.0
                    left_atom.draw.angle = math.radians(90)
                    right_atom.draw.angle = math.radians(-90)
                self.create_next_bond(straight_atom, atom, previous_angle + straight_atom.draw.angle)
                self.create_next_bond(left_atom, atom, previous_angle + left_atom.draw.angle)
                self.create_next_bond(right_atom, atom, previous_angle + right_atom.draw.angle)
            elif len(neighbours) == 4:
                subgraph_1_size = self.get_subgraph_size(neighbours[0], {atom})
                subgraph_2_size = self.get_subgraph_size(neighbours[1], {atom})
                subgraph_3_size = self.get_subgraph_size(neighbours[2], {atom})
                subgraph_4_size = self.get_subgraph_size(neighbours[3], {atom})
                atom_1 = neighbours[0]
                atom_2 = neighbours[1]
                atom_3 = neighbours[2]
                atom_4 = neighbours[3]
                if subgraph_2_size > subgraph_1_size and subgraph_2_size > subgraph_3_size\
                        and subgraph_2_size > subgraph_4_size:
                    atom_1 = neighbours[1]
                    atom_2 = neighbours[0]
                elif subgraph_3_size > subgraph_1_size and subgraph_3_size > subgraph_2_size\
                        and subgraph_3_size > subgraph_4_size:
                    atom_1 = neighbours[2]
                    atom_2 = neighbours[0]
                    atom_3 = neighbours[1]
                elif subgraph_4_size > subgraph_1_size and subgraph_4_size > subgraph_2_size\
                        and subgraph_4_size > subgraph_3_size:
                    atom_1 = neighbours[3]
                    atom_2 = neighbours[0]
                    atom_3 = neighbours[1]
                    atom_4 = neighbours[2]
                atom_1.draw.angle = math.radians(-36)
                atom_2.draw.angle = math.radians(36)
                atom_3.draw.angle = math.radians(-108)
                atom_4.draw.angle = math.radians(108)
                self.create_next_bond(atom_1, atom, previous_angle + atom_1.draw.angle)
                self.create_next_bond(atom_2, atom, previous_angle + atom_2.draw.angle)
                self.create_next_bond(atom_3, atom, previous_angle + atom_3.draw.angle)
                self.create_next_bond(atom_4, atom, previous_angle + atom_4.draw.angle)
                

    def restore_ring_information(self) -> None:
        bridged_rings: List[Ring] = self.get_bridged_rings()
        self.rings: List[Ring] = []
        self.ring_overlaps: List[RingOverlap] = []
        for ring in bridged_rings:
            for subring in ring.subrings:
                self.original_rings[subring.id].center = subring.center
        for ring in self.original_rings:
            for i, atom in enumerate(ring.members):
                positioned_atom = self.atom_nr_to_atom[atom.nr]
                ring.members[i] = positioned_atom
            self.set_ring_center(ring)
            self.rings.append(ring)
        for ring_overlap in self.original_ring_overlaps:
            self.ring_overlaps.append(ring_overlap)
        for atom in self.structure.graph:
            atom.draw.restore_rings()
    
    @staticmethod
    def bond_is_rotatable(bond: Bond) -> bool:
        if bond.atom_1.draw.rings and \
                bond.atom_2.draw.rings and \
                len(set(bond.atom_1.draw.rings).intersection(set(bond.atom_2.draw.rings))) > 0:
            return False
        if bond.type != 'single':
            if bond.chiral:
                return False
            if len(bond.atom_1.drawn_neighbours) > 1 and len(bond.atom_2.drawn_neighbours) > 1:
                return False
        chiral = False
        for bond_1 in bond.atom_1.bonds:
            if bond_1.chiral:
                chiral = True
                break
        for bond_2 in bond.atom_2.bonds:
            if bond_2.chiral:
                chiral = True
                break
        if chiral:
            return False
        if bond.chiral_symbol:
            return False
        
        return True
    
    @staticmethod
    def can_rotate_around_bond(bond: Bond) -> bool:
        if bond.type != 'single':
            return False
        if len(bond.atom_1.drawn_neighbours) == 1 or len(bond.atom_2.drawn_neighbours) == 1:
            return False
        if bond.atom_1.draw.rings and \
                bond.atom_2.draw.rings and \
                len(set(bond.atom_1.draw.rings).intersection(set(bond.atom_2.draw.rings))) > 0:
            return False
        return True
    
    def resolve_primary_overlaps(self) -> None:
        overlaps: List[dict[str, Union[Atom, List[int], List[Atom]]]] = []
        resolved_atoms: dict[Atom, bool] = {}
        for atom in self.structure.graph:
            if atom.draw.is_drawn:
                resolved_atoms[atom] = False
        for ring in self.rings:
            for atom in ring.members:
                if resolved_atoms[atom]:
                    continue
                resolved_atoms[atom] = True
                if not atom._adjacent_to_stereobond():
                    non_ring_neighbours = self.get_non_ring_neighbours(atom)
                    if len(non_ring_neighbours) > 1 or (len(non_ring_neighbours) == 1 and len(atom.draw.rings) == 2):
                        overlaps.append({'common': atom,
                                         'rings': atom.draw.rings,
                                         'vertices': non_ring_neighbours})
        for overlap in overlaps:
            branches_to_adjust = overlap['vertices']
            rings = overlap['rings']
            root = overlap['common']
            if len(branches_to_adjust) == 2:
                atom_1, atom_2 = branches_to_adjust
                if not atom_1.draw.is_drawn or not atom_2.draw.is_drawn:
                    continue
                angle = (2 * math.pi - self.id_to_ring[rings[0]].get_angle()) / 6.0
                self.rotate_subtree(atom_1, root, angle, root.draw.position)
                self.rotate_subtree(atom_2, root, -angle, root.draw.position)
                total, sorted_scores, atom_to_score = self.get_overlap_score()
                subtree_overlap_atom_1_1, _ = self.get_subtree_overlap_score(atom_1, root, atom_to_score)
                subtree_overlap_atom_2_1, _ = self.get_subtree_overlap_score(atom_2, root, atom_to_score)
                total_score = subtree_overlap_atom_1_1 + subtree_overlap_atom_2_1
                self.rotate_subtree(atom_1, root, -2.0 * angle, root.draw.position)
                self.rotate_subtree(atom_2, root, 2.0 * angle, root.draw.position)
                total, sorted_scores, atom_to_score = self.get_overlap_score()
                subtree_overlap_atom_1_2, _ = self.get_subtree_overlap_score(atom_1, root, atom_to_score)
                subtree_overlap_atom_2_2, _ = self.get_subtree_overlap_score(atom_2, root, atom_to_score)
                total_score_2 = subtree_overlap_atom_1_2 + subtree_overlap_atom_2_2
                if total_score_2 > total_score:
                    self.rotate_subtree(atom_1, root, 2.0 * angle, root.draw.position)
                    self.rotate_subtree(atom_2, root, -2.0 * angle, root.draw.position)
            elif len(branches_to_adjust) == 1:
                if len(rings) == 2:
                    pass
    
    def resolve_secondary_overlaps(self, sorted_scores: List[Tuple[float, Atom]]) -> None:
        for score, atom in sorted_scores:
            if score > self.options.overlap_sensitivity:
                if len(atom.drawn_neighbours) <= 1:
                    if atom.drawn_neighbours and atom.drawn_neighbours[0]._adjacent_to_stereobond():
                        continue
                    closest_atom = self.get_closest_atom(atom)
                    drawn_neighbours = closest_atom.drawn_neighbours
                    if len(drawn_neighbours) <= 1:
                        if not closest_atom.draw.previous_position:
                            closest_position = drawn_neighbours[0].draw.position
                        else:
                            closest_position = closest_atom.draw.previous_position
                    else:
                        if not closest_atom.draw.previous_position:
                            closest_position = drawn_neighbours[0].draw.position
                        else:
                            closest_position = closest_atom.draw.position
                    if not atom.draw.previous_position:
                        atom_previous_position = atom.drawn_neighbours[0].draw.position
                    else:
                        atom_previous_position = atom.draw.previous_position
                    atom.draw.position.rotate_away_from_vector(closest_position, atom_previous_position,
                                                               math.radians(20))
    
    def get_atom_nr_to_atom(self) -> None:
        self.atom_nr_to_atom = {}
        for atom in self.structure.graph:
            self.atom_nr_to_atom[atom.nr] = atom
    
    def get_subtree_overlap_score(self, root: Atom, root_parent: Atom,
                                  atom_to_score: Dict[Atom, float]) -> Tuple[float, Vector]:
        score = 0.0
        center = Vector(0, 0)
        count = 0
        for atom in self.traverse_substructure(root, {root_parent}):
            subscore = atom_to_score[atom]
            if subscore > self.options.overlap_sensitivity:
                score += subscore
                count += 1
            position = atom.draw.position.copy()
            position.multiply_by_scalar(subscore)
            center.add(position)
        if score:
            center.divide(score)
        if count == 0:
            count = 1
        return score / count, center
    
    def get_overlap_score(self) -> Tuple[float, List[Tuple[float, Atom]], Dict[Atom, float]]:
        total = 0.0
        overlap_scores = {}
        for atom in self.drawn_atoms:
            overlap_scores[atom] = 0.0
        for i, atom_1 in enumerate(self.drawn_atoms):
            for j in range(i + 1, len(self.drawn_atoms)):
                atom_2 = self.drawn_atoms[j]
                distance = Vector.subtract_vectors(atom_1.draw.position, atom_2.draw.position).get_squared_length()
                if distance < self.options.bond_length_squared:
                    weight = (self.options.bond_length - math.sqrt(distance)) / self.options.bond_length
                    total += weight
                    overlap_scores[atom_1] += weight
                    overlap_scores[atom_2] += weight
        sorted_overlaps: List[Tuple[float, Atom]] = []
        for atom in self.drawn_atoms:
            sorted_overlaps.append((overlap_scores[atom], atom))
        sorted_overlaps.sort(key=lambda x: x[0], reverse=True)
        return total, sorted_overlaps, overlap_scores
    
    @staticmethod
    def get_non_ring_neighbours(atom: Atom) -> List[Atom]:
        non_ring_neighbours: List[Atom] = []
        for neighbour in atom.drawn_neighbours:
            nr_overlapping_rings = len(set(atom.draw.rings).intersection(set(neighbour.draw.rings)))
            if nr_overlapping_rings == 0 and not neighbour.draw.is_bridge:
                non_ring_neighbours.append(neighbour)
        return non_ring_neighbours
    
    def rotate_subtree(self, root: Atom, root_parent: Atom, angle: float, center: Vector):
        for atom in self.traverse_substructure(root, {root_parent}):
            atom.draw.position.rotate_around_vector(angle, center)
            for anchored_ring in atom.draw.anchored_rings:
                if anchored_ring.center:
                    anchored_ring.center.rotate_around_vector(angle, center)
    
    def rotate_structure(self, angle, midpoint: Optional[Vector] = None) -> None:
        if midpoint is None:
            midpoint: Vector = self.get_average_position()
        for atom in self.drawn_atoms:
            atom.draw.position.rotate_around_vector(angle, midpoint)
            for anchored_ring in atom.draw.anchored_rings:
                if anchored_ring.center:
                    anchored_ring.center.rotate_around_vector(angle, midpoint)
    
    def traverse_substructure(self, atom: Atom, visited: Set[Atom]) -> Generator[Atom, None, None]:
        yield atom
        visited.add(atom)
        for neighbour in atom.drawn_neighbours:
            if neighbour not in visited:
                yield from self.traverse_substructure(neighbour, visited)
    
    def get_subgraph_size(self, atom: Atom, masked_atoms: Set[Atom]) -> int:
        masked_atoms.add(atom)
        for neighbour in atom.drawn_neighbours:
            if neighbour not in masked_atoms:
                self.get_subgraph_size(neighbour, masked_atoms)
        return len(masked_atoms) - 1
    
    @staticmethod
    def get_last_atom_with_angle(atom: Atom) -> Optional[Atom]:
        parent_atom: Optional[Atom] = atom.draw.previous_atom
        angle: float = parent_atom.draw.angle
        while parent_atom and not angle:
            parent_atom = parent_atom.draw.previous_atom
            angle = parent_atom.draw.angle
        return parent_atom
    
    def create_ring(self, ring: Ring, center: Optional[Vector] = None, start_atom: Optional[Atom] = None,
                    previous_atom: Optional[Atom] = None) -> None:
        if ring.positioned:
            return
        if center is None:
            center = Vector(0, 0)
        ordered_neighbour_ids: List[int] = ring.get_ordered_neighbours(self.ring_overlaps)
        starting_angle: float = 0
        if start_atom:
            starting_angle = Vector.subtract_vectors(start_atom.draw.position, center).angle()
        ring_size: int = len(ring.members)
        radius: float = Polygon.find_polygon_radius(self.options.bond_length, ring_size)
        angle: float = Polygon.get_central_angle(ring_size)
        ring.central_angle = angle
        if start_atom not in ring.members:
            if start_atom:
                start_atom.draw.positioned = False
            start_atom = ring.members[0]
            ring.positioned = True
            self.set_ring_center(ring)
            center = ring.center
            for subring in ring.subrings:
                self.set_ring_center(subring)
        else:
            ring.set_member_positions(self.structure, start_atom, previous_atom, center, starting_angle, radius, angle)
        ring.positioned = True
        ring.center = center
        for neighbour_id in ordered_neighbour_ids:
            neighbour = self.id_to_ring[neighbour_id]
            if neighbour.positioned:
                continue
            atoms = list(RingOverlap.get_vertices(self.ring_overlaps, ring.id, neighbour.id))
            if len(atoms) == 2:
                ring.fused = True
                neighbour.fused = True
                atom_1 = atoms[0]
                atom_2 = atoms[1]
                midpoint = Vector.get_midpoint(atom_1.draw.position, atom_2.draw.position)
                normals = Vector.get_normals(atom_1.draw.position, atom_2.draw.position)
                normals[0].normalise()
                normals[1].normalise()
                apothem = Polygon.get_apothem_from_side_length(self.options.bond_length,
                                                               len(neighbour.members))
                normals[0].multiply_by_scalar(apothem)
                normals[1].multiply_by_scalar(apothem)
                normals[0].add(midpoint)
                normals[1].add(midpoint)
                next_center = normals[0]
                distance_to_center_1 = Vector.subtract_vectors(center, normals[0]).get_squared_length()
                distance_to_center_2 = Vector.subtract_vectors(center, normals[1]).get_squared_length()
                if distance_to_center_2 > distance_to_center_1:
                    next_center = normals[1]
                position_1 = Vector.subtract_vectors(atom_1.draw.position, next_center)
                position_2 = Vector.subtract_vectors(atom_2.draw.position, next_center)
                if position_1.get_clockwise_orientation(position_2) == 'clockwise':
                    if not neighbour.positioned:
                        self.create_ring(neighbour, next_center, atom_1, atom_2)
                else:
                    if not neighbour.positioned:
                        self.create_ring(neighbour, next_center, atom_2, atom_1)
            elif len(atoms) == 1:
                ring.spiro = True
                neighbour.spiro = True
                atom = atoms[0]
                next_center = Vector.subtract_vectors(center, atom.draw.position)
                next_center.invert()
                next_center.normalise()
                distance_to_center = Polygon.find_polygon_radius(self.options.bond_length, len(neighbour.members))
                next_center.multiply_by_scalar(distance_to_center)
                next_center.add(atom.draw.position)
                if not neighbour.positioned:
                    self.create_ring(neighbour, next_center, atom)
        for atom in ring.members:
            for neighbour in atom.drawn_neighbours:
                if neighbour.draw.positioned:
                    continue
                atom.draw.connected_to_ring = True
                self.create_next_bond(neighbour, atom, 0.0)
    
    @staticmethod
    def set_ring_center(ring: Ring) -> None:
        total: Vector = Vector(0, 0)
        for atom in ring.members:
            total.add(atom.draw.position)
        total.divide(len(ring.members))
        ring.center = total
    
    def get_current_centre_of_mass(self) -> Vector:
        total = Vector(0, 0)
        count = 0
        for atom in self.structure.graph:
            if atom.draw.positioned:
                total.add(atom.draw.position)
                count += 1
        total.divide(count)
        return total
    
    @staticmethod
    def atoms_are_in_same_ring(atom_1: Atom, atom_2: Atom) -> bool:
        for ring_id_1 in atom_1.draw.rings:
            for ring_id_2 in atom_2.draw.rings:
                if ring_id_1 == ring_id_2:
                    return True
        return False
    
    def define_rings(self) -> None:
        rings = self.structure.cycles.find_sssr()
        if not rings:
            return
        for ring_members in rings:
            ring = Ring(ring_members)
            self.add_ring(ring)
            for atom in ring_members:
                structure_atom = self.structure.get_atom(atom)
                structure_atom.draw.rings.append(ring.id)
        for i, ring_1 in enumerate(self.rings[:-1]):
            for ring_2 in self.rings[i + 1:]:
                ring_overlap = RingOverlap(ring_1, ring_2)
                if len(ring_overlap.atoms) > 0:
                    self.add_ring_overlap(ring_overlap)
        for ring in self.rings:
            neighbouring_rings = find_neighbouring_rings(self.ring_overlaps, ring.id)
            ring.neighbouring_rings = neighbouring_rings
        for ring in self.rings:
            anchor = ring.members[0]
            if ring not in anchor.draw.anchored_rings:
                anchor.draw.anchored_rings.append(ring)
        self.backup_ring_info()
        while True:
            ring_id: int = -1
            for ring in self.rings:
                if self.is_part_of_bridged_ring(ring.id) and not ring.bridged:
                    ring_id: int = ring.id
            if ring_id == -1:
                break
            ring: Ring = self.id_to_ring[ring_id]
            involved_ring_ids: Union[List[int], set[int]] = []
            self.get_bridged_ring_subrings(ring.id, involved_ring_ids)
            involved_ring_ids = set(involved_ring_ids)
            self.has_bridged_ring = True
            self.create_bridged_ring(involved_ring_ids)
            for involved_ring_id in involved_ring_ids:
                involved_ring = self.id_to_ring[involved_ring_id]
                self.remove_ring(involved_ring)
        bridged_systems = find_bridged_systems(self.rings, self.ring_overlaps)
        if bridged_systems and not self.has_bridged_ring:
            self.has_bridged_ring = True
            for bridged_system in bridged_systems:
                involved_ring_ids = set(bridged_system)
                self.create_bridged_ring(involved_ring_ids)
                for involved_ring_id in involved_ring_ids:
                    involved_ring = self.id_to_ring[involved_ring_id]
                    self.remove_ring(involved_ring)
    
    # def hide_hydrogens(self) -> None:
    #     hidden: List[Atom] = []
    #     exposed: List[Atom] = []
    #     self.structure.refresh_structure()
    #     for atom in self.structure.graph:
    #         if atom.type != 'H':
    #             continue
    #         elif atom.charge != 0:
    #             continue
    #         neighbour = atom.neighbours[0]
    #         neighbour.draw.has_hydrogen = True
    #         atom.draw.is_drawn = False
    #         hidden.append(atom)
    #         if len(neighbour.draw.rings) < 2 and neighbour.draw.bridged_ring is None and \
    #                 neighbour.draw.bridged_ring is not None and len(neighbour.draw.original_rings) < 2:
    #             atom.draw.is_drawn = False
    #             neighbour.draw.has_hydrogen = True
    #             hidden.append(atom)
    #         else:
    #             exposed.append(atom)
    #     for atom in self.structure.graph:
    #         atom.set_drawn_neighbours()
    #         if atom.type == 'O':
    #             pass
    #     self.drawn_bonds = []
    #     for bond_nr, bond in self.structure.bonds.items():
    #         if bond.atom_1.draw.is_drawn and bond.atom_2.draw.is_drawn:
    #             self.drawn_bonds.append(bond)
    #     self.drawn_atoms = self.structure.get_drawn_atoms()
    def hide_hydrogens(self) -> None:
        hidden: List[Atom] = []
        exposed: List[Atom] = []
        self.structure.refresh_structure()
        for atom in self.structure.graph:
            if atom.type != 'H':
                continue
            elif atom.charge != 0:
                continue
            neighbour = atom.neighbours[0]
            neighbour.draw.has_hydrogen = True
            # if self.should_hide_hydrogen(atom, neighbour):
            #     atom.draw.is_drawn = False
            #     hidden.append(atom)
            # else:
            #     exposed.append(atom)
            exposed.append(atom)
        for atom in self.structure.graph:
            atom.set_drawn_neighbours()
        self.drawn_bonds = []
        for bond_nr, bond in self.structure.bonds.items():
            if bond.atom_1.draw.is_drawn and bond.atom_2.draw.is_drawn:
                self.drawn_bonds.append(bond)
        self.drawn_atoms = self.structure.get_drawn_atoms()

    def should_hide_hydrogen(self, hydrogen_atom: Atom, neighbour_atom: Atom) -> bool:
        # Здесь должна быть логика, которая определяет, должен ли атом водорода быть скрыт
        # В этом примере мы просто возвращаем False, что означает, что атом водорода не будет скрыт
        # Вам нужно будет добавить здесь свою логику для определения, когда атом водорода должен быть скрыт
        return False
    
    def get_bridged_rings(self) -> List[Ring]:
        bridged_rings: List[Ring] = []
        for ring in self.rings:
            if ring.bridged:
                bridged_rings.append(ring)
        return bridged_rings
    
    def get_bridged_ring_subrings(self, ring_id: int, involved_ring_ids: List[int]) -> None:
        involved_ring_ids.append(ring_id)
        ring = self.id_to_ring[ring_id]
        for neighbour_id in ring.neighbouring_rings:
            if neighbour_id not in involved_ring_ids and neighbour_id != ring_id and \
                    rings_connected_by_bridge(self.ring_overlaps, ring_id, neighbour_id):
                self.get_bridged_ring_subrings(neighbour_id, involved_ring_ids)
   
    def backup_ring_info(self):
        self.original_rings = []
        for ring in self.rings:
            self.original_rings.append(ring.copy())
        self.original_ring_overlaps = []
        for ring_overlap in self.ring_overlaps:
            self.original_ring_overlaps.append(ring_overlap.copy())
        for atom in self.structure.graph:
            atom.draw.original_rings = []
            for ring in atom.draw.rings:
                atom.draw.original_rings.append(ring)

    def add_ring(self, ring):
        ring.id = self.ring_id_tracker
        self.rings.append(ring)
        self.id_to_ring[ring.id] = ring
        self.ring_id_tracker += 1

    def add_ring_overlap(self, ring_overlap):
        ring_overlap.id = self.ring_overlap_id_tracker
        self.ring_overlaps.append(ring_overlap)
        self.ring_overlap_id_tracker += 1

    def is_part_of_bridged_ring(self, ring_id):
        for ring_overlap in self.ring_overlaps:
            if ring_overlap.involves_ring(ring_id) and ring_overlap.is_bridge():
                return True
        return False
    
    def get_closest_atom(self, atom):
        minimal_distance = 9999999
        closest_atom = None
        for atom_2 in self.drawn_atoms:
            if atom == atom_2:
                continue
            squared_distance = atom.draw.position.get_squared_distance(atom_2.draw.position)
            if squared_distance < minimal_distance:
                minimal_distance = squared_distance
                closest_atom = atom_2
        return closest_atom
    
    