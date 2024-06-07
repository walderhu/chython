
def hide_hydrogens(self) -> None:
if not self.options.draw_hydrogens:
    for atom in self.structure.graph:


        if atom.type!= 'H':
            continue

        neighbour = atom.neighbours[0]
        neighbour.draw.has_hydrogen = True
        atom.draw.is_drawn = False
        if len(neighbour.draw.rings) < 2 and neighbour.draw.bridged_ring is None or \
        (neighbour.draw.bridged_ring is not None and len(neighbour.draw.rings) < 2):
            atom.draw.is_drawn = False
        else:
            atom.draw.is_drawn = True 
    self.structure.refresh_structure()
    self.drawn_atoms = [atom for atom in self.structure.graph if atom.draw.is_drawn]
    self.drawn_bonds = [bond for bond in self.structure.bonds.values()\
                        if bond.atom_1.draw.is_drawn and bond.atom_2.draw.is_drawn]







def hide_hydrogens(self) -> None:
    if not self.options.draw_hydrogens:
        for atom in self.structure.graph:
            if atom.type!= 'H':
                continue

            # Hydrogens should have only one neighbour, so just take the first
            # Also set has_hydrogen true on connected atom
            neighbour = atom.neighbours[0]
            neighbour.draw.has_hydrogen = True

            # Check if the neighbour is not a stereocenter or is part of less than 2 rings or is a bridged ring with less than 2 original rings
            if not neighbour.atom._adjacent_to_stereobond or len(neighbour.draw.rings) < 2 and neighbour.draw.bridged_ring is None or \
                    (neighbour.draw.bridged_ring is not None and len(neighbour.draw.original_rings) < 2):
                atom.draw.is_drawn = False
            else:
                atom.draw.is_drawn = True

        self.structure.refresh_structure()
        self.drawn_atoms = [atom for atom in self.structure.graph if atom.draw.is_drawn]
        self.drawn_bonds = [bond for bond in self.structure.bonds.values()
                           if bond.atom_1.draw.is_drawn and bond.atom_2.draw.is_drawn]