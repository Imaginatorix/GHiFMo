from Chromosome import *
import random
import json

class Mutate():
    def __init__(self, head_atom: Atom, ring_manager: Ring_Manager):
        self.head_atom = head_atom
        self.ring_manager = ring_manager
        self.atoms = []
        self.bonds = []
        self.recount()

    def recount(self):
        # Traverse the entire molecule and place the bonds, atoms, and rings into an array
        # Making each one of them have an equal chance of being chosen in the mutation
        frontier = [self.head_atom]
        while frontier:
            current_atom = frontier.pop(0)
            # Add the atom to the list
            self.atoms.append(current_atom)
            # Add the bond to the list
            for branch in current_atom.branches:
                self.bonds.append(branch)
                if type(branch.atom) != Ring_Closure:
                    frontier.append(branch.atom)


    def OpenRing(self):
        # Remove the ring bond and change it to H
        # Check if it needs re-rooting
        rings = list(self.ring_manager.rings.values())
        if len(rings) < 1:
            return self

        while len(rings) != 0:
            ring_group = rings.pop(random.randint(0, len(rings)-1))
            # Picked a removed group
            if not ring_group.is_closed():
                continue

            # Pick random atom in members
            ring_opening = random.randint(0, len(ring_group.members)-1)
            # If original ring_opening is chosen
            if ring_opening == len(ring_group.members)-1:
                # Same as ring removal
                self.ring_manager.remove_ring(ring_group.identity)
            else:
                # for i in ring_group.members:
                #     print("ATOM:", i.value)
                # print(ring_group.top_hierarchy_index, ring_opening)
                # Find the relationship between chosen ring_opening and its +1 position
                # Check if it's the highest (i.e. in the ring group)
                if ring_group.top_hierarchy_index == ring_opening:
                    parent = ring_group.members[ring_opening]
                    child = ring_group.members[ring_opening+1]
                # Check if +1 is parent
                elif ring_group.members[ring_opening].from_branch.from_atom is ring_group.members[ring_opening+1]:
                    parent = ring_group.members[ring_opening+1]
                    child = ring_group.members[ring_opening]
                # Else +1 is a child
                else:
                    parent = ring_group.members[ring_opening]
                    child = ring_group.members[ring_opening+1]

                # Remove their connection
                parent.branches.remove(child.from_branch)
                child.from_branch.from_atom = None
                child.from_branch.atom = None
                child.from_branch = None

                # Remove ring so we can re-root it
                ring_group_closure_bond = ring_group.from_branch_a.category
                self.ring_manager.remove_ring(ring_group.identity)
                # Re-root the rest of the rings to the highest atom
                # If bond removed is on right side of top
                if ring_opening >= ring_group.top_hierarchy_index:
                    # Re-root everything on the left of the ring opening
                    # Connect right-most until ring_opening+1
                    for index in range(len(ring_group.members)-1, ring_opening+1, -1):
                        parent = ring_group.members[index]
                        child = ring_group.members[index-1]

                        parent_bond = parent.from_branch.category
                        # Remove previous connection
                        child.branches.remove(parent.from_branch)
                        parent.from_branch.from_atom = None
                        parent.from_branch.atom = None
                        parent.from_branch = None
                        # Add new connection
                        parent.add_branches([Branch(parent_bond, child)])

                    # Connect left-most to right-most
                    parent = ring_group.members[0]
                    child = ring_group.members[-1]
                    parent.add_branches([Branch(ring_group_closure_bond, child)])

                # If bond removed is on left side of top
                else:
                    # Re-root everything on the right of the ring opening
                    # Connect left-most until ring_opening
                    for index in range(0, ring_opening):
                        parent = ring_group.members[index]
                        child = ring_group.members[index+1]

                        parent_bond = parent.from_branch.category
                        # Remove previous connection
                        child.branches.remove(parent.from_branch)
                        parent.from_branch.from_atom = None
                        parent.from_branch.atom = None
                        parent.from_branch = None
                        # Add new connection
                        parent.add_branches([Branch(parent_bond, child)])

                    # Connect right-most to left-most
                    parent = ring_group.members[-1]
                    child = ring_group.members[0]
                    parent.add_branches([Branch(ring_group_closure_bond, child)])
            break
        return self

    def CloseRing(self):
        # Randomly pick two atoms (with Hydrogens) and connect
        # Chosen atom must not be vicinal (adjacent) or geminal (same)
        if len(self.atoms) <= 2:
            return self

        atoms = self.atoms.copy()
        atom_a = None
        atom_b = None
        while (len(atoms) != 0) and (not atom_a or not atom_b):
            atom = atoms.pop(random.randint(0, len(atoms)-1)) # Removes geminality
            # If atom has hydrogen, finalize it as atom_a or atom_b
            if atom.valence_remain() > 0:
                if atom_a is None:
                    atom_a = atom
                else:
                    # Check if it is vicinal
                    # Check all branches of atom_a (except for the ring closure)
                    flag = False
                    for branch in atom_a.branches:
                        if type(branch.atom) != Ring_Closure and branch.atom is atom:
                            flag = True
                            break
                    if flag:
                        continue

                    # Check all branches of atom (except for the ring closure)
                    for branch in atom.branches:
                        if type(branch.atom) != Ring_Closure and branch.atom is atom_a:
                            flag = True
                            break
                    if flag:
                        continue

                    atom_b = atom

        # Connect
        if atom_a and atom_b:
            new_ring = self.ring_manager.add_ring()
            atom_a.add_branches([Branch("single", new_ring)])
            atom_b.add_branches([Branch("single", new_ring)])
        return self


    def ExchangeVicinal(self):
        # The purpose of this is to extend or contract rings
        # Exchange hydrogen with a singly bond substructure at an adjacent position
        if len(self.atoms) <= 2:
            return self

        atoms = self.atoms.copy()
        while len(atoms) != 0:
            atom = atoms.pop(random.randint(0, len(atoms)-1))
            if not atom.valence_remain() > 0:
                continue

            vicinals = []
            # Look at adjacent atoms
            # Below
            for branch_ in atom.branches:
                if type(branch_.atom) != Ring_Closure:
                    for branch in branch_.atom.branches:
                        if branch.category == "single" and type(branch.atom) != Ring_Closure:
                            vicinals.append(branch)
            # Above
            if not atom.from_branch is None:
                for branch in atom.from_branch.from_atom.branches:
                    if branch.category == "single" and type(branch.atom) != Ring_Closure and not branch.atom is atom:
                        vicinals.append(branch)

            # Exchange
            if vicinals:
                picked_branch = random.choice(vicinals)
                picked_branch.from_atom.branches.remove(picked_branch)
                atom.add_branches([picked_branch])

                # If invalid molecule then turn ring closure to Hydrogen
                if not Atom.is_valid(self.head_atom, self.ring_manager):
                    # Find ring closures that were moved, and then turn them to hydrogen
                    frontier = [picked_branch.atom]
                    while frontier:
                        cur_atom = frontier.pop(0)
                        for branch in cur_atom.branches:
                            if type(branch.atom) == Ring_Closure:
                                self.ring_manager.remove_ring(branch.atom.identity)
                            else:
                                frontier.append(branch.atom)
                break

        return self


    def CrossOver(self, mutate_head_atom):
        # Randomly pick two non-ring fragments of molecules and exchange
        if len(self.bonds) < 1 or len(mutate_head_atom.bonds) < 1:
            return self

        bonds_1 = self.bonds.copy()
        bonds_2 = mutate_head_atom.bonds.copy()
        matching_bond = False
        while len(bonds_1) != 0 and len(bonds_2) != 0:
            bond_1 = bonds_1.pop(random.randint(0, len(bonds_1)-1))
            # If parent of bond is part of ring
            if self.ring_manager.part_of_ring(bond_1.from_atom):
                continue

            while len(bonds_2) != 0:
                bond_2 = bonds_2.pop(random.randint(0, len(bonds_2)-1))
                # If parent of bond is part of ring
                if self.ring_manager.part_of_ring(bond_2.from_atom) or bond_2.category != bond_1.category:
                    continue
                matching_bond = True
                break

            if not matching_bond:
                return self

            # Exchange rings in ring_manager
            # From self to mutate_bond_atom
            # BFS to look for rings being exchanged
            frontier = [bond_1.atom]
            while frontier:
                atom = frontier.pop(0)
                if type(atom) == Ring_Closure:
                    # Remove ring from old ring_manager
                    self.ring_manager.remove_ring(atom.identity)
                    # Add ring to new ring_manager
                    mutate_head_atom.ring_manager.add_ring(ring_closure_object=atom)
                else:
                    # Add branches
                    for branch in atom.branches:
                        frontier.append(branch.atom)
            # From mutate bond atom to self
            # BFS to look for rings being exchanged
            frontier = [bond_2.atom]
            while frontier:
                atom = frontier.pop(0)
                if type(atom) == Ring_Closure:
                    # Remove ring from old ring_manager
                    mutate_head_atom.ring_manager.remove_ring(atom.identity)
                    # Add ring to new ring_manager
                    self.ring_manager.add_ring(ring_closure_object=atom)
                else:
                    # Add branches
                    for branch in atom.branches:
                        frontier.append(branch.atom)

            # Exchange bond_1 and bond_2
            self_atom = bond_1.from_atom
            other_atom = bond_2.from_atom
            # Remove connections
            self_atom.branches.remove(bond_1)
            bond_1.from_atom = None
            other_atom.branches.remove(bond_2)
            bond_2.from_atom = None
            # Reconnect
            self_atom.add_branches([bond_2])
            other_atom.add_branches([bond_1])

            break
        return self


    def MutateAtom(self):
        # Pick an atom and change it to a new one with the same valence (or +-1 still making the molecule valid)
        # All supported atoms so far in the mutation
        if len(self.atoms) < 1:
            return self

        possible_atoms = {
            1: ["F", "Cl", "Br", "I"],
            2: ["O", "S", "Se", "Hg"],
            3: ["B", "N", "Al", "P", "As"],
            4: ["C", "Si"]
        }

        atom = self.atoms.pop(random.randint(0, len(self.atoms)-1))
        atom_valence = atom.atom_valence()
        possible_mutations = possible_atoms[atom_valence]
        if atom_valence < 4:
            possible_mutations += possible_atoms[atom_valence+1]
        if atom_valence > 1 and atom.valence_remain() > 0:
            possible_mutations += possible_atoms[atom_valence-1]

        atom.value = element(random.choice(possible_mutations))
        return self


    def MutateFragment(self):
        # Pick a fragment and change it to a new one with the same length in the database
        if len(self.atoms) < 4:
            return self

        with open("data/fragments_transformations_frequency.json", "r") as f:
            all_fragments_frequency = json.load(f)
        with open("data/fragments_transformations_length.json", "r") as f:
            all_fragments_length = json.load(f)

        # Pick random length
        chosen_length = float("inf")
        mol_length = self.get_length()

        fragment_lengths = list(all_fragments_length.keys())
        while len(fragment_lengths) != 0 and chosen_length >= mol_length:
            chosen_length = int(fragment_lengths.pop(random.randint(0, len(fragment_lengths)-1)))

        if chosen_length >= mol_length:
            return self

        fragment_weights = [all_fragments_frequency[fragment] for fragment in all_fragments_length[str(chosen_length)]]
        chosen_fragment = random.choice(list(filter(lambda x: x.count("*") == 2, all_fragments_length[str(chosen_length)])))#, weights=fragment_weights)[0]

        # TODO: Fix valence issues
        try:
            chosen_fragment_hier = atom_from_smiles(chosen_fragment)
        except Exception:
            return self

        fragment_head = chosen_fragment_hier[0]
        fragment_ring_manager = chosen_fragment_hier[1]

        print(chosen_fragment)

        # If one attachment point
        if chosen_fragment.count("*") == 1:
            # Find an atom that is <length>+1 away from leaf
            # Random walk until non-ring leaf
            visited = []
            frontier = [self.head_atom]
            while frontier:
                atom = frontier.pop(0)

                # Add random branches to frontier
                atom_branches = atom.branches.copy()
                branch = None
                while len(atom_branches) != 0:
                    # Pick random branch
                    branch = atom_branches.pop(random.randint(0, len(atom_branches)-1))
                    if type(branch.atom) != Ring_Closure and not branch.atom in visited:
                        frontier.append(branch.atom)
                        visited.append(branch.atom)
                        break

            if len(visited) < chosen_length+1:
                return self

            fragment_end = visited[len(visited)-chosen_length-1]
            if not self.ring_manager.part_of_ring(fragment_end):
                # Pick len-length-1 from visited as atom
                branched_atom = visited[len(visited)-chosen_length]
                attached_branch = None
                # Exchange with chosen fragment
                for branch in fragment_end.branches:
                    if branch.atom is branched_atom:
                        attached_branch = branch
                        break

                # Remove branch
                fragment_end.branches.remove(attached_branch)

                # Add all rings in ring manager to molecule
                frontier = [fragment_head.branches[0].atom]
                while frontier:
                    cur_atom = frontier.pop(0)
                    if type(cur_atom) == Ring_Closure:
                        if cur_atom.identity in fragment_ring_manager.rings:
                            # Remove ring from old ring_manager
                            fragment_ring_manager.remove_ring(cur_atom.identity)
                            # Add ring to new ring_manager
                            self.ring_manager.add_ring(ring_closure_object=cur_atom)
                    else:
                        for branch in cur_atom.branches:
                            frontier.append(branch.atom)

                # Connect fragment to branch
                fragment_end.add_branches(fragment_head.branches)

        # If two attachment point
        elif chosen_fragment.count("*") == 2:
            # Change fragment to same length
            fragment_end_a = None
            fragment_end_b = None
            # Pick higher end
            atoms = self.atoms.copy()
            while len(atoms) != 0:
                atom = atoms.pop(random.randint(0, len(atoms)-1))
                # Must not be part of ring
                if self.ring_manager.part_of_ring(atom):
                    continue
                fragment_end_a = atom
                break

            # Pick lower end using random walk
            # BFS until same length, make sure end is not part of ring
            length = -1
            visited = []
            frontier = [fragment_end_a]
            while frontier and not fragment_end_a is None:
                atom = frontier.pop(0)
                # Add branches to frontier
                atom_branches = atom.branches.copy()
                while len(atom_branches) != 0:
                    branch = atom_branches.pop(random.randint(0, len(atom.branches)-1))
                    if type(branch.atom) != Ring_Closure:
                        frontier.append(branch.atom)
                        visited.append(branch.atom)
                        break

                length += 1
                if length == chosen_length+1:
                    fragment_end_b = visited[-1]
                    break

            # Both are identified
            if not fragment_end_a is None and not fragment_end_b is None:
                # Identify lower attachment point using BFS
                lower_attachment_point = None
                frontier = [fragment_head]
                while frontier and lower_attachment_point is None:
                    atom = frontier.pop(0)
                    for branch in atom.branches:
                        if type(branch.atom) != Ring_Closure:
                            if branch.atom.value.symbol == "*":
                                lower_attachment_point = atom
                                break
                        frontier.append(branch.atom)

                # If for some reason, attachment point is none, do nothing
                if lower_attachment_point is None:
                    return self
                # Exchange fragments
                # Connect higher end
                # Remove traversed branch
                fragment_end_a.branches.remove(visited[0].from_branch)

                # Add all rings in ring manager to molecule
                frontier = [fragment_head.branches[0].atom]
                while frontier:
                    cur_atom = frontier.pop(0)
                    if type(cur_atom) == Ring_Closure:
                        if cur_atom.identity in fragment_ring_manager.rings:
                            # Remove ring from old ring_manager
                            fragment_ring_manager.remove_ring(cur_atom.identity)
                            # Add ring to new ring_manager
                            self.ring_manager.add_ring(ring_closure_object=cur_atom)
                    else:
                        # Add branches
                        for branch in cur_atom.branches:
                            frontier.append(branch.atom)

                # Connect fragment to branch
                fragment_end_a.add_branches(fragment_head.branches)


                # Connect lower end
                # Remove the parent from the lower end
                branches_lower = fragment_end_b.from_branch.from_atom.branches
                for branch in branches_lower:
                    branch.from_atom = None

                # Remove lost rings from the original fragment
                ## Loop visited, and try to identify what rings the intermediate atoms are part of
                for atom in visited:
                    ring_groups = self.ring_manager.part_of_ring(atom, show_groups=True)
                    if len(ring_groups) != 0:
                        # Remove ring from old ring_manager
                        self.ring_manager.remove_ring(atom.identity)

                # Add fragment
                for branch in lower_attachment_point.branches.copy():
                    if type(branch.atom) != Ring_Closure and branch.atom.value.symbol == "*":
                        lower_attachment_point.branches.remove(branch)
                lower_attachment_point.add_branches(branches_lower)

        return self


    def MutateBranch(self):
        # Replace branch with a fragment in the library
        if len(self.atoms) < 1:
            return self

        with open("data/fragments_transformations_branches.json", "r") as f:
            all_fragments_branches = json.load(f)

        atoms = self.atoms.copy()
        while len(atoms) != 0:
            atom = atoms.pop(random.randint(0, len(atoms)-1))
            # Make sure atom is not part of ring and can take more bonds
            if self.ring_manager.part_of_ring(atom) or len(atom.branches) == 0:
                continue

            # Pick random fragment to connect
            fragment_smiles = random.choices(list(all_fragments_branches.keys()), weights=list(all_fragments_branches.values()))[0]

            # TODO: Fix valence issues
            try:
                fragment = atom_from_smiles(fragment_smiles)
            except Exception:
                continue

            # Remove random branch
            atom.branches.pop(random.randint(0, len(atom.branches)-1))

            fragment_head = fragment[0]
            fragment_ring_manager = fragment[1]

            # Add all rings in ring manager to molecule
            frontier = [fragment_head.branches[0].atom]
            while frontier:
                cur_atom = frontier.pop(0)
                if type(cur_atom) == Ring_Closure:
                    if cur_atom.identity in fragment_ring_manager.rings:
                        # Remove ring from old ring_manager
                        fragment_ring_manager.remove_ring(cur_atom.identity)
                        # Add ring to new ring_manager
                        self.ring_manager.add_ring(ring_closure_object=cur_atom)
                else:
                    # Add branches
                    for branch in cur_atom.branches:
                        frontier.append(branch.atom)

            # Connect fragment to branch
            atom.add_branches(fragment_head.branches)
            break

        return self


    def MutateBond(self):
        # Remove/add Hydrogen to increase/decrease bond
        if len(self.bonds) < 1:
            return self

        bonds = self.bonds.copy()
        while len(bonds) != 0:
            bond = bonds.pop(random.randint(0, len(bonds)-1))
            # If picked bond is single, increase
            if bond.category == "single":
                possible_mutation = []
                # Check if both connecting atoms have hydrogen
                # If yes, then randomly pick which bond increase mutation can happen
                max_increase = min(bond.atom.valence_remain(), bond.from_atom.valence_remain())
                if max_increase >= 1:
                    possible_mutation.append("double")
                if max_increase >= 2:
                    possible_mutation.append("triple")

                if not possible_mutation:
                    continue

                bond.category = random.choice(possible_mutation)

            # If picked bond is double, have a chance to either increase/decrease
            elif bond.category == "double":
                possible_mutation = ["single"]
                # Check if both connecting atoms have hydrogen
                # If yes, include triple
                if bond.atom.valence_remain() and bond.from_atom.valence_remain():
                    possible_mutation.append("triple")

                bond.category = random.choice(possible_mutation)

            # If picked bond is triple, decrease
            elif bond.category == "triple":
                # Random, either reduce to single or double bond
                bond.category = random.choice(["single", "double"])

            break

        return self


    def get_length(self):
        # Identify starting point, as much as possible should be an attachment point, else head atom
        starting_atom = self.head_atom
        for atom in self.atoms:
            if atom.value.symbol == "*":
                starting_atom = atom
                break

        # BFS search for length until end of molecule or next attachment point
        length = 0
        visited = []
        frontier = [[starting_atom]]
        found_next_attachment = False
        first = True
        while frontier:
            layer = frontier.pop(0)
            next_layer = []
            for atom in layer:
                if atom in visited:
                    continue
                visited.append(atom)

                if type(atom) != Atom:
                    continue

                if atom.value.symbol == "*" and length != 0:
                    # Reached an attachment point
                    found_next_attachment = True
                    break
                
                # Add parent
                if not atom.from_branch is None and not atom.from_branch.from_atom is None and not atom.from_branch.from_atom in visited:
                    next_layer.append(atom.from_branch.from_atom)
                # Add branches to frontier
                for branch in atom.branches:
                    if type(branch.atom) != Ring_Closure and not branch.atom in visited:
                        next_layer.append(branch.atom)
                    elif type(branch.atom) == Ring_Closure:
                        # Add other atom of Ring_Closure
                        other_atom = branch.atom.atom_a
                        if other_atom is atom:
                            other_atom = branch.atom.atom_b
                        next_layer.append(other_atom)

            if found_next_attachment:
                break

            if first:
                first = False
            else:
                length += 1

            if len(next_layer) != 0:
                frontier.append(next_layer)

        return length






