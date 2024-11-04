from Chromosome import *
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolAlign
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem.Draw import MolDrawing, DrawingOptions
import py3Dmol
import numpy as np
import random
import json
from urllib.request import urlopen
from bs4 import BeautifulSoup

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


    def randomMutate(self):
        mutation_type = random.choices(
                        ["openRing", "closeRing", "vicinal", "atom", "fragment", "branch", "bond", "scaffold"], weights=
                        [      0.08,        0.08,      0.18,   0.23,       0.05,     0.02,   0.18,       0.18])[0]
        
        i = 0
        modified = False
        while i < 5:
            try:
                if mutation_type == "scaffold":
                    self.head_atom, self.ring_manager = self.ScaffoldHop()
                    modified = True
                    break
                elif mutation_type == "openRing" and self.OpenRing():
                    modified = True
                    break
                elif mutation_type == "closeRing" and self.CloseRing():
                    modified = True
                    break
                elif mutation_type == "vicinal" and self.ExchangeVicinal():
                    modified = True
                    break
                elif mutation_type == "atom" and self.MutateAtom():
                    modified = True
                    break
                elif mutation_type == "fragment" and self.MutateFragment():
                    modified = True
                    break
                elif mutation_type == "branch" and self.MutateBranch():
                    modified = True
                    break
                elif mutation_type == "bond" and self.MutateBond():
                    modified = True
                    break
            except Exception as e:
                print(e)
            i += 1

        self.recount()
        return modified


    def ScaffoldHop(self):
        smiles = atom_to_smiles(self.head_atom)
        start_mol = Chem.MolFromSmiles(smiles)
        # Your initial scaffold
        scaffold = MurckoScaffold.GetScaffoldForMol(start_mol)
        scaffold_with_R = Chem.ReplaceSidechains(start_mol, scaffold)
        substituents_with_R = Chem.ReplaceCore(start_mol, scaffold)
        # Gather scaffold hop suggestions from website
        # https://peter-ertl.com/cgi/skeys/skeys.py?smiles=c1ccc(-c2ccccc2OC2CC2)nc1&coloring=tclass
        page = urlopen(f"https://peter-ertl.com/cgi/skeys/skeys.py?smiles={Chem.MolToSmiles(scaffold)}&coloring=tclass").read()
        script = BeautifulSoup(page, features="html.parser").findAll('script')[2].decode()
        data = script.split("var hits=")[1].split(";")[0].replace(",}", "}")
        data = json.loads(data)
        # Scaffold suggestion
        suggestion = Chem.MolFromSmiles(random.choice(data)["smiles"])

        # Original attachment point indexes
        scaffold_attachment_points = []
        for atom in scaffold_with_R.GetAtoms():
            if atom.GetSymbol() == "*":
                scaffold_attachment_points.append(atom.GetIdx())

        # Align molecules
        NUM_CONFS = 250
        p = AllChem.ETKDGv2()
        p.verbose = True
        scaffold_ = Chem.AddHs(scaffold)
        suggestion_ = Chem.AddHs(suggestion)
        AllChem.EmbedMultipleConfs(scaffold_, NUM_CONFS, p)
        AllChem.EmbedMultipleConfs(suggestion_, NUM_CONFS, p)
        mmff_params = [AllChem.MMFFGetMoleculeProperties(mol) for mol in [scaffold_, suggestion_]]
        mmff_ref_param = mmff_params[0]
        mmff_prob_params = mmff_params[1]
        # Align scaffold B to scaffold A using O3A
        tempscore = []
        for cid in range(NUM_CONFS):
            alignment = rdMolAlign.GetO3A(suggestion_, scaffold_, mmff_prob_params, mmff_ref_param, cid, 0)
            alignment.Align()
            tempscore.append(alignment.Score())
        best = np.argmax(tempscore)

        # Find atoms closest to attachment points
        # Get the 3D coordinates of the attachment points in scaffold A
        conf_A = scaffold_.GetConformer(0)
        coords_A = np.array([conf_A.GetAtomPosition(pt) for pt in scaffold_attachment_points])
        # Get the 3D coordinates of all atoms in scaffold B
        conf_B = suggestion_.GetConformer(int(best))
        coords_B = np.array([conf_B.GetAtomPosition(i) for i in range(suggestion_.GetNumAtoms()) if suggestion_.GetAtoms()[i].GetSymbol() != "H"])
        # Find the closest atoms in scaffold B to the attachment points in scaffold A
        suggestion_attachment_points = []
        for coord_A in coords_A:
            distances = np.linalg.norm(coords_B - coord_A, axis=1)
            suggestion_attachment_points.append(np.argmin(distances))

        # Add attachment points to suggestion
        suggestion = Chem.RemoveHs(suggestion)
        # Create an editable version of the molecule
        editable_mol = Chem.EditableMol(suggestion)
        attach_points = []
        # Add the generic atom (use * for generic atom) and connect it
        for i, idx in enumerate(suggestion_attachment_points):
            attachment_point = editable_mol.AddAtom(Chem.Atom(0))  # Dummy atom (*)
            attach_points.append(attachment_point)
            editable_mol.AddBond(int(idx), attachment_point, order=Chem.BondType.SINGLE)
        # Create the new molecule
        suggestion_with_R = editable_mol.GetMol()
        for i, point in enumerate(attach_points):
            suggestion_with_R.GetAtomWithIdx(point).SetIsotope(i+1)
        Chem.SanitizeMol(suggestion_with_R)

        # Connect the attachment points together to create final molecule
        # Create an editable version of the molecule
        mol_to_connect = Chem.RemoveHs(Chem.CombineMols(suggestion_with_R, substituents_with_R))
        R_groups_to_connect = {}
        editable_final = Chem.EditableMol(mol_to_connect)
        # Indexes to connect to respective attachment points
        for atom in mol_to_connect.GetAtoms():
            if atom.GetSymbol() == "*":
                if not atom.GetIsotope() in R_groups_to_connect:
                    R_groups_to_connect[atom.GetIsotope()] = []
                R_groups_to_connect[atom.GetIsotope()].append(atom.GetIdx())
        # Connect indexes
        for attachment_points in R_groups_to_connect.values():
            # Assume length is two
            if len(attachment_points) != 2:
                continue
            a_ = mol_to_connect.GetAtomWithIdx(attachment_points[0]).GetNeighbors()
            b_ = mol_to_connect.GetAtomWithIdx(attachment_points[1]).GetNeighbors()
            for a in a_:
                for b in b_:
                    editable_final.AddBond(a.GetIdx(), b.GetIdx(), order=Chem.BondType.SINGLE)
        # Remove atoms
        while True:
            broken = False
            for atom in editable_final.GetMol().GetAtoms():
                if atom.GetSymbol() == "*":
                    editable_final.RemoveAtom(atom.GetIdx())
                    broken = True
                    break
            if not broken:
                break
        # Create the new molecule
        final = editable_final.GetMol()

        return atom_from_smiles(Chem.MolFromSmiles(final))[0]


    def OpenRing(self):
        # Remove the ring bond and change it to H
        # Check if it needs re-rooting
        rings = list(self.ring_manager.rings.values())
        if len(rings) < 1:
            return True

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
        return True

    def CloseRing(self):
        # Randomly pick two atoms (with Hydrogens) and connect
        # Chosen atom must not be vicinal (adjacent) or geminal (same)
        if len(self.atoms) <= 2:
            return True

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
        return True


    def ExchangeVicinal(self):
        # The purpose of this is to extend or contract rings
        # Exchange hydrogen with a singly bond substructure at an adjacent position
        if len(self.atoms) <= 2:
            return True

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

        return True


    def CrossOver(self, mutate_head_atom):
        # Randomly pick two non-ring fragments of molecules and exchange
        if len(self.bonds) < 1 or len(mutate_head_atom.bonds) < 1:
            return True

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
                return False

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
        return True


    def MutateAtom(self):
        # Pick an atom and change it to a new one with the same valence (or +-1 still making the molecule valid)
        # All supported atoms so far in the mutation
        if len(self.atoms) < 1:
            return True

        possible_atoms = {
            1: ["F", "Cl", "Br", "I"],
            2: ["O", "S", "Se", "Hg"],
            3: ["B", "N", "Al", "P", "As"],
            4: ["C", "Si"]
        }

        possible_atoms_weight = {
            1: [3, 5, 3, 3],
            2: [5, 3, 1, 1],
            3: [3, 5, 1, 3, 1],
            4: [5, 1]
        }

        atom = self.atoms.pop(random.randint(0, len(self.atoms)-1))
        atom_valence = atom.atom_valence()
        possible_mutations = possible_atoms[atom_valence]
        possible_mutations_weights = possible_atoms_weight[atom_valence]
        if atom_valence < 4:
            possible_mutations += possible_atoms[atom_valence+1]
            possible_mutations_weights += possible_atoms_weight[atom_valence+1]
        if atom_valence > 1 and atom.valence_remain() > 0:
            possible_mutations += possible_atoms[atom_valence-1]
            possible_mutations_weights += possible_atoms_weight[atom_valence-1]

        atom.value = element(random.choices(possible_mutations, weights=possible_mutations_weights)[0])
        return True


    def MutateFragment(self):
        # Pick a fragment and change it to a new one with the same length in the database
        if len(self.atoms) < 4:
            return True

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
            return False

        fragment_weights = [all_fragments_frequency[fragment] for fragment in all_fragments_length[str(chosen_length)]]
        chosen_fragment = random.choices(all_fragments_length[str(chosen_length)], weights=fragment_weights)[0]

        # TODO: Fix valence issues
        try:
            chosen_fragment_hier = atom_from_smiles(chosen_fragment)
        except Exception:
            return False

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
                return False

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
                    return False
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

        return True


    def MutateBranch(self):
        # Replace branch with a fragment in the library
        if len(self.atoms) < 1:
            return True

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

        return True


    def MutateBond(self):
        # Remove/add Hydrogen to increase/decrease bond
        if len(self.bonds) < 1:
            return True

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

        return True


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






