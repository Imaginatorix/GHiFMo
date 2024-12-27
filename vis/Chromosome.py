from mendeleev import element
from rdkit import Chem

class Attachment_Point():
    def __init__(self):
        self.symbol = "*"
        self.period = 4

    def nvalence(self):
        return 4

class Atom():
    def __init__(self, value, branches=None, from_branch=None):
        if value == "*":
            self.value = Attachment_Point()
        else:
            self.value = element(value)

        self.from_branch = from_branch

        if branches is None:
            branches = []
        self.branches = branches

        self.valid_valence(error=True)

    def valid_valence(self, error=False):
        if self.valence_used() > self.atom_valence():
            if error:
                raise ValueError("Invalid molecule configuration. Too many valence electrons used.")
            return False
        return True

    def valence_remain(self):
        return self.atom_valence() - self.valence_used()

    def valence_used(self):
        # Branch connecting to it
        valence = 0
        if self.from_branch is not None:
            valence = self.from_branch.get_valence()

        # Extending branches
        for branch in self.branches:
            valence += branch.get_valence()

        return valence

    def atom_valence(self):
        # def allowed_valence(period):
        #     if period == 1:
        #         return 2
        #     elif period == 2:
        #         return 8
        #     elif period == 3:
        #         return 18
        #     else:
        #         return 32

        return min(self.value.nvalence(), 8-self.value.nvalence())

    def add_branches(self, branches):
        # Add branches to self.branch
        self.branches += branches
        for branch in branches:
            # Add from_atom in each branch
            branch.from_atom = self

        # Connect ring to atom
        for branch in branches:
            if type(branch.atom) == Ring_Closure:
                branch.atom.connect(self)

        self.valid_valence(error=True)
        return self

    def show(self, loc=None):
        # Print self
        if loc is None:
            loc = ["H"]
            print(self.value.symbol)

        # Print Connections
        connections = []
        for branch in self.branches:
            branch2str = {"single": "-", "double": "=", "triple": "#"}
            if type(branch.atom) == Atom:
                connections.append(branch2str[branch.category]+branch.atom.value.symbol)
            elif type(branch.atom) == Ring_Closure:
                connections.append(branch2str[branch.category]+branch.atom.symbol)

        if len(connections) != 0:
            print(f"{loc} {' '.join(connections)}")

        # Print next level
        for i, branch in enumerate(self.branches):
            if type(branch.atom) == Atom:
                branch.atom.show(loc+[i+1])

    @staticmethod
    def is_valid(molecule, ring_manager):
        # Check valence
        if not molecule.valid_valence():
            return False
        # Check rings
        for ring in ring_manager.rings.values():
            # Unmatched rings
            if not ring.is_closed():
                return False
            # Different ring bonds
            if ring.from_branch_a.category != ring.from_branch_b.category:
                return False
            # Valid ring size
            ring.render_members()
            if not ring.valid_ring_size():
                return False
        return True


class Branch():
    def __init__(self, category, atom, from_atom=None):
        self.category = category
        if self.category not in ["single", "double", "triple"]:
            raise ValueError("Branch must either be 'single', 'double', or 'triple'")

        if type(atom) == Ring_Closure:
            if atom.from_branch_a:
                if not atom.from_branch_a.category is self.category:
                    raise ValueError("Ring Closure bond must be the same")
                atom.from_branch_b = self
                self.from_atom = from_atom
                self.atom = atom
            else:
                atom.from_branch_a = self
                self.from_atom = from_atom
                self.atom = atom
        else:
            atom.from_branch = self
            self.from_atom = from_atom
            self.atom = atom

    def get_valence(self):
        valence = {
            "single": 1,
            "double": 2,
            "triple": 3
            # "dative": 2 # 2 or 0
        }
        return valence[self.category]


class Ring_Manager():
    def __init__(self):
        self.highest_number = None # Always increasing, how many rings were formed in this specific molecule?
        self.rings = {}

    def add_ring(self, identity=None, ring_closure_object=None):
        if identity is None:
            identity = self.last_untaken_ring()

        if ring_closure_object is None:
            self.rings[identity] = Ring_Closure(identity)
        else:
            ring_closure_object.identity = identity
            self.rings[identity] = ring_closure_object
        self.highest_number += 1

        return self.rings[identity]

    def last_untaken_ring(self):
        if self.highest_number is None:
            self.highest_number = len(self.rings)+1
        return self.highest_number

    def ring(self, identity):
        if not identity in self.rings:
            self.add_ring(identity)
        return self.rings[identity]

    def part_of_ring(self, atom, show_groups=False):
        if show_groups:
            ring_groups = []
            for ring in self.rings.values():
                if ring.is_closed() and atom in ring.render_members():
                    ring_groups.append(ring)
            return ring_groups

        for ring in self.rings.values():
            if ring.is_closed() and atom in ring.render_members():
                return True
        return False

    def remove_ring(self, identity):
        self.rings[identity].remove()
        del self.rings[identity]

class Ring_Closure():
    def __init__(self, identity):
        self.identity = identity
        self.symbol = f"R[{identity}]"
        self.from_branch_a = None
        self.from_branch_b = None
        self.atom_a = None
        self.atom_b = None
        self.members = None
        self.top_heirarchy = None

    def valence_remain(self):
        return min(self.atom_a.valence_remain(), self.atom_b.valence_remain())

    def connect(self, atom):
        if self.atom_a is None:
            self.atom_a = atom
            return
        if self.atom_b is None:
            self.atom_b = atom
            # Both are already connected, check if it's neither a 1- or 2-membered ring
            self.render_members()
            self.valid_ring_size(error=True)
            return

        raise ValueError("Ring already fully connected")

    def remove(self):
        # Remove connection from previous atom to branch
        self.atom_a.branches.remove(self.from_branch_a)
        self.from_branch_a.from_atom = None
        self.atom_b.branches.remove(self.from_branch_b)
        self.from_branch_b.from_atom = None
        self.atom_a = None
        self.atom_b = None
        # Remove connection from branch to this
        self.from_branch_a.atom = None
        self.from_branch_a = None
        self.from_branch_b.atom = None
        self.from_branch_b = None

    def is_closed(self):
        return self.atom_a and self.atom_b

    def render_members(self):
        # Look at first ring closure and then BFS until the other ring closure is found
        if not self.is_closed():
            raise ValueError("Cannot render members of unclosed ring closure")

        # BFS
        prev = {self.atom_a: -1}
        frontier = [self.atom_a]
        while frontier:
            atom = frontier.pop(0)
            if atom is self.atom_b:
                break

            for branch in atom.branches:
                if not branch.atom in prev and type(branch.atom) != Ring_Closure:
                    prev[branch.atom] = atom
                    frontier.append(branch.atom)

            if atom.from_branch:
                if not atom.from_branch.from_atom in prev:
                    prev[atom.from_branch.from_atom] = atom
                    frontier.append(atom.from_branch.from_atom)

        # Backtrack
        previous_atom = self.atom_b
        self.top_hierarchy_index = 0
        self.members = []
        while previous_atom != -1:
            self.members.append(previous_atom)
            if previous_atom.from_branch and prev[previous_atom] is previous_atom.from_branch.from_atom:
                # It's a parent
                self.top_hierarchy_index = len(self.members)
            previous_atom = prev[previous_atom]
        return self.members

    def valid_ring_size(self, error=False):
        # 1-ring
        if len(self.members) == 1:
            if error:
                raise ValueError("Cannot create ring closure with a 1-membered ring")
            return False
        # 2-ring
        elif len(self.members) == 2:
            if error:
                raise ValueError("Cannot create ring closure with a 2-membered ring")
            return False

        return True
        # # If they ends are the same, it is 1-membered
        # if self.atom_a is self.atom_b:
        #     if error:
        #         raise ValueError("Cannot create ring closure with a 1-membered ring")
        #     return False
        # # If the they are connected in some other way (aside from the ring closures), it is 2-membered
        # # Check all branches of atom_a (except for the ring closure)
        # for branch in self.atom_a.branches:
        #     if type(branch.atom) != Ring_Closure and branch.atom is self.atom_b:
        #         if error:
        #             raise ValueError("Cannot create ring closure with a 2-membered ring")
        #         return False
        # # Check all branches of atom_b (except for the ring closure)
        # for branch in self.atom_b.branches:
        #     if type(branch.atom) != Ring_Closure and branch.atom is self.atom_a:
        #         if error:
        #             raise ValueError("Cannot create ring closure with a 2-membered ring")
        #         return False
        # return True

def atom_to_smiles(node, is_branch=False, ring_tags=None):
    branch2str = {
        "single": "",
        "double": "=",
        "triple": "#"
    }

    if ring_tags is None:
        ring_tags = {}

    smiles = ""
    if is_branch:
        smiles += branch2str[node.category]
        node = node.atom
    smiles += node.value.symbol

    # Start with attaching all rings
    branches = node.branches.copy()
    for branch in branches.copy():
        if type(branch.atom) == Ring_Closure:
            # Add bond
            smiles += branch2str[branch.category]
            # Add ring closure if it doesn't exist
            if not branch.atom.identity in ring_tags:
                ring_tags[branch.atom.identity] = len(ring_tags)+1

            ring_tag = str(ring_tags[branch.atom.identity])
            # If >= double digit
            if len(ring_tag) >= 2:
                ring_tag = "$"+ring_tag

            smiles += ring_tag
            branches.remove(branch)

    # All excepted first is a substituent
    for i in range(len(branches)-1, 0, -1):
        smiles += f"({atom_to_smiles(branches[i], True, ring_tags)})"
    # Add first as a main backbone
    if branches:
        smiles += f"{atom_to_smiles(branches[0], True, ring_tags)}"

    if not is_branch:
        # Turn to canonical smiles through rdkit
        smiles = Chem.MolToSmiles(Chem.MolFromSmiles(smiles))

    return smiles

def atom_from_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    # Kekulize to have non-aromatic bonds
    Chem.KekulizeIfPossible(mol)
    ring_manager = Ring_Manager()
    atoms = [Atom(atom.GetSymbol()) for atom in mol.GetAtoms()]
    connections = Chem.GetAdjacencyMatrix(mol, useBO=True)

    int2bond = {
        1: "single",
        2: "double",
        3: "triple"
    }

    # Length of adjacency matrix
    n = len(connections)

    # For ring detection
    branched_to = []
    for iy in range(n):
        for ix in range(iy+1, n):
            if connections[iy][ix] == 0:
                continue

            # If ring
            if atoms[ix] in branched_to:
                new_ring = ring_manager.add_ring()
                atoms[iy].add_branches([Branch(int2bond[connections[iy][ix]], new_ring)])
                atoms[ix].add_branches([Branch(int2bond[connections[iy][ix]], new_ring)])
            else:
                atoms[iy].add_branches([Branch(int2bond[connections[iy][ix]], atoms[ix])])
                branched_to.append(atoms[ix])

    return (atoms[0], ring_manager)


if __name__ == "__main__":
    # S=c1nc[nH]c2nc[nH]c12
    mercaptopurine_ring = Ring_Manager()
    mercaptopurine =                                                Atom("S").add_branches([
                                                   Branch("double", Atom("C").add_branches([
    Branch("single", mercaptopurine_ring.ring(1)), Branch("single", Atom("N").add_branches([
                                                   Branch("double", Atom("C").add_branches([
                                                   Branch("single", Atom("N").add_branches([
                                                   Branch("single", Atom("C").add_branches([
    Branch("double", mercaptopurine_ring.ring(2)), Branch("single", Atom("N").add_branches([
                                                   Branch("double", Atom("C").add_branches([
                                                   Branch("single", Atom("N").add_branches([
                                                   Branch("single", Atom("C").add_branches([
    Branch("single", mercaptopurine_ring.ring(1)),
    Branch("double", mercaptopurine_ring.ring(2))
    ]))
    ]))
    ]))
    ]))
    ]))
    ]))
    ]))
    ]))
    ]))
    ])

    mercaptopurine.show()

    # print(Atom.is_valid(mercaptopurine, mercaptopurine_ring))




