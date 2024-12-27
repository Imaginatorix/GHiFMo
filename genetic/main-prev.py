from Chromosome import *
from GeneticMutation import *
import copy

# S=c1nc[nH]c2nc[nH]c12
# mercaptopurine_ring = Ring_Manager()
# mercaptopurine =                                                Atom("S").add_branches([
#                                                 Branch("double", Atom("C").add_branches([
# Branch("single", mercaptopurine_ring.ring(1)), Branch("single", Atom("N").add_branches([
#                                                 Branch("double", Atom("C").add_branches([
#                                                 Branch("single", Atom("N").add_branches([
#                                                 Branch("single", Atom("C").add_branches([
# Branch("double", mercaptopurine_ring.ring(2)), Branch("single", Atom("N").add_branches([
#                                                 Branch("double", Atom("C").add_branches([
#                                                 Branch("single", Atom("N").add_branches([
#                                                 Branch("single", Atom("C").add_branches([
# Branch("single", mercaptopurine_ring.ring(1)),
# Branch("double", mercaptopurine_ring.ring(2))
# ]))
# ]))
# ]))
# ]))
# ]))
# ]))
# ]))
# ]))
# ]))
# ])

# mercaptopurine.show()

# C_string_ring = Ring_Manager()
# C_string =       Atom("C").add_branches([
# Branch("single", Atom("N").add_branches([
# Branch("single", Atom("O").add_branches([
# Branch("single", Atom("Cl").add_branches([
# ]))
# ]))
# ]))
# ])

# C_string_ring_copy = Ring_Manager
# C_string_copy = copy.deepcopy(C_string)

# C_string.show()

# a = Mutate(C_string, C_string_ring)
# b = Mutate(C_string_copy, C_string_ring_copy)

# print("Close Ring")
# a.CloseRing()
# a.recount()
# a.head_atom.show()

# print()
# print("Mutate Bond")
# a.MutateBond()
# a.recount()
# a.head_atom.show()

# print()
# print("Mutate Atom")
# a.MutateAtom()
# a.recount()
# a.head_atom.show()

# print()
# print("Exchange Vicinal")
# a.ExchangeVicinal()
# a.recount()
# a.head_atom.show()

# print()
# print("Open Ring")
# a.OpenRing()
# a.recount()
# a.head_atom.show()

# print()
# print("(a)")
# a.head_atom.show()
# print("(b)")
# b.head_atom.show()
# print("Crossover")
# a.CrossOver(b)
# print("(a)")
# a.head_atom.show()
# print("(b)")
# b.head_atom.show()
# a.recount()

# print()
# print("MutateBranch")
# a.MutateBranch()
# a.head_atom.show()

# print()
# print("MutateFragment")
# a.MutateFragment()
# a.recount()
# a.head_atom.show()


# atom_from_smiles("C1=CC=CC=C1")[0].show()
# atom_from_smiles("[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O")[0].show()
# atom_from_smiles("*C(=O)Cc1c(C)n(*)c2ccc(OC)cc12")[0].show()

# print(atom_to_smiles(atom_from_smiles("*C(=O)Cc1c(C)n(*)c2ccc(OC)cc12")[0]))


# a.randomMutate()

a = atom_from_smiles("S=c1nc[nH]c2nc[nH]c12")
b = atom_from_smiles("Nc2nc(=S)c1[nH]cnc1[nH]2")[0]

a_mutate = Mutate(a[0], a[1])
a_mutate.head_atom.show()
a_mutate.ScaffoldHop()[0].show()

