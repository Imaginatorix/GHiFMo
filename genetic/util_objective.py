import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
import psi4
from Chromosome import *
from util_clusters import get_cluster
from util_binding_affinity_predictor import *
from util_admet_selenium_extraction import automated_admet
import copy

def get_pareto_ranking(fitness_scores):
    def dominates(sol1, sol2):
        return all(x <= y for x, y in zip(sol1, sol2)) and any(x < y for x, y in zip(sol1, sol2))
    
    num_solutions = len(fitness_scores)
    pareto_levels = np.zeros(num_solutions, dtype=int)  # Initialize ranks

    for i in range(num_solutions):
        for j in range(num_solutions):
            if i != j:
                if dominates(fitness_scores[j], fitness_scores[i]):
                    pareto_levels[i] += 1  # Increase rank if solution i is dominated by solution j

    return pareto_levels


def get_scores(population, fitness_scores):
    # (Pareto Ranking, Cluster)
    scores = list(zip(get_pareto_ranking(fitness_scores), get_cluster(population)))
    return scores


def get_admet(molecules):
    molecules = [atom_to_smiles(molecule.head_atom) for molecule in molecules]
    return automated_admet(molecules)

def calculate_homo_lumo(geometry):
    # Initialize Psi4 molecule
    psi4.geometry(f"""
    {geometry}
    """)

    # Set Psi4 options
    psi4.set_options({
        "basis": "sto-3g",  # Basis set
        "scf_type": "df",  # Density fitting for speed
        "reference": "rhf"  # Restricted HF for closed-shell systems,
    })
    psi4.set_memory(5e8)
    psi4.set_num_threads(4)
    
    # Perform DFT calculation with B3LYP
    energy, wavefunction = psi4.energy("b3lyp", return_wfn=True)

    # Extract orbital energies
    eps = wavefunction.epsilon_a_subset("AO", "ALL").np

    nocc = wavefunction.nalpha()  # Number of alpha electrons (occupied orbitals)
    if nocc == len(eps):
        return 0

    # HOMO is the last occupied orbital, and LUMO is the first unoccupied orbital
    homo = eps[nocc - 1]
    lumo = eps[nocc]
    
    return (lumo-homo)*27.2114

def smiles_to_geometry(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)  # Add hydrogens
    AllChem.EmbedMolecule(mol)  # Generate initial 3D coordinates
    AllChem.UFFOptimizeMolecule(mol)  # Optimize geometry with UFF force field

    # Extract atomic symbols and 3D coordinates
    conf = mol.GetConformer()
    geometry = []
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        geometry.append(f"{atom.GetSymbol()} {pos.x:.4f} {pos.y:.4f} {pos.z:.4f}")
    return "\n".join(geometry)

def get_hlg(molecules):
    smiles_strings = [atom_to_smiles(molecule.head_atom) for molecule in molecules]
    geometries = [smiles_to_geometry(smiles_string) for smiles_string in smiles_strings]
    homo_lumo_gap = [calculate_homo_lumo(geometry) for geometry in geometries]
    return homo_lumo_gap

def get_fitness(molecules):
    # Get Admet and Get SA are merged
    admet_props = get_admet(molecules).values.tolist()

    # Predicted binding affinity (e.g., lower values are better for binding affinity)
    binding_affinity = get_binding_affinity(molecules)
    # Predicted HOMO-LUMO Gap using Quantum Chemistry
    homo_lumo_gap = get_hlg(molecules)

    for i in range(len(admet_props)):
        admet_props[i].append(binding_affinity[i])
        admet_props[i].append(homo_lumo_gap[i])

    # columns_to_extract = [
    #                       "Lipinski", "MW",               # Physiological
    #                       "pgp_inh", "pgp_sub",           # A
    #                       "logVDss",                      # D
    #                       "CYP2D6-inh", "CYP2D6-sub",     # M
    #                       "cl-plasma",                    # E
    #                       "hERG",                         # T
    #                       "Synth" # + Binding Affinity, Homo-Lumo Gap
    #                      ]

    def distance_to_ideal(value, ideal):
        if ideal == float("+inf"):
            return value
        elif ideal == float("-inf"):
            return -value
        else:
            return -abs(value - ideal)

    ideal_values = [
                    0, 350,                 # Physiological
                    0, 0,                   # A
                    10.02,                  # D
                    0, 0,                   # M
                    0,                      # E
                    0,                      # T
                    0, float("-inf"), 6     # Other
                ]

    distance = copy.deepcopy(admet_props)
    for i in range(len(distance)):
        for index in range(len(distance[i])):
            # Distance formula
            distance[i][index] = distance_to_ideal(distance[i][index], ideal_values[index])
    admet_props = [tuple(admet_prop) for admet_prop in admet_props]
    distance = [tuple(dist) for dist in distance]

    return admet_props, distance

