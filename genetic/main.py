import os
from Chromosome import *
from GeneticMutation import *
from util_binding_affinity_predictor import TanimotoKernel, smiles_to_fingerprint
from util_mutation import select_parents, mutation, crossover
from util_objective import get_fitness, get_scores
import pickle

POPULATION_SIZE = 150
NUM_PARENTS = POPULATION_SIZE
MUT_RATE = 1
CROSS_RATE = 1
GENERATIONS = 1000

def init_population():
    # Initialize a population of random molecular structures
    # population_chromosomes = [atom_from_smiles("Nc1nc(-c2cccs2)c(NC(=O)c2ccccc2)s1"), atom_from_smiles("CC(C)OC(=O)CC(=O)CSc1nc2c(cc1C#N)CCC2"),
    #         atom_from_smiles("N#CC(C#N)=c1ccc2c(c1)NC1(CCCCC1)N=2"), atom_from_smiles("CCc1n[nH]c2c1/C(=N\O)CC(c1ccccc1)C2")]

    # molecules = [mercaptopurine, tioguanine, methotrexate, cytarabine]
    molecules = ["S=c1nc[nH]c2nc[nH]c12", "Nc2nc(=S)c1[nH]cnc1[nH]2", "O=C([C@H](CCC(O)=O)NC(C1=CC=C(N(CC2=CN=C(N=C(N)N=C3N)C3=N2)C)C=C1)=O)O", "O=C1/N=C(/N)\C=C/N1[C@@H]2O[C@@H]([C@@H](O)[C@@H]2O)CO"]
    population_chromosomes = [atom_from_smiles(smiles_string) for smiles_string in molecules]
    return [Mutate(head_atom, ring_manager) for head_atom, ring_manager in population_chromosomes]


def population_from_smiles(smiles):
    population_chromosomes = [atom_from_smiles(smiles_string) for smiles_string in smiles]
    return [Mutate(head_atom, ring_manager) for head_atom, ring_manager in population_chromosomes]


# Genetic Algorithm for molecular structures
def genetic_algorithm(generations, mut_rate, cross_rate, num_parents):
    # Create log.pkl if it doesn't exist
    log_path = "./history/log.pkl"
    if not os.path.exists(log_path):
        with open(log_path, "wb") as f:
            p = pickle.Pickler(f)
            p.dump([])
    
    # Read log.pkl
    with open(log_path, "rb") as f:
        up = pickle.Unpickler(f)
        history = up.load()

    if len(history) == 0:
        # Initialize population and archive
        population = init_population()
        pareto_archive = []
    else:
        population = history[-1]["population_molecules"]
        pareto_archive = history[-1]["pareto_archive"]

    for generation in range(len(history), generations):
        print("History")
        print(history)
        
        print("Running generation", generation)
        # Get scores
        print("Getting fitness scores...")
        if len(history) >= 1:
            fitness_scores = history[-1]["scores"]
        else:
            fitness_scores = get_fitness(population)

        print("Ranking molecules:")
        scores = get_scores(population, fitness_scores)
        print(scores)

        # Parent Selection [Tournament Selection]
        print("Selecting parents...")
        parents = select_parents(population, scores, num_parents)

        # Mutate
        print("Mutating...")
        offsprings_mutate = mutation(population, mut_rate)
        # Crossover
        print("Breeding...")
        offsprings_crossover = crossover(parents, cross_rate)

        # New population
        new_population_unfiltered = list(set(population + offsprings_mutate + offsprings_crossover + pareto_archive))

        # Filter all invalid through rdkit
        print("New population :OO")
        new_population = []
        unique_smiles = []
        for molecule in new_population_unfiltered:
            try:
                smiles_string = atom_to_smiles(molecule.head_atom)
                if smiles_string in unique_smiles:
                    continue
                unique_smiles.append(smiles_string)

                mol = Chem.MolFromSmiles(smiles_string)
                mol = Chem.AddHs(mol)  # Add hydrogens
                AllChem.EmbedMolecule(mol)  # Generate initial 3D coordinates
                AllChem.UFFOptimizeMolecule(mol)  # Optimize geometry with UFF force field
                fp = smiles_to_fingerprint(smiles_string)
                if fp is None:
                    continue

                new_population.append(molecule)
            except Exception as e:
                print(e)
                pass

        print(len(new_population))

        # Get scores
        print("Filtering new population...")
        print("> Getting fitness scores...")
        new_fitness_scores = get_fitness(new_population)
        print("> Ranking molecules:")
        new_scores = get_scores(new_population, new_fitness_scores)

        pareto_archive = []
        for i in range(len(population)):
            if scores[i][1] == 0:
                pareto_archive.append(population[i])

        # Look-up table {Mutate: fitness}
        new_population_fitness_look_up = dict(zip(new_population, new_fitness_scores))
        # Look-up table {Mutate: (rank, cluster)}
        new_population_look_up = dict(zip(new_population, new_scores))
        # Create a dictionary {cluster: [Mutate]}
        new_population_clusters = {}
        # Create a dictionary {cluster: average rank}
        new_population_cluster_rank = {}

        # Populate dictionaries
        for individual, score in zip(new_population, new_scores):
            if not score[1] in new_population_clusters:
                new_population_clusters[score[1]] = []
            if not score[1] in new_population_cluster_rank:
                new_population_cluster_rank[score[1]] = []

            new_population_clusters[score[1]].append(individual)
            new_population_cluster_rank[score[1]].append(score[0])

        # Take sort [Mutate] by rank
        for cluster in new_population_clusters:
            new_population_clusters[cluster] = sorted(new_population_clusters[cluster], key=lambda x: new_population_look_up[x][0])

        # Take average of cluster_rank
        for cluster in new_population_cluster_rank:
            ranks = new_population_cluster_rank[cluster]
            new_population_cluster_rank[cluster] = sum(ranks)/len(ranks)

        # Sort clusters, so higher ranking (lower numbers) would be picked first
        sorted_clusters = sorted(list(new_population_cluster_rank.keys()), key=lambda x: new_population_cluster_rank[x])

        # Population picking until POPULATION_SIZE
        print("> Filtering...")
        new_population_ = []
        new_population_smiles = []
        new_population_fitness = []
        loops = 0
        while len(new_population_) < POPULATION_SIZE:
            error = 0
            for cluster in sorted_clusters:
                if loops >= len(new_population_clusters[cluster]):
                    error += 1
                    continue

                individual = new_population_clusters[cluster][loops]

                canon_identity = atom_to_smiles(individual.head_atom)
                if not canon_identity in new_population_smiles:
                    new_population_.append(individual)
                    new_population_smiles.append(canon_identity)
                    new_population_fitness.append(new_population_fitness_look_up[individual])

                if len(new_population_) >= POPULATION_SIZE:
                    break

            if error >= len(sorted_clusters):
                break
            loops += 1

        # Set new population to population
        print("Next gen population:")
        population = new_population_.copy()
        print(len(population))

        # Print best binding affinity
        print(f"Generation {generation}: {new_population_smiles}")

        # Save history to log.pkl
        print("Logging to history...")
        history.append({
            "population_molecules": population,
            "population_smiles": new_population_smiles,
            "pareto_archive": [atom_to_smiles(individual.head_atom) for individual in pareto_archive],
            "scores": new_population_fitness
        })
        with open(log_path, "wb") as f:
            p = pickle.Pickler(f)
            p.dump(history)

    print(population)


if __name__ == "__main__":
    genetic_algorithm(GENERATIONS, MUT_RATE, CROSS_RATE, NUM_PARENTS)




