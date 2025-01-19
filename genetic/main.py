import os
from Chromosome import *
from GeneticMutation import *
from util_binding_affinity_predictor import TanimotoKernel, smiles_to_fingerprint, process_smiles
from util_mutation import select_parents, mutation, crossover
from util_objective import get_fitness, get_scores
import pickle
import numpy as np
from sklearn.preprocessing import StandardScaler
import joblib

print("IMPORTS COMPLETE")

POPULATION_SIZE = 150
NUM_PARENTS = POPULATION_SIZE
MUT_RATE = 1
CROSS_RATE = 1
GENERATIONS = 1000

def init_population():
    molecules = ["O=C2Nc1c(ccnc1N(c3ncccc23)C4CC4)C", "O=C1Nc2ccc(Cl)cc2[C@@](C#CC2CC2)(C(F)(F)F)O1"]
    population_chromosomes = [atom_from_smiles(smiles_string) for smiles_string in molecules]
    return [Mutate(head_atom, ring_manager) for head_atom, ring_manager in population_chromosomes]

def population_from_smiles(smiles):
    population_chromosomes = [atom_from_smiles(smiles_string) for smiles_string in smiles]
    return [Mutate(head_atom, ring_manager) for head_atom, ring_manager in population_chromosomes]

def evaluate_fitness_with_binding_affinity(population, model, scaler):
    smiles_list = [atom_to_smiles(individual.head_atom) for individual in population]
    valid_smiles, fingerprints = process_smiles(smiles_list)
    
    fingerprints_scaled = scaler.transform(fingerprints)
   
    predicted_affinities = model.predict(fingerprints_scaled)
    
    return predicted_affinities

def genetic_algorithm(generations, mut_rate, cross_rate, num_parents, model, scaler):
    log_path = "./history/log.pkl"
    if not os.path.exists(log_path):
        with open(log_path, "wb") as f:
            p = pickle.Pickler(f)
            p.dump([])

    with open(log_path, "rb") as f:
        up = pickle.Unpickler(f)
        history = up.load()

    if len(history) == 0:
        population = init_population()
        pareto_archive = []
    else:
        population = history[-1]["population_molecules"]
        pareto_archive = history[-1]["pareto_archive"]

    for generation in range(len(history), generations):
        print("Running generation", generation)
        
        # Get fitness using binding affinity prediction
        predicted_affinities = evaluate_fitness_with_binding_affinity(population, model, scaler)
        fitness_scores = predicted_affinities
        
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
        print("Filtering new population...")
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

        # Evaluate fitness scores for the new population
        predicted_affinities = evaluate_fitness_with_binding_affinity(new_population, model, scaler)
        print("> Ranking new molecules:")
        new_scores = get_scores(new_population, predicted_affinities)

        pareto_archive = []
        for i in range(len(population)):
            if scores[i][1] == 0:
                pareto_archive.append(population[i])

        # Population clustering and ranking
        new_population_fitness_look_up = dict(zip(new_population, predicted_affinities))
        new_population_look_up = dict(zip(new_population, new_scores))
        
        new_population_clusters = {}
        new_population_cluster_rank = {}

        for individual, score in zip(new_population, new_scores):
            if not score[1] in new_population_clusters:
                new_population_clusters[score[1]] = []
            if not score[1] in new_population_cluster_rank:
                new_population_cluster_rank[score[1]] = []

            new_population_clusters[score[1]].append(individual)
            new_population_cluster_rank[score[1]].append(score[0])

        for cluster in new_population_clusters:
            new_population_clusters[cluster] = sorted(new_population_clusters[cluster], key=lambda x: new_population_look_up[x][0])

        for cluster in new_population_cluster_rank:
            ranks = new_population_cluster_rank[cluster]
            new_population_cluster_rank[cluster] = sum(ranks) / len(ranks)

        sorted_clusters = sorted(list(new_population_cluster_rank.keys()), key=lambda x: new_population_cluster_rank[x])

        new_population_ = []
        new_population_smiles = []
        new_population_fitness = []
        new_population_fitness_scores = []
        loops = 0
        while len(new_population_) < POPULATION_SIZE:
            error = 0
            for cluster in sorted_clusters:
                if loops >= len(new_population_clusters[cluster]):
                    error += 1
                    continue

                individual = new_population_clusters[cluster][loops]
                canon_identity = atom_to_smiles(individual.head_atom)
                if canon_identity not in new_population_smiles:
                    new_population_.append(individual)
                    new_population_smiles.append(canon_identity)
                    new_population_fitness.append(new_population_fitness_look_up[individual])
                    new_population_fitness_scores.append(new_population_look_up[individual])

                if len(new_population_) >= POPULATION_SIZE:
                    break

            if error >= len(sorted_clusters):
                break
            loops += 1

        # Set new population to population
        print("Next gen population:")
        population = new_population_.copy()
        print(len(population))

        # Save history to log.pkl
        print("Logging to history...")
        history.append({
            "population_molecules": population,
            "population_smiles": new_population_smiles,
            "pareto_archive": [atom_to_smiles(individual.head_atom) for individual in pareto_archive],
            "fitness": new_population_fitness,
            "fitness_scores": new_population_fitness_scores
        })
        with open(log_path, "wb") as f:
            p = pickle.Pickler(f)
            p.dump(history)

    print(population)

if __name__ == "__main__":
 
    model_filepath = ".genetic\trained_model.pkl"
    scaler_filepath = ".genetic\scaler.pkl"
    model = joblib.load(model_filepath)
    scaler = joblib.load(scaler_filepath)

    genetic_algorithm(GENERATIONS, MUT_RATE, CROSS_RATE, NUM_PARENTS, model, scaler)




