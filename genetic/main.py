import os
import time
import glob
import pandas as pd
import numpy as np
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.keys import Keys
import random
import copy
from rdkit import DataStructs
from rdkit.Chem import AllChem
from sklearn.metrics import silhouette_score, pairwise_distances
from Chromosome import *
from GeneticMutation import *
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from binding_affinity_predictor import *
from admet_selenium_extraction import automated_admet
import json

POPULATION_SIZE = 1000
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


def select_parents(population, scores, num_parents):
    # Select a parent based on scores [Tournament Selection]
    clustered = {}
    for i in range(len(population)):
        if not scores[i] in clustered:
            clustered[scores[i]] = []
        clustered[scores[i]].append((population[i], scores[i]))

    parents = []
    while len(parents) != num_parents:
        for cluster, individual_scores in clustered.items():
            parent1, rank1 = random.choice(individual_scores)
            parent2, rank2 = random.choice(individual_scores)

            # Better gets selected
            if rank1 > rank2:
                parents.append(parent2)
            else:
                parents.append(parent1)

            if len(parents) == num_parents:
                break
    return parents


def crossover(parents, cross_rate):
    if len(parents) % 2 != 0:
        parents = parents[:-1]
        if len(parents) == 0:
            return []

    offsprings = []

    random.shuffle(parents)
    parents1 = parents[:len(parents) // 2]
    parents2 = parents[len(parents) // 2:]
    for parent1_orig, parent2_orig in zip(parents1, parents2):
        if random.random() < cross_rate:
            parent1 = copy.deepcopy(parent1_orig)
            parent2 = copy.deepcopy(parent2_orig)

            tries = 5
            i = 0
            while i < tries:
                # Perform crossover (swap parts of the molecule) [try at least <tries> times]
                modified = parent1.CrossOver(parent2)
                if modified:
                    # Append offsprings
                    offsprings.append(parent1)
                    offsprings.append(parent2)
                    break
                i += 1
    return offsprings


def mutation(molecules, mut_rate):
    offsprings = []
    for molecule in molecules:
        if random.random() < mut_rate:
            molecule_to_mutate = copy.deepcopy(molecule)
            modified = molecule_to_mutate.randomMutate()
            if modified:
                offsprings.append(molecule_to_mutate)
    return offsprings


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


def get_cluster(molecules):
    # Create similarity matrix
    fpgen = AllChem.GetMorganGenerator()
    molecules_fingerprint = [fpgen.GetFingerprint(Chem.MolFromSmiles(atom_to_smiles(molecule.head_atom))) for molecule in molecules]

    # Manually calculate the Tanimoto distance matrix
    n = len(molecules_fingerprint)
    dist_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            tanimoto_sim = DataStructs.TanimotoSimilarity(molecules_fingerprint[i], molecules_fingerprint[j])
            dist_matrix[i, j] = 1 - tanimoto_sim
            dist_matrix[j, i] = dist_matrix[i, j]  # Matrix is symmetric

    condensed_dist_matrix = squareform(dist_matrix)
    # Cluster the similarity matrix while identfying optimal number of clusters
    cluster_linkage = linkage(condensed_dist_matrix, method='average')
    max_clusters = min(len(molecules)//2, 50)
    # Calculate silhouette score for each number of clusters
    silhouette_scores = []
    for n_clusters in range(2, max_clusters + 1):
        labels = fcluster(cluster_linkage, n_clusters, criterion='maxclust')
        score = silhouette_score(dist_matrix, labels, metric='precomputed')
        silhouette_scores.append(score)
    # Optimal number of clusters is the one with the maximum silhouette score
    optimal_clusters = np.argmax(silhouette_scores) + 2
    print(f'Optimal number of clusters: {optimal_clusters}')

    # Cluster the molecules
    labels = fcluster(cluster_linkage, optimal_clusters, criterion='maxclust')
    return labels


def get_scores(population, fitness_scores):
    # (Pareto Ranking, Cluster)
    scores = list(zip(get_pareto_ranking(fitness_scores), get_cluster(population)))
    return scores


def get_admet(molecules):
    molecules = [atom_to_smiles(molecule.head_atom) for molecule in molecules]
    return automated_admet(molecules)


def get_binding_affinity(molecules):
    smiles = [atom_to_smiles(molecule.head_atom) for molecule in molecules]
    model_filepath = "./genetic/trainer.pkl"
    predicted_affinity = get_binding_affinities(smiles, model_filepath)
    
    return predicted_affinity


def get_fitness(molecules):
    # Get Admet and Get SA are merged
    admet_props = get_admet(molecules).values.tolist()

    # Predicted binding affinity (e.g., lower values are better for binding affinity)
    binding_affinity = get_binding_affinity(molecules)
    for i in range(len(admet_props)):
        admet_props[i].append(binding_affinity[i])

    for i in range(len(admet_props)):
        for index in [0, 2, 4, 6, 8, 9, 11, 12]:
            admet_props[i][index] = -admet_props[i][index]
    admet_props = [tuple(admet_prop) for admet_prop in admet_props]

    return admet_props


# Genetic Algorithm for molecular structures
def genetic_algorithm(generations, mut_rate, cross_rate, num_parents):
    # Create log.json if it doesn't exist
    log_path = "./history/log.json"
    if not os.path.exists(log_path):
        with open(log_path, "w") as f:
            json.dump([], f)
    
    # Read log.json
    with open(log_path, "r") as f:
        history = json.load(f)

    if len(history) == 0:
        # Initialize population and archive
        population = init_population()
        pareto_archive = []
    else:
        population = population_from_smiles(history[-1]["population"])
        pareto_archive = history[-1]["pareto_archive"]

    for generation in range(len(history), generations):
        print("Running generation", generation)
        # Get scores
        fitness_scores = get_fitness(population)
        scores = get_scores(population, fitness_scores)

        # Parent Selection [Tournament Selection]
        parents = select_parents(population, scores, num_parents)

        # Mutate
        offsprings_mutate = mutation(population, mut_rate)
        # Crossover
        offsprings_crossover = crossover(parents, cross_rate)

        # New population
        new_population_unfiltered = list(set(population + offsprings_mutate + offsprings_crossover + pareto_archive))

        # Filter all invalid through rdkit
        new_population = []
        for molecule in new_population_unfiltered:
            try:
                Chem.MolFromSmiles(atom_to_smiles(molecule.head_atom))
                new_population.append(molecule)
            except Exception as e:
                print(e)
                pass

        # Get scores
        new_fitness_scores = get_fitness(new_population)
        new_scores = get_scores(new_population, new_fitness_scores)

        pareto_archive = []
        for i in range(len(population)):
            if scores[i][1] == 0:
                pareto_archive.append(population[i])

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
        new_population_ = []
        loops = 0
        while len(new_population_) < POPULATION_SIZE:
            error = 0
            for cluster in sorted_clusters:
                if loops >= len(new_population_clusters[cluster]):
                    error += 1
                    continue

                individual = new_population_clusters[cluster][loops]
                new_population_.append(individual)

                if len(new_population_) >= POPULATION_SIZE:
                    break

            if error >= len(sorted_clusters):
                break
            loops += 1

        # Set new population to population
        population = new_population_.copy()

        # Print best binding affinity
        print(f"Generation {generation}: {[atom_to_smiles(individual.head_atom) for individual in population]}")

        # Save history to log.json
        history.append({
            "population": [atom_to_smiles(individual.head_atom) for individual in population],
            "pareto_archive": [atom_to_smiles(individual.head_atom) for individual in pareto_archive]
        })
        with open(log_path, "w") as f:
            json.dump(history, f)

    print(population)


if __name__ == "__main__":
    genetic_algorithm(GENERATIONS, MUT_RATE, CROSS_RATE, NUM_PARENTS)




