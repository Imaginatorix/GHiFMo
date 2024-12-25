import numpy as np
from Chromosome import *
from util_clusters import get_cluster
from util_binding_affinity_predictor import *
from util_admet_selenium_extraction import automated_admet

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
