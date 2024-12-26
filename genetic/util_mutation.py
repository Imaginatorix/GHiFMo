import random
import copy

def select_parents(population, scores, num_parents):
    # Select a parent based on scores [Tournament Selection]
    clustered = {}
    for i in range(len(population)):
        if not scores[i] in clustered:
            clustered[scores[i]] = []
        clustered[scores[i]].append((population[i], scores[i]))

    parents = []
    while len(parents) < min(len(population), num_parents):
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

