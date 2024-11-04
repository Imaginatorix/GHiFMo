import random
import copy
from Chromosome import *
from GeneticMutation import *
from rdkit import DataStructs
from sklearn.metrics.pairwise import pairwise_distances
from rdkit.Chem import AllChem
from sklearn.metrics import silhouette_score
from scipy.cluster.hierarchy import linkage, fcluster
import numpy as np
import pandas as pd

POPULATION_SIZE = 5
NUM_PARENTS = 5
MUT_RATE = 0.25
CROSS_RATE = 1
GENERATIONS = 200
GP_PATH = 'trained_model.pkl'
gp = joblib.load(GP_PATH)
morgan_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


# Create a simple molecular structure as a genome (analogous to random_genome in basic GA)
def random_molecule():
    C_string_ring = Ring_Manager()
    C_string = Atom("C").add_branches([  # Define the molecular structure
        Branch("single", Atom("N").add_branches([
            Branch("single", Atom("O").add_branches([
                Branch("single", Atom("Cl").add_branches([]))
            ]))
        ]))
    ])
    return C_string

def init_population():
    # Initialize a population of random molecular structures
    # return [random_molecule() for _ in range(pop_size)]
    pass


def select_parents(population, scores, num_parents):
    # Select a parent based on scores [Tournament Selection]
    clustered = {}
    for i in range(len(population)):
        if not scores[i][1] in clustered:
            clustered[scores[i][1]] = []
        clustered[scores[i][1]].append((population[i], scores[i][0]))

    parents = []
    while len(parents) != num_parents:
        for cluster, score in clustered.items():
            parent1 = random.choice(score)
            parent2 = random.choice(score)

            if parent1[1] > parent2[1]:
                parents.append(parent1[0])
            else:
                parents.append(parent2[0])

    return parents


def crossover(parents, cross_rate):
    if parents % 2 != 0:
        raise Exception("Invalid number of parents")

    offsprings = []

    random.shuffle(parents)
    parents1 = parents[:len(parents) // 2]
    parents2 = parents[len(parents) // 2:]
    for parent1_orig, parent2_orig in zip(parents1, parents2):
        if random.random() < cross_rate:
            parent1 = copy.deepcopy(parent1_orig)
            parent2 = copy.deepcopy(parent2_orig)

            i = 0
            while i < 5:
                # Perform crossover (swap parts of the molecule) [try at least n times]
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


def smiles_to_fingerprint(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        fp = morgan_gen.GetFingerprint(mol)
        return np.array(fp)
    return None


def get_binding_affinity(molecule):
    smiles = atom_to_smiles(molecule)
    fingerprint = smiles_to_fingerprint(smiles)

    if fingerprint is None:
        print(f"Invalid molecule: {smiles}")
        return float('-inf')  # Assign a very low score value if invalid

    predicted_affinity = gp.predict([fingerprint])[0]
    return predicted_affinity

    # Get Admet and Get SA are merged. "predicted_affinity" values are directly plugged into the "binding_affinity" in 'def fitness()' for fitness calculation


def get_admet(file_path='C:\\A. Personal Files\\ReSearch\\A. Admet\\selenium\\1234.csv'):
    # Load the CSV file
    try:
        df = pd.read_csv(file_path)
    except FileNotFoundError:
        return "File not found. Please check the file path and name."
    
    # Define the columns to extract
    columns_to_extract = ['smiles', 'Lipinski', 'PPB', 'logVDss', 'CYP3A4-inh', 'CYP3A4-sub', 
                          'CYP2D6-inh', 'CYP2D6-sub', 'cl-plasma', 't0.5', 'DILI', 'hERG', 'Synth']
    
    # Check if all columns exist in the dataset
    missing_columns = [col for col in columns_to_extract if col not in df.columns]
    if missing_columns:
        return f"The following columns are missing in the dataset: {missing_columns}"
    
    # Extract the relevant columns
    extracted_data = df[columns_to_extract]
    
    # Save extracted data to a new CSV file (Optional)
    output_file_path = 'extracted_data.csv'
    extracted_data.to_csv(output_file_path, index=False)
    
    return extracted_data

    # Call the function to test (Optional)
    # extracted_data = get_admet()
    # print(extracted_data)


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
    molecules_fingerprint = [fpgen.GetFingerprint(atom_to_smiles(molecule)) for molecule in molecules]
    dist_matrix = pairwise_distances(molecules_fingerprint, metric=lambda x, y: 1-DataStructs.TanimotoSimilarity(x, y))

    # Cluster the similarity matrix while identfying optimal number of clusters
    cluster_linkage = linkage(dist_matrix, method='average', metric='precomputed')
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


def get_fitness(molecule):
    # Extract ADMET-related properties from the molecule and store them in a tuple
    extracted_data = (
        molecule.get('Lipinski', 0.5),    # Lipinski's rule of five
        molecule.get('PPB', 0.5),         # Plasma Protein Binding
        molecule.get('logVDss', 0.5),     # Volume of Distribution
        molecule.get('CYP3A4-inh', 0.5),  # CYP3A4 inhibition
        molecule.get('CYP3A4-sub', 0.5),  # CYP3A4 substrate
        molecule.get('CYP2D6-inh', 0.5),  # CYP2D6 inhibition
        molecule.get('CYP2D6-sub', 0.5),  # CYP2D6 substrate
        molecule.get('cl-plasma', 0.5),   # Plasma clearance
        molecule.get('t0.5', 0.5),        # Half-life
        molecule.get('DILI', 0.5),        # Drug-Induced Liver Injury
        molecule.get('hERG', 0.5),        # hERG inhibition (cardiotoxicity risk)
        molecule.get('Synth', 0.5)        # Synthetic accessibility
    )
    
    # Predicted binding affinity (e.g., lower values are better for binding affinity)
    binding_affinity = get_binding_affinity(molecule)

    # Set weights for each property
    weights = (
        0.1,  # Lipinski
        0.1,  # PPB
        0.1,  # logVDss
        0.1,  # CYP3A4-inh
        0.1,  # CYP3A4-sub
        0.1,  # CYP2D6-inh
        0.1,  # CYP2D6-sub
        0.1,  # cl-plasma
        0.1,  # t0.5
        0.05, # DILI
        0.05, # hERG
        0.05  # Synth
    )

    # Compute the fitness score using the extracted_data tuple and weights
    fitness_score = sum(value * weight for value, weight in zip(extracted_data, weights))
    fitness_score += (1 / binding_affinity) * 0.2  # Adjusted for binding affinity weight

    return fitness_score


def get_scores(fitness_scores):
    # (Pareto Ranking, Cluster)
    scores = list(zip(get_pareto_ranking(fitness_scores), get_cluster(fitness_scores)))
    return scores


# Genetic Algorithm for molecular structures
def genetic_algorithm(generations, mut_rate, cross_rate, num_parents):
    # Initialize population and archive
    population = init_population()
    # Get scores
    fitness_scores = get_fitness(population)
    scores = get_scores(fitness_scores)

    for generation in range(generations):
        # Parent Selection [Tournament Selection]
        parents = select_parents(population, scores)

        # Mutate
        offsprings_mutate = mutation(parents, mut_rate)
        # Crossover
        offsprings_crossover = crossover(parents, cross_rate)

        # New population
        new_population = population + offsprings_mutate + offsprings_crossover

        # Get scores
        new_fitness_scores = get_fitness(new_population)
        new_scores = get_scores(new_fitness_scores)

        # Population reduction until POPULATION_SIZE
        new_population_ = []
        new_scores_ = []
        sorted_new_population = sorted(dict(zip(new_population, new_scores)).items(), key=lambda x: x[1])

        # Set new population to population
        population = new_population_.copy()
        scores = new_scores_.copy()

        # Print best binding affinity
        print(f"Generation {generation}: Best Fitness = {sorted_new_population[0][1]}")

    print(population)


if __name__ == "__main__":
    genetic_algorithm(GENERATIONS, MUT_RATE, CROSS_RATE, NUM_PARENTS)




