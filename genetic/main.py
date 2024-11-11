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
from rdkit.Chem import AllChem, rdFingerprintGenerator
from sklearn.metrics import silhouette_score, pairwise_distances
from Chromosome import *
from GeneticMutation import *
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
import joblib
from automated_admet import automated_admet  
from admet_selenium_extraction import automated_admet_extraction
from binding_affinity_extractor import get_binding_affinities
POPULATION_SIZE = 5
NUM_PARENTS = 5
MUT_RATE = 0.25
CROSS_RATE = 1
GENERATIONS = 200


def binding_affinity(molecules):
    smiles = [atom_to_smiles(molecule) for molecule in molecules]
    model_filepath = 'genetic\trainer.pkl'
    predicted_affinity = get_binding_affinities(smiles, model_filepath)
    
    return predicted_affinity
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


def smiles_to_fingerprint(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        fp = morgan_gen.GetFingerprint(mol)
        return np.array(fp)
    return None


def init_population():
    # Initialize a population of random molecular structures
    # return [atom_from_smiles("S=c1nc[nH]c2nc[nH]c12")]
    population_chromosomes = [atom_from_smiles("Nc1nc(-c2cccs2)c(NC(=O)c2ccccc2)s1"), atom_from_smiles("CC(C)OC(=O)CC(=O)CSc1nc2c(cc1C#N)CCC2"),
            atom_from_smiles("N#CC(C#N)=c1ccc2c(c1)NC1(CCCCC1)N=2"), atom_from_smiles("CCc1n[nH]c2c1/C(=N\O)CC(c1ccccc1)C2")]

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


def smiles_to_fingerprint(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
    except Exception:
        mol = None

    if mol is not None:
        fp = morgan_gen.GetFingerprint(mol)
        return np.array(fp)
    return None

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


def get_admet():
    # Define paths and parameters
    smiles_file_path = r"C:\A. Personal Files\ReSearch\Final\download\smiles.xlsx"
    chromedriver_path = r"C:\A. Personal Files\ReSearch\Final\chromedriver-win64\chromedriver.exe"
    download_folder = r"C:\A. Personal Files\ReSearch\Final\download"
    batch_size = 15
    columns_to_extract = [
        'Lipinski', 'PPB', 'logVDss', 'CYP3A4-inh', 'CYP3A4-sub', 
        'CYP2D6-inh', 'CYP2D6-sub', 'cl-plasma', 't0.5', 'DILI', 'hERG', 'Synth'
    ]
    
    # Load SMILES data
    try:
        df = pd.read_excel(smiles_file_path)
        smiles_list = df['SMILES'].tolist()
    except FileNotFoundError:
        print("SMILES file not found. Please check the file path.")
        return None
    except KeyError:
        print("The specified 'SMILES' column was not found in the file.")
        return None

    # Set up Selenium options
    options = webdriver.ChromeOptions()
    options.add_experimental_option("detach", True)
    prefs = {"download.default_directory": download_folder}
    options.add_experimental_option("prefs", prefs)

    # Initialize the WebDriver
    service = Service(chromedriver_path)
    driver = webdriver.Chrome(service=service, options=options)
    
    try:
        # Open ADMET Lab website
        driver.get("https://admetlab3.scbdd.com/server/screening")
        WebDriverWait(driver, 10).until(EC.presence_of_element_located((By.ID, "profile-tab"))).click()

        # Process SMILES in batches
        for i in range(0, len(smiles_list), batch_size):
            smiles_batch = smiles_list[i:i + batch_size]
            smiles_string = "\n".join(smiles_batch)

            # Enter SMILES data into the input field
            search = WebDriverWait(driver, 10).until(
                EC.presence_of_element_located((By.ID, "exampleFormControlTextarea1"))
            )
            search.clear()
            search.send_keys(smiles_string)
            time.sleep(3)

            # Submit the batch for processing
            submit_button = WebDriverWait(driver, 10).until(
                EC.element_to_be_clickable((By.CSS_SELECTOR, "button.btn.btn-success"))
            )
            submit_button.click()

            # Wait for results to load and download the results
            WebDriverWait(driver, 10).until(EC.presence_of_element_located((By.CLASS_NAME, "text-center")))
            download_button = WebDriverWait(driver, 10).until(
                EC.element_to_be_clickable((By.CLASS_NAME, "btn-outline-success"))
            )
            download_button.click()

            # Wait for download to complete
            time.sleep(5)
            driver.get("https://admetlab3.scbdd.com/server/screening")
            WebDriverWait(driver, 10).until(EC.presence_of_element_located((By.ID, "profile-tab"))).click()

    finally:
        # Close the driver after all downloads are complete
        driver.quit()

    # Merge downloaded CSV files
    csv_files = glob.glob(os.path.join(download_folder, '*.csv'))
    if not csv_files:
        print("No CSV files were downloaded.")
        return None

    merged_df = pd.concat((pd.read_csv(f) for f in csv_files), ignore_index=True)
    merged_output_path = os.path.join(download_folder, 'merged_output.csv')
    merged_df.to_csv(merged_output_path, index=False)
    print("Files have been merged into:", merged_output_path)

    # Clean up the downloaded CSV files to avoid redundancy
    for f in csv_files:
        os.remove(f)

    # Check if all columns exist in the merged dataset
    missing_columns = [col for col in columns_to_extract if col not in merged_df.columns]
    if missing_columns:
        print(f"The following columns are missing in the dataset: {missing_columns}")
        return None

    # Extract the relevant columns
    extracted_data = merged_df[columns_to_extract].values
    print("Extracted Data:", extracted_data)
    return extracted_data

# Run the ADMET automation and data extraction
get_admet()


def get_fitness(molecule):
    # Get Admet and Get SA are merged. "predicted_affinity" values are directly plugged into the "binding_affinity" in 'def fitness()' for fitness calculation
    admet_props = get_admet().tolist()
    admet_props = [tuple(admet_prop) for admet_prop in admet_props]

    # Predicted binding affinity (e.g., lower values are better for binding affinity)
    # binding_affinity = get_binding_affinity(molecule)

    return admet_props


# Genetic Algorithm for molecular structures
def genetic_algorithm(generations, mut_rate, cross_rate, num_parents):
    # Initialize population and archive
    population = init_population()
    # Get scores
    fitness_scores = get_fitness(population)
    scores = get_scores(population, fitness_scores)

    for generation in range(generations):
        # Parent Selection [Tournament Selection]
        parents = select_parents(population, scores, num_parents)

        # Mutate
        offsprings_mutate = mutation(parents, mut_rate)
        # Crossover
        offsprings_crossover = crossover(parents, cross_rate)

        # New population
        new_population_unfiltered = population + offsprings_mutate + offsprings_crossover

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

        # Population reduction until POPULATION_SIZE
        new_population_ = []
        new_scores_ = []
        sorted_new_population = sorted(dict(zip(new_population, new_scores)).items(), key=lambda x: x[1])
        for individual, score in sorted_new_population[:POPULATION_SIZE]:
            new_population_.append(individual)
            new_scores_.append(score)

        # Set new population to population
        population = new_population_.copy()
        scores = new_scores_.copy()

        # Print best binding affinity
        print(f"Generation {generation}: {population}")

    print(population)


if __name__ == "__main__":
    genetic_algorithm(GENERATIONS, MUT_RATE, CROSS_RATE, NUM_PARENTS)




