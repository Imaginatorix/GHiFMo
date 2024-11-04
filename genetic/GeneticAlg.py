import random
import copy
import numpy as np
from Chromosome import *
from GeneticMutation import *
import pandas as pd
from rdkit import Chem
from rdkit.chem import rdFingerprintGenerator
import joblib

pop_size = 5
gen_length = 20  # placeholder for molecular complexity
mut_rate = 0.01  # probability of mutation happening (1% in this case)
cross_rate = 0.07  # probability of doing crossover instead of returning parents
generations = 200  # number of processes

model_filepath = 'trainer.pkl'
gp = joblib.load(model_filepath)

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

def molecule_to_smiles(molecule):
    return molecule.to_smiles()  # Assuming each molecule has a `to_smiles()` method
    
def smiles_to_fingerprint(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        fp = morgan_gen.GetFingerprint(mol)
        return np.array(fp)
    return None

# Initialize a population of random molecular structures
def init_population():
    # return [random_molecule() for _ in range(pop_size)]
    pass

# Updated fitness function: calculate the length of the molecule using the provided get_length function
def fitness(molecule):
    # Assuming each molecule has the method `get_length()` to return the molecular string length
    return Mutate(molecule, Ring_Manager()).get_length()

# Select a parent based on fitness values
def select_parent(population, fitness_values):
    total_fitness = sum(fitness_values)
    pick = random.uniform(0, total_fitness)
    current = 0
    
    for individual, fitness_value in zip(population, fitness_values):
        current += fitness_value
        if current > pick:
            return individual

# Crossover two molecules by swapping branches or other structural elements
def crossover(parent1, parent2):
    if random.random() < cross_rate:
        parent1_mutate = Mutate(parent1, Ring_Manager())
        parent2_mutate = Mutate(parent2, Ring_Manager())
        
        # Perform crossover (swap parts of the molecule)
        parent1_mutate.CrossOver(parent2_mutate)
        
        # Return the resulting offspring
        return parent1_mutate.head_atom, parent2_mutate.head_atom
    else:
        return parent1, parent2

# Mutate a molecule using chemical mutation operations
def mutate(molecule):
    molecule_mutate = Mutate(molecule, Ring_Manager())
    
    # if random.random() < mut_rate:
    mutation_type = random.choice(['bond', 'atom', 'ring', 'branch'])
    
    i = 0
    while i <= 5:
        try:
            if mutation_type == 'bond':
                molecule_mutate.MutateBond()
            elif mutation_type == 'atom':
                molecule_mutate.MutateAtom()
            elif mutation_type == 'ring':
                if random.choice([True, False]):
                    molecule_mutate.CloseRing()
                else:
                    molecule_mutate.OpenRing()
            elif mutation_type == 'branch':
                molecule_mutate.MutateBranch()
            break
        except Exception:
            continue

    return molecule_mutate.head_atom

def get_binding_affinity():
    smiles = molecule_to_smiles(molecule)
    fingerprint = smiles_to_fingerprint(smiles)

    if fingerprint is None:
        print(f"Invalid molecule: {smiles}")
        return float('-inf')  # Assign a very low score value if invalid

    predicted_affinity = gp.predict([fingerprint])[0]
    return predicted_affinity

# Get Admet and Get SA are merged. "predicted_affinity" values are directly plugged into the "binding_affinity" in 'def fitness()' for fitness calculation

# Change file path as to where the stored ADMET and SA values are
def get_admet(file_path='C:\A. Personal Files\ReSearch\A. Admet\selenium\merged_excel.xlsx'):
    try:
        df = pd.read_excel(file_path)
    except FileNotFoundError:
        return "File not found. Please check the file path and name."
    
    # Define the columns to extract#No CaCo2 and HIA
    columns_to_extract = ['smiles' , 'Lipinski','PPB', 'logVDss', 'CYP3A4-inh', 'CYP3A4-sub', 
                          'CYP2D6-inh', 'CYP2D6-sub', 'cl-plasma', 't0.5', 'DILI', 'hERG', 'Synth']
     
    # Check if all columns exist in the dataset
    missing_columns = [col for col in columns_to_extract if col not in df.columns]
    if missing_columns:
        return f"The following columns are missing in the dataset: {missing_columns}"
    
        # Extract the relevant columns
        extracted_data = df[columns_to_extract]
        # No output yet. Not stored in .csv nor .xlsx file. 
   
        return extracted_data

    # get_sa() # Uncomment this line to run it in a local environment
    
    extracted_data = get_admet()

    #Optional
    print(extracted_data)

    # When conversion is needed, just uncomment:
    # extracted_columns = get_admet()
    # extracted_columns.to_excel('extracted_data.xlsx', index=False)

def fitness(molecule):
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

    
# Genetic Algorithm for molecular structures
def genetic_algorithm():

    # Initialize population
    population = init_population()
    # Get Fitness

    for generation in range(generations):
        # Parent Selection [Tournament Selection]

        # Mutate

        # Crossover

        # Scaffold Hop

        # New population

        # Pareto Ranking, Filtration, and Clustering [Pareto Archive]

        # Set new population to population





        # for molecule in population:
        #     print(molecule.show())

        new_population = population
        for i in range(100):
            molecule_to_mutate = random.choice(population)
            new_population.append(mutate(molecule_to_mutate))

        print(new_population, "AHHH")

        # for _ in range(pop_size // 2):
        #     parent1 = select_parent(population, fitness_values)
        #     parent2 = select_parent(population, fitness_values)

        #     offspring1, offspring2 = crossover(copy.deepcopy(parent1), copy.deepcopy(parent2))
        #     new_population.extend([mutate(offspring1), mutate(offspring2)])

        fitness_values = [fitness(molecule) for molecule in population]
        mol_fitness = dict(zip(new_population, fitness_values))
        best_fitness = max(fitness_values)
        print(f"Generation {generation}: Best Fitness = {best_fitness}")

        sorted_pop = sorted(mol_fitness.items(), key=lambda x: x[1])
        print(sorted_pop)

    best_index = fitness_values.index(max(fitness_values))
    best_solution = population[best_index]
    print('Best Solution:')
    best_solution.show()  # Assuming show() prints the molecular structure
    print(f'Best Fitness: {fitness(best_solution)}')


if __name__ == "__main__":
    genetic_algorithm()




