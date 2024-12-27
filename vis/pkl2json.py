from Chromosome import *
from GeneticMutation import *
import pickle
import json

from_path = "../history/log.pkl"
to_path = "../history/log.json"

with open(from_path, "rb") as f:
    up = pickle.Unpickler(f)
    history = up.load()

# Remove population molecules
for gen in history:
    del gen["population_molecules"]

with open(to_path, "w") as f:
    json.dump(history, f)

