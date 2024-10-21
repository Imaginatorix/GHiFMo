import json
import sys
sys.path.append("genetic")
from GeneticMutation import *
from Chromosome import *


if __name__ == "__main__":
    with open("data/fragments_transformations_frequency.json", "r") as f:
        all_fragments = json.load(f)

    lengths = {}
    # Get length of fragments and then place them in a list (with length as key and the fragments as values)
    for fragment in all_fragments.keys():
        print(fragment)
        # TODO: Fix valence issues
        try:
            fragment_hier = atom_from_smiles(fragment)
        except Exception:
            continue

        length = Mutate(fragment_hier[0], fragment_hier[1]).get_length()
        if not length in lengths:
            lengths[length] = []
        
        lengths[length].append(fragment)

    # Save
    with open("data/fragments_transformations_length.json", "w") as f:
        json.dump(lengths, f)

