import csv
import json
import os
from rdkit import Chem
from rdkit.Chem import Recap
from collections import Counter

def write2json(location, batch, batch_name):
    with open(location+"-"+batch_name+".json", "w") as f:
        json.dump(batch, f)


if __name__ == "__main__":
    BATCH_SIZE = 1000
    SOURCE_LOCATION = "data/cleaned.csv"
    DATABASE_LOCATION = "data/fragments"

    # If database does not exist yet
    if not os.path.exists(DATABASE_LOCATION+".txt"):
        with open(DATABASE_LOCATION+".txt", "w") as f:
            f.write("0")


    with open(SOURCE_LOCATION, "r") as f:
        smiles = csv.reader(f)
        next(smiles) # Skip header

        # Total entries
        total_counts = 15913 # Hard-coded... but you can get it from len(smiles.readlines()). I did this cause I need it to be batching (cause of how many it is and how weak my laptop is)

        # Skip already done
        with open(DATABASE_LOCATION+".txt", "r") as f:
            last_batch = int(f.readline())
            for _ in range(last_batch):
                next(smiles)

        while last_batch < total_counts:
            fragments = []

            # By batch
            i = 0
            while last_batch < total_counts and i < BATCH_SIZE:
                mol = next(smiles)[0]
                print(last_batch, mol)
                m = Chem.MolFromSmiles(mol)
                res = Recap.RecapDecompose(m)
                leaves = res.GetLeaves()
                fragments += list(leaves.keys())

                i += 1
                last_batch += 1

            # Save fragments to json
            write2json(DATABASE_LOCATION, dict(Counter(fragments)), str(last_batch))
            # Log assay number
            with open(DATABASE_LOCATION+".txt", "w") as f:
                f.write(str(last_batch))

            print(f"{last_batch}/{total_counts} molecules fragmented...")



