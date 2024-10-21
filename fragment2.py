import json
import os
from collections import Counter

if __name__ == "__main__":
    all_fragments_counter = Counter({})
    for file in os.listdir("data"):
        if file.startswith("fragments-") and file.endswith(".json"):
            with open("data/" + file) as f:
                all_fragments_counter += Counter(json.load(f))

    with open("data/fragmentsAll.json", "w") as f:
        json.dump(dict(all_fragments_counter), f)
