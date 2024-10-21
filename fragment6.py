import json

if __name__ == "__main__":
    with open("data/fragments_transformations_frequency.json", "r") as f:
        all_fragments = json.load(f)

    # Find all fragments with only one attachment point
    all_fragments = dict(filter(lambda x: x[0].count("*") != 2, all_fragments.items()))

    # Save
    with open("data/fragments_transformations_branches.json", "w") as f:
        json.dump(all_fragments, f)


