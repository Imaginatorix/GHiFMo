import json

if __name__ == "__main__":
    with open("data/fragmentsAll.json", "r") as f:
        all_fragments = json.load(f)

    required_frequency = 2
    # Remove molecules below the required frequency
    all_fragments = dict(filter(lambda x: x[1] >= 2, all_fragments.items()))
    # Remove unsupported molecules
    all_fragments = dict(filter(lambda x: not "+" in x[0], all_fragments.items()))
    # Remove unsupported molecules
    all_fragments = dict(filter(lambda x: not "-" in x[0], all_fragments.items()))
    # Remove fragments with more than two attachment points
    all_fragments = dict(filter(lambda x: x[0].count("*") <= 2, all_fragments.items()))

    # Save
    with open("data/fragments_transformations_frequency.json", "w") as f:
        json.dump(all_fragments, f)


