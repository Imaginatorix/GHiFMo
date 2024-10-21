import json
import csv

if __name__ == "__main__":
    with open("data/fragmentsAll.json", "r") as f:
        all_fragments = json.load(f)

    # Place frequency in csv to graph
    with open('data/frequency.csv', 'w', newline='') as f:
        for i in all_fragments.values():
            writer = csv.writer(f)
            writer.writerow([i])


