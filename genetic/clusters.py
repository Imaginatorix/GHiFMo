def select_balanced_population(mol_data, population_size):
    """
    Selects molecules from each cluster, prioritizing lower ranks.
    
    Parameters:
    mol_data (list of tuples): Each entry is a tuple in the form (Mol, (rank, cluster)).
    population_size (int): The size of the new population to select.

    Returns:
    list: Selected molecules with equal representation from each cluster, sorted by rank.
    """
    
    # Sort by rank within each cluster
    from collections import defaultdict
    
    cluster_dict = defaultdict(list)
    
    # Group molecules by cluster and sort by rank
    for mol, (rank, cluster) in mol_data:
        cluster_dict[cluster].append((rank, mol))
    
    # Sort each cluster list by rank
    for cluster in cluster_dict:
        cluster_dict[cluster].sort(key=lambda x: x[0])  # Sort by rank ascending
    
    # Select molecules equally from each cluster
    selected_population = []
    cluster_keys = list(cluster_dict.keys())
    
    i = 0
    while len(selected_population) < population_size:
        cluster_key = cluster_keys[i % len(cluster_keys)]
        if cluster_dict[cluster_key]:
            selected_population.append(cluster_dict[cluster_key].pop(0)[1])
        i += 1

    return selected_population


# Sample usage
mol_data = [
    # (Mol, (rank, cluster)) examples
    ('MolA', (1, 'Cluster1')),
    ('MolB', (2, 'Cluster1')),
    ('MolC', (1, 'Cluster2')),
    ('MolD', (2, 'Cluster2')),
    ('MolE', (1, 'Cluster3')),
    # Add more molecules with rank and cluster information as needed
]

population_size = 5  # Desired size of the selected population
selected_population = select_balanced_population(mol_data, population_size)

print("Selected Population:", selected_population)
