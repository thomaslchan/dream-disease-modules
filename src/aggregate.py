import sys, os, operator, time, math, datetime, pickle, argparse, json
import networkx as nx
import numpy as np
from sklearn.cluster import KMeans

from impute_methods import *
from load_data import *

NUM_GENES = 21115


def resize_networks(networks, nodelists, total_nodes=NUM_GENES):
    """
    Resizes networks' distance matrices to fit in a total_nodes by 
    total_nodes matrix. Returns a list of matrices as Numpy arrays.

    Args:
    -----------------------------------------------------------------
    - networks: A list of distance matrices for the networks
    - nodelists: A list of node-lists (which are integer lists)
    - total_nodes: An integer for the total nodes (default=NUM_GENES)
    """
    start = time.time()
    matrices = []
    for network, nodelist in zip(networks, nodelists):
        num_nodes = len(nodelist)
        matrix = np.zeros((total_nodes, total_nodes))
        for j in range(0, num_nodes):
            for k in range(0, num_nodes):
                matrix[nodelist[j]][nodelist[k]] = network[j][k]
        matrices.append(matrix)
    end = time.time()
    print("Resized networks in " + str(end - start) + " seconds")
    return matrices


def aggregate_dsds(matrices, impute='mean'):
    """
    Aggregates distance matrices using some imputation method and returns
    a single matrix as a Numpy array.

    Args:
    -----------------------------------------------------------------
    - matrices: List of matrices as Numpy arrays
    - impute: Keyword for an impute function in impute_methods
    """
    impute_func = impute + '_impute'
    imputed_matrices = eval(impute_func)(matrices)
    return np.mean(imputed_matrices, axis=0)

def cluster(k, agg_matrix, gene_ids):
    """
    Clusters an aggregated distance matrix with K-Means clustering and returns
    the clusters as lists of HGNC gene symbols.

    Args:
    -----------------------------------------------------------------
    - k: Number of clusters, integer
    - agg_matrix: Distance matrix, Numpy array
    - gene_ids: List of gene IDs
    """
    k_means = KMeans(n_clusters=k)


    start = time.time()
    km_labels = k_means.fit(agg_matrix).labels_
    km_clusters = [[] for x in range(k)]
    [km_clusters[label].append(gene_ids[gene_id]) for label, gene_id
                                     in zip(km_labels, gene_ids)]
    end = time.time()
    print("Clustered aggregated network in " + str(end - start) + " seconds")

    return km_clusters

def convert_to_py2_pickle(path):
    """
    Converts pickle of python3 protocol to a pickle of python2 protocol.

    Args:
    -----------------------------------------------------------------
    - path: The path of the folder containing python3 pickle files
    """
    start = time.time()
    for filename in sorted(os.listdir(path)):
        print(filename)
        start = time.time()
        p3_pickle = pickle.load(open(os.path.join(path, filename), 'rb'))
        pickle.dump(p3_pickle, open(filename, 'wb'), protocol=2)
    end = time.time()
      
    print("Converted pickles in " + str(int(end-start)) + " seconds.\n\n")

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Multigraph Clustering \
                                                of Genes and Proteins')
    parser.add_argument('--impute', type=str, default='mean_local',
                        help="Imputation method")
    parser.add_argument('--n_clusters', type=int, default=250,
                        help="Number of clusters")
    args = parser.parse_args()
   
    nodelists = load_nodelists()
    networks = load_DSDs()
    gene_ids = gene_id_dict(GENEFILE)

    # Resize, aggregate, cluster, score, and output
    matrices = resize_networks(networks, nodelists)
    print("Resized networks")

    impute_method = args.impute
    agg = aggregate_dsds(matrices, impute=impute_method)
    print("Aggregated matrices")
    
    clusters = cluster(args.n_clusters, agg, gene_ids)
    print("Finished clustering")
    with open(impute_method + '.pkl', 'wb') as f:
        pickle.dump(clusters, f)

