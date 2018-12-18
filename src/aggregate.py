import sys, os, operator, time, math, datetime, pickle, argparse
import networkx as nx
import numpy as np
# from funcassociate.client import _fc as fc
from sklearn.cluster import KMeans, SpectralClustering

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
    [km_clusters[label].append(gene_id) for label, gene_id
                                     in zip(km_labels, gene_ids)]
    end = time.time()
    print("Clustered aggregated network in " + str(end - start) + " seconds")

    # print("Clustered aggregated network in " + str(end - start) + " seconds")
    return km_clusters

def score(clusters, pval_cutoff=0.05):
    """
    Calculates score for clustering of genes by statistical significance of
    given p value cutoff using FuncAssociate and returns percentage of enriched
    clusters and GO terms for each cluster

    Args:
    -----------------------------------------------------------------
    - clusters: List of clusters
    - pval_cutoff: Integer for p value cutoff
    """
    start = time.time()
    score = 0
    cluster_terms = []
    for i in range(0, len(clusters)):
        try:
            ret = fc.functionate(query=clusters[i], species="Homo sapiens",
            namespace="hgnc_symbol", cutoff=(1-pval_cutoff))["over"]

            # Result pval is the product of the the log (base 10) of the odds
            # ratio and the number of genes in the query that have the attribute
            pval = ret[0][4] * (ret[0][1])

            # Cluster terms is appending the attribute ID
            cluster_terms.append(ret[0][7])
        except:
            pval = 1
            cluster_terms.append("")
    if (pval < pval_cutoff):
        score += 1
    end = time.time()
    print("Scored clusters in " + str(end - start) + " seconds")
    return score/float(len(clusters)), cluster_terms

def to_file(clusters, impute_method):
    """
    Writes impute method and clusters with corresponding GO terms to a file

    Args:
    -----------------------------------------------------------------
    - clusters: List of clusters
    - impute_method: String of impute method used to generate clusters
    """
    with open("output.txt", "w+") as file:
        file.write("Imputation Technique: " + impute_method + "\n")
        for i in range(0, len(clusters)):
            file.write("Cluster {}:\n".format(i+1))
            for j in range(0, len(clusters[i])):
                file.write("{}\t".format(clusters[i][j]))
                file.write("\n")

# def test_small():
#     size_list = [400, 480, 430]
#     num_nodes = 500
#     networks = [np.random.randint(2, 16, (s, s)).astype(float) for s in size_list]
#     [np.fill_diagonal(n, 0) for n in networks]

#     nodelists = [list(np.random.choice(range(num_nodes), s, replace=False)) for s in size_list]

    
if __name__ == '__main__':
    # Load data
    parser = argparse.ArgumentParser(description='Multigraph Clustering \
                                                of Genes and Proteins')
    parser.add_argument('--impute', type=str, default='mean_local',
                        help="Imputation method")
    parser.add_argument('--n_clusters', type=int, default=250,
                        help="Number of clusters")
    args = parser.parse_args()
   
    # nodelists = get_nodelists(read_graphs())
    # with open ('nodes.pkl', 'wb') as f:
    #     pickle.dump(nodelists, f)
    # print(nodelists[0])

    networks = load_DSDs(DSDs)
    nodelists = read_nodelists(NODELISTS)
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

    score_val, cluster_terms = score(clusters)
    print("Score: " + str(score_val) + "\n")
    to_file(clusters, impute_method)

    #CONVERT
    # max = pickle.load(open('../results/max.pkl', 'rb'))
    # pickle.dump(max, open('MAX.pkl', 'wb'), protocol=2)

    # max_local = pickle.load(open('../results/max_local.pkl', 'rb'))
    # pickle.dump(max_local, open('MAX_LOCAL.pkl', 'wb'), protocol=2)

    # mean = pickle.load(open('../results/mean.pkl', 'rb'))
    # pickle.dump(mean, open('MEAN.pkl', 'wb'), protocol=2)

    # mean_local = pickle.load(open('../results/mean_local.pkl', 'rb'))
    # pickle.dump(mean_local, open('MEAN_LOCAL.pkl', 'wb'), protocol=2)

    # min = pickle.load(open('../results/min.pkl', 'rb'))
    # pickle.dump(min, open('MIN.pkl', 'wb'), protocol=2)

    # min_local = pickle.load(open('../results/min_local.pkl', 'rb'))
    # pickle.dump(min_local, open('MIN_LOCAL.pkl', 'wb'), protocol=2)

    # zero = pickle.load(open('../results/zero.pkl', 'rb'))
    # pickle.dump(zero, open('ZERO.pkl', 'wb'), protocol=2)

    #LOAD
    # max = pickle.load(open('MAX.pkl', 'rb'))
    # max_score, cluster_terms = score(max)
    # print('MAX: ' + str(max_score))

    # max_local = pickle.load(open('MAX_LOCAL.pkl', 'rb'))
    # max_local_score, cluster_terms = score(max_local)
    # print('MAX_LOCAL: ' + str(max_local_score))

    # mean = pickle.load(open('MEAN.pkl', 'rb'))
    # mean_score, cluster_terms = score(mean)
    # print('MEAN: ' + str(mean_score))

    # mean_local = pickle.load(open('MEAN_LOCAL.pkl', 'rb'))
    # mean_local_score, cluster_terms = score(mean_local)
    # print('MEAN_LOCAL: ' + str(mean_local_score))

    # min = pickle.load(open('MIN.pkl', 'rb'))
    # min_score, cluster_terms = score(min)
    # print('MIN: ' + str(min_score))

    # min_local = pickle.load(open('MIN_LOCAL.pkl', 'rb'))
    # min_score, cluster_terms = score(min_local)
    # print('MIN_LOCAL: ' + str(min_score))

    # zero = pickle.load(open('ZERO.pkl', 'rb'))
    # zero_score, cluster_terms = score(zero)
    # print('MIN_LOCAL: ' + str(zero_score))
