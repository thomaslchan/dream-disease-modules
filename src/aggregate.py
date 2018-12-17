import sys, os, operator, time, math, datetime, pickle

import networkx as nx
import numpy as np
# from funcassociate.client import _fc as fc
from sklearn.cluster import KMeans

from impute_methods import *
from load_data import *

NUM_GENES = 21115



# Takes a list of np matrices and a list of nodes of the same length and resizes
# the matrices to the same size, where each gene is the same index in all
# matrices
def resize_networks(networks, nodelists):
    start = time.time()
    matrices = []
    for network, nodelist in zip(networks, nodelists):
        num_nodes = len(nodelist)
        matrix = np.zeros((NUM_GENES, NUM_GENES))
        for j in range(0, num_nodes):
            for k in range(0, num_nodes):
                matrix[nodelist[j]][nodelist[k]] = network[j][k]
        matrices.append(matrix)
    end = time.time()
    print("Resized networks in " + str(end - start) + " seconds")
    return matrices
    return adjs

def aggregate_dsds(matrices, impute='mean'):
    impute_func = impute + '_impute'
    imputed_matrices = eval(impute_func)(matrices)
    return np.mean(imputed_matrices, axis=0)


# Takes a list of np matrices and a weight vector and sums all matrices together
# scaled by the weights in the vector and added to a regularizer value.
def aggregate(adjs, wvec):
    start = time.time()
    agg = np.full_like(adjs[0], (wvec[-1]))
    for i in range(0, len(adjs)):
        agg = np.add(agg, np.multiply(adjs[i], (wvec[i])))
    end = time.time()
    print("Aggregated networks in " +
        str(datetime.timedelta(seconds=int(end-start))))
    return agg

# Takes a k value, an aggregated matrix, and a dictionary of gene ids. Uses
# k-means clustering on the aggregated matrix to produce clusters of genes
def cluster(k, agg_matrix, gene_ids):
    start = time.time()
    k_means = KMeans(n_clusters=k)
    labels = k_means.fit(agg_matrix).labels_
    clusters = [[] for x in range(k)]
    [clusters[label].append(gene_id) for label, gene_id
                                     in zip(labels, gene_ids)]
    end = time.time()
    print("Clustered aggregated network in " + str(end - start) + " seconds")
    return clusters

# Takes a p value and a set of clusters and scores each cluster using Func
# Associate. Return the percentage of clusters that are found to be enriched, as
# well as the best GO term of those enriched clusters
def score(pval_cutoff, clusters):
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

# Takes a set of clusters, their GO terms, a score, and the weight vector used
# to produce these results and writes this information to file
# def to_file(clusters, cluster_terms, scoreVal, wvec):
#     with open("output.txt", "w+") as file:
#   file.write("{}%\tWeights: {}\tRegularization: {}\n"
#             .format(scoreVal * 100, wvec[:-1], wvec[-1]))
#   for i in range(0, len(clusters)):
#       file.write("Cluster {}: {}\n".format(i+1, cluster_terms[i]))
#       for j in range(0, len(clusters[i])):
#       file.write("{}\t".format(clusters[i][j]))
#       file.write("\n")
#     end = time.time()

def to_file(clusters, impute_method):
    with open("output.txt", "w+") as file:
        file.write("Imputation Technique: " + impute_method)
        for i in range(0, len(clusters)):
            file.write("Cluster {}:\n".format(i+1))
            for j in range(0, len(clusters[i])):
                file.write("{}\t".format(clusters[i][j]))
                file.write("\n")


    
if __name__ == '__main__':
    # Load data
    parser = argparse.ArgumentParser(description='Multigraph Clustering \
                                                of Genes and Proteins')
    parser.add_argument('--impute', type=str, default='mean_local',
                        help="Imputation method")
    parser.add_argument('--n_clusters', type=int, default=250,
                        help="Number of clusters")
    args = parser.parse_args()
    
    networks = load_networks(DSDs)
    nodelists = read_nodelists(NODELISTS)
    gene_ids = gene_id_dict(GENEFILE)

    # Parameters
    # wvec = [.2, .2, .1, .1, .3, .1, 1.5]
    pval_cutoff = .05

    # Resize, aggregate, cluster, score, and output
    adjs = resize_networks(networks, nodelists)
    # agg = aggregate(adjs, wvec)
    impute_method = args.impute
    agg = aggregate_dsds(adjs, impute=impute_method)
    clusters = cluster(args.n_clusters, agg, gene_ids)
    with open(impute_method + '.pkl', 'wb') as f:
        pickle.dump(clusters, f)
    #scoreVal, cluster_terms = score(pval_cutoff, clusters)
    # to_file(clusters, impute_method)
