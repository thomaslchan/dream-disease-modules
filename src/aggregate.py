import sys, os, operator, time, math, datetime
import networkx as nx
import numpy as np

# from funcassociate.client import _fc as fc
from sklearn.cluster import KMeans

DSDs = '../data/networks/DSDs'
NODELISTS = '../data/nodelists'
GENEFILE = '../data/ids/gene_ids.txt'

# Takes a path and loads all npy files in that directory into np matrices
def load_networks(path):
    start = time.time()
    adjs = []
    for filename in sorted(os.listdir(path)):
        adj = np.load(os.path.join(path, filename))
        adjs.append(adj)
    end = time.time()
    print("Loaded networks in " +
        str(datetime.timedelta(seconds=int(end-start))))
    return adjs

# Takes a path and reads all nodelists files in that directory into lists
def read_nodelists(path):
    nodelists = []
    for filename in sorted(os.listdir(path)):
        with open(os.path.join(path, filename), "r") as file:
            nodelist = []
            for line in file.readlines():
                nodelist.append(int(line.strip()))
            nodelists.append(nodelist)
    return nodelists

# Takes a filename that corresponds to a list of key value pairs and creates a
# dictionary of gene ids 
def gene_id_dict(filename):
    gene_ids = {}
    with open(filename, "r") as file:
        for line in file.readlines():
            d = line.strip().split("\t")
            gene_ids[int(d[1])] = d[0]
    return gene_ids

# Takes a list of np matrices and a list of nodes of the same length and resizes
# the matrices to the same size, where each gene is the same index in all
# matrices
def resize_networks(networks, nodelists):
    start = time.time()
    len_networks = len(networks)
    adjs = []
    for i in range(0, len_networks):
        len_nodelist = len(networks[i])
        adj = np.zeros((21115, 21115))
        for j in range(0, len_nodelist):
            for k in range(0, len_nodelist):
                adj[j][k] = networks[i][j][k]
        adjs.append(adj)
    end = time.time()
    print("Resized networks in " +
        str(datetime.timedelta(seconds=int(end-start))))
    return adjs

def aggregate_dsds(matrices, impute='mean'):
    impute_func = impute + '_impute'
    imputed_matrices = eval(impute_func)(matrices)
    return np.mean(imputed_matrices, axis=0)

def mean_impute(matrices):
    vals = [np.mean(matrix) for matrix in matrices]
    return impute(matrices, vals)

def min_impute(matrices):
    vals = [np.min(matrix[np.nonzero(matrix)]) for matrix in matrices]
    return impute(matrices, vals)

def max_impute(matrices):
    vals = [np.max(matrix) for matrix in matrices]
    return impute(matrices, vals)

def median_impute(matrices):
    vals = [np.median(matrix) for matrix in matrices]
    return impute(matrices, vals)

def preprocess(matrices):
    rows, num_matrices = len(matrices[0]), len(matrices)
    flat = np.array(matrices).reshape(rows**2, num_matrices)
    flat[flat == 0] = np.nan
    return flat

def mean_local_impute(matrices, global_impute=mean_impute):
    flat = preprocess(matrices)
    mean = np.nanmean(flat, axis=1).reshape(rows, rows)
    vals = global_impute([np.nan_to_num(mean)])
    return local_impute(matrices, vals)

def median_local_impute(matrices, global_impute=mean_impute):
    flat = preprocess(matrices)
    mean = np.nanmedian(flat, axis=1).reshape(rows, rows)
    vals = global_impute([np.nan_to_num(mean)])
    return local_impute(matrices, vals)

def min_local_impute(matrices, global_impute=mean_impute):
    flat = preprocess(matrices)
    mean = np.nanmin(flat, axis=1).reshape(rows, rows)
    vals = global_impute([np.nan_to_num(mean)])
    return local_impute(matrices, vals)

def max_local_impute(matrices, global_impute=mean_impute):
    flat = preprocess(matrices)
    mean = np.nanmax(flat, axis=1).reshape(rows, rows)
    vals = global_impute([np.nan_to_num(mean)])
    return local_impute(matrices, vals)

def local_impute(matrices, vals):
    return [np.copyto(matrix, vals, where=(matrix==0)) for matrix in matrices]


def impute(matrices, imputed_values):
    for matrix, val in zip(matrices, imputed_values):
        matrix[matrix == 0] = val
        np.fill_diagonal(matrix, 0)
    return matrices

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
def cluster(k, agg, gene_ids):
    start = time.time()
    clustering = KMeans(n_clusters=k)
    labels = clustering.fit(agg).labels_
    clusters = d = [[] for x in xrange(k)]
    for i in range(0, labels.size):
           clusters[labels[i]].append(gene_ids[i])
    end = time.time()
    print("Clustered aggregated network in " +
        str(datetime.timedelta(seconds=int(end-start))))
    return clusters

# Takes a p value and a set of clusters and scores each cluster using Func
# Associate. Return the percentage of clusters that are found to be enriched, as
# well as the best GO term of those enriched clusters
def score(pval_cutoff, clusters):
    start = time.time()
    score = 0
    cluster_terms = []
    func = fc()
    for i in range(0, len(clusters)):
        try:
            ret = func.functionate(query=clusters[i], species="Homo sapiens",
            namespace="hgnc_symbol", cutoff=(1-pval_cutoff))["over"]
            p = ret[0][4]*(ret[0][1])
            cluster_terms.append(ret[0][7])
        except:
            p = 1
            cluster_terms.append("")
    if (p < pval_cutoff):
        score += 1
    end = time.time()
    print("Scored clusters in " +
        str(datetime.timedelta(seconds=int(end-start))))
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

def to_file(clusters, cluster_terms, impute_method):
    with open("output.txt", "w+") as file:
        file.write("Imputation Technique: " + impute_method)
        for i in range(0, len(clusters)):
            file.write("Cluster {}:\n".format(i+1))
            for j in range(0, len(clusters[i])):
                file.write("{}\t".format(clusters[i][j]))
                file.write("\n")


    
if __name__ == '__main__':
    # Load data
    networks = load_networks(DSDs)
    nodelists = read_nodelists(NODELISTS)
    gene_ids = gene_id_dict(GENEFILE)

    # Parameters
    wvec = [.2, .2, .1, .1, .3, .1, 1.5]
    k = 200
    pval_cutoff = .05

    # Resize, aggregate, cluster, score, and output
    adjs = resize_networks(networks, nodelists)
    # agg = aggregate(adjs, wvec)
    impute_method = 'mean_local'
    agg = aggregate_dsds(adjs, impute=impute_method)
    clusters = cluster(k, agg, gene_ids)
    #scoreVal, cluster_terms = score(pval_cutoff, clusters)
    to_file(clusters, impute_method)
