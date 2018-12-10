from __future__ import division
import sys, os, operator, time, math, datetime
import networkx as nx
import numpy as np
import mpmath as mp

from funcassociate.client import _fc as fc
from numpy.linalg import inv
from scipy.spatial.distance import pdist, squareform
from sklearn.cluster import KMeans

DSDs = '../data/networks/DSDs'
NODELISTS = '../data/nodelists'
GENEFILE = '../data/ids/gene_ids.txt'
GOFILE = '../data/go/go_terms.txt'

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

def read_nodelists(path):
    nodelists = []
    for filename in sorted(os.listdir(path)):
        with open(os.path.join(path, filename), "r") as file:
            nodelist = []
            for line in file.readlines():
                nodelist.append(int(line.strip()))
            nodelists.append(nodelist)
    return nodelists

def gene_id_dict(filename):
    gene_ids = {}
    with open(filename, "r") as file:
        for line in file.readlines():
            d = line.strip().split("\t")
            gene_ids[int(d[1])] = d[0]
    return gene_ids

def resize_networks(networks, nodelists):
    start = time.time()
    len_networks = len(networks)
    adjs = []
    for i in range(0, len_networks):
        len_nodelist = len(nodelists[i])
        adj = np.zeros((21115, 21115))
        for j in range(0, len_nodelist):
            for k in range(0, len_nodelist):
                adj[j][k] = networks[i][j][k]
        adjs.append(adj)
    end = time.time()
    print("Resized networks in " +
        str(datetime.timedelta(seconds=int(end-start))))
    return adjs


def aggregate(adjs, wvec):
    start = time.time()
    agg = np.full_like(adjs[0], (wvec[-1]))
    for i in range(0, len(adjs)):
        agg = np.add(agg, np.multiply(adjs[i], (wvec[i])))
    end = time.time()
    print("Aggregated networks in " +
        str(datetime.timedelta(seconds=int(end-start))))
    return agg

def cluster(k, adj, gene_ids):
    start = time.time()
    clustering = KMeans(n_clusters=k)
    labels = clustering.fit(adj).labels_
    clusters = d = [[] for x in xrange(k)]
    for i in range(0, labels.size):
           clusters[labels[i]].append(gene_ids[i])
    end = time.time()
    print("Clustered aggregated network in " +
        str(datetime.timedelta(seconds=int(end-start))))
    return clusters

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
            print(ret[0])
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

def to_file(clusters, cluster_terms, scoreVal, wvec):
    start = time.time()
    with open("output.txt", "w+") as file:
	file.write("{}%\tWeights: {}\tRegularization: {}\n"
            .format(scoreVal * 100, wvec[:-1], wvec[-1]))
	for i in range(0, len(clusters)):
	    file.write("Cluster {}: {}\n".format(i+1, cluster_terms[i]))
	    for j in range(0, len(clusters[i])):
		file.write("{}\t".format(clusters[i][j]))
	    file.write("\n")
    end = time.time()
    print("Wrote clusters in " +
        str(datetime.timedelta(seconds=int(end-start))))
	

if __name__ == '__main__':
    networks = load_networks(DSDs)
    nodelists = read_nodelists(NODELISTS)
    gene_ids = gene_id_dict(GENEFILE)

    wvec = [.2, .2, .1, .1, .3, .1, 1.5]
    k = 300
    pval_cutoff = .05
    adjs = resize_networks(networks, nodelists)
    agg = aggregate(adjs, wvec)
    clusters = cluster(k, agg, gene_ids)
    scoreVal, cluster_terms = score(pval_cutoff, clusters)
    to_file(clusters, cluster_terms, scoreVal, wvec)

"""
def go_term_dict(filename):
    go_terms= {}
    with open(filename, "r") as file:
        for line in file.readlines():
            d = line.strip().split("\t")
            t = d[1].split(" ")
            g = []
            for i in range(0, len(t)):
                g.append(map(int, t[i].split(",")))
            go_terms[d[0]] = g
    return go_terms
"""
