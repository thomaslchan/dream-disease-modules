from __future__ import division
import sys, os, operator, time

import networkx as nx
import numpy as np
from numpy.linalg import inv
from scipy.spatial.distance import pdist, squareform
from sklearn.cluster import SpectralClustering
import gseapy as gp
from prettytable import PrettyTable

DATA = '../data/networks/DSDs'
GENEIDS = '../data/gene_ids.txt'

def load_data(path):
    adjs = []
    for filename in sorted(os.listdir(path)):
        start = time.time()
        adj = np.load(os.path.join(path, filename))
        adjs.append(adj)
        end = time.time()
        print("Loaded network " + filename + " in " 
                + str(float(end-start)) + " seconds.\n")
    return adjs

def gene_id_dict(filename):
    gene_ids = {}
    with open(filename, "r") as file:
	for line in file.readlines():
	    d = line.strip().split("\t")
	    gene_ids[int(d[1])] = d[0]
    return gene_ids

def score(adj, gene_ids):
    labels = SpectralClusterting().fit_predict(adj)
    clusters = {}
    for i in range(0, labels.size):
	try:
	    clusters[labels[i]] = clusters[labels[i]].append(gene_ids[i])
	except:
	    clusters[labels[i]] = []
    scores = []
    for k,v in clusters.items():
	enr = gp.enrichr(gene_list=v)
	scores.append(enr)

def show_eda_table(graphs):
    cols = ['Nodes', 'Edges', 'Size of Largest Component', 
            'Mean Edge Weight', 'Std. Dev. Edge Weight']
    t = PrettyTable(cols)
    for G in graphs:
        info = []
        info += [nx.number_of_nodes(G), nx.number_of_edges(G)]
        try:
            info.append(len(max(nx.connected_components(G), key=len)))
        except:
            info.append(len(max(nx.strongly_connected_components(G), key=len)))
        edge_weights = [d['weight'] for (n1, n2, d) in G.edges(data=True)]
        info += [np.mean(edge_weights), np.std(edge_weights)]
        t.add_row(info)
    print(t)


def sparsify(G, threshold):
    for u,v,data in G.edges(data=True):
        if data['weight'] < threshold:
            G.remove_edge(u,v)
    return G

if __name__ == '__main__':
    adjs = load_data(DATA)
    gene_ids = gene_id_dict(GENEIDS)
    wvec = [0, 0, 0, 0, 0, 0]
    w0 = 0
