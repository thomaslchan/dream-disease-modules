from __future__ import division
import sys, os, operator, time

import networkx as nx
import numpy as np
from numpy.linalg import inv
from scipy.spatial.distance import pdist, squareform
from sklearn.cluster import SpectralClustering
import gseapy as gp
from prettytable import PrettyTable

DATA = '../data/networks/original/'
GENEIDS = '../data/gene_ids.txt'

def load_data(path):
    graphs = []
    for filename in sorted(os.listdir(path)):
        print(filename)
        start = time.time()
        G = nx.read_weighted_edgelist(os.path.join(path, filename))
        end = time.time()
        graphs.append(G)
        print("Loaded network " + filename + " in " 
                + str(float(end-start)) + " seconds.\n")
    return graphs

def gene_id_dict(filename):
    gene_ids = {}
    with open(filename, "r") as file:
        for line in file.readlines():
            d = line.strip().split("\t")
            gene_ids[int(d[1])] = d[0]
    return gene_ids

# def score(adj, gene_ids):
#     labels = SpectralClusterting().fit_predict(adj)
#     clusters = {}
#     for i in range(0, labels.size):
#         try:
#             clusters[labels[i]] = clusters[labels[i]].append(gene_ids[i])
#         except:
#             clusters[labels[i]] = []
#     scores = []
#     for k,v in clusters.items():
#         enr = gp.enrichr(gene_list=v)
#         scores.append(enr)
#     return scores

# def step(adjs, wvec, dwvec):
#     agg = np.full_like(adjs[0], (w0[-1] + dw0[-1]))
#     for i in range(0, len(adj)):
#         agg = np.add(agg, np.multiply(adj[i], (wvec[i] + dwvec[i])))
#     return agg

def learn(adjs, wvec, gene_ids):
    stepsize = .1
    while(true):
        maxscore = 0
        maxvec = wvec
        for i in range(0, len(wvec)):
            downvec = wvec[:][i] -= stepsize
            downstep = step(adjs, wvec, downvec)
            downscore = score(downstep, gene_ids)
            upvec = wvec[:][i] += stepsize
            upstep = step(adjs, wvec, upvec)
            upscore = score(upstep, gene_ids)

            new_max = max(maxscore, downscore, upscore) == downscore)
            if (new_max == downscore):
                maxscore = downscore
                maxvec = downvec
            elif (new_max == upscore):
                maxscore = upscore
                maxvec = upvec

        if (wvec == maxvec):
            break
        else:
            wvec = max_vec

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
    graphs = load_data(DATA)
    # gene_ids = gene_id_dict(GENEIDS)
    # wvec = np.zeros(7)
    # wvec = learn(adjs, wvec, gene_ids)
    # print(wvec)
    show_eda_table(graphs)