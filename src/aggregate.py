from __future__ import division
import sys, os, operator, time, math
import networkx as nx
import numpy as np
import mpmath as mp

from numpy.linalg import inv
from scipy.spatial.distance import pdist, squareform
from sklearn.cluster import KMeans
from prettytable import PrettyTable

DSDs = '../data/networks/DSDs'
NODELISTS = '../data/nodelists'
GENEFILE = '../data/ids/gene_ids.txt'
GOFILE = '../data/go/go_terms.txt'

def load_networks(path):
    adjs = []
    for filename in sorted(os.listdir(path)):
        start = time.time()
        adj = np.load(os.path.join(path, filename))
        adjs.append(adj)
        end = time.time()
        print("Loaded network " + filename + " in "
                + str(float(end-start)) + " seconds.")
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

def resize_networks(networks, nodelists):
    len_networks = len(networks)
    adjs = []
    for i in range(0, len_networks):
        start = time.time()
        len_nodelist = len(nodelists[i])
        adj = np.zeros((21115, 21115))
        for j in range(0, len_nodelist):
            for k in range(0, len_nodelist):
                adj[j][k] = networks[i][j][k]
        end = time.time()
        print("Resized network " + str(i+1) + " in "
                + str(float(end-start)) + " seconds.")
        adjs.append(adj)
    return adjs


def aggregate(adjs, wvec):
    agg = np.full_like(adjs[0], (wvec[-1]))
    for i in range(0, len(adjs)):
        start = time.time()
        agg = np.add(agg, np.multiply(adjs[i], (wvec[i])))
        end = time.time()
        print("Aggregated network " + str(i+1) + " in "
                + str(float(end-start)) + " seconds.")
    return agg

def gcn(agg, go_terms, gene_ids):
    start = time.time()
    end = time.time()
    print("Condensed aggregated network in "
                + str(float(end-start)) + " seconds.")
    return agg

def cluster(k, adj, gene_ids):
    start = time.time()
    clustering = KMeans(n_clusters=k)
    labels = clustering.fit(adj).labels_
    clusters = d = [[] for x in xrange(k)]
    for i in range(0, labels.size):
           clusters[labels[i]].append(gene_ids[i])
    end = time.time()
    print("Clustered aggregated network in "
                + str(float(end-start)) + " seconds.")
    return clusters

def sp(x, y):
    ret = 1
    while x < y:
        x += 1
        ret *= x
    return ret

def calc_pval(w, x, y, z):
    w = mp.mpf(w)
    x = mp.mpf(x)
    y = mp.mpf(y)
    z = mp.mpf(z)
    return (sp(z, y) * (sp (x, y) / float((sp (0, z)) * (sp (w-x, w-x+z)))))

def score(pval_cutoff, clusters, go_terms):
    num_genes = 21115
    score = 0
    num_clusters = len(clusters)
    for i in range(0, num_clusters):
        start = time.time()
        score_dict = {}
        cluster_size = len(clusters[i])
        for j in range(0, cluster_size):
            try:
                go = go_terms[clusters[i][j]]
            except KeyError:
                go = []
            for k in range (0, len(go)):
                try:
                    num = score_dict[go[k][0]][0]
                except KeyError:
                    score_dict[go[k][0]] = [0, 0]
                    num = 0
                score_dict[go[k][0]] = [num + 1 , go[k][1]]
        score_vec = list(score_dict.items())
        len_score_vec = len(score_vec)
        pval_vec = []
        for j in range(0, len(score_vec)):
            in_cluster = score_vec[j][1][0]
            in_total = score_vec[j][1][1]
            pval = calc_pval(num_genes, cluster_size, in_total, in_cluster)
            if (pval < pval_cutoff/float(len_score_vec)):
                score += 1
                break
        end = time.time()
        print("Scored cluster " + str(i+1) + " in "
                + str(float(end-start)) + " seconds.")
    return score/float(num_clusters)

def to_file(clusters, scoreVal):
    with open("output.txt", "w+") as file:
	file.write("{}\n".format(scoreVal))
	for i in range(0, len(clusters)):
	    file.write("Cluster {}:\n".format(i+1))
	    for g in cluster[i]:
		file.write("{}\t".format(g))
	    file.write("\n")
	

if __name__ == '__main__':
    networks = load_networks(DSDs)
    nodelists = read_nodelists(NODELISTS)
    gene_ids = gene_id_dict(GENEFILE)
    go_terms = go_term_dict(GOFILE)

    adjs = resize_networks(networks, nodelists)

    wvec = [1, 1, 1, 1, 1, 1, 0]
    k = 2
    pval_cutoff = .05
    agg = aggregate(adjs, wvec)
    embedding = gcn(agg, go_terms, gene_ids)
    clusters = cluster(k, embedding, gene_ids)
    scoreVal = score(pval_cutoff, clusters, go_terms)


"""

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
"""
