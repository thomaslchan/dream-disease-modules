from __future__ import division
import sys, os, operator, time

import networkx as nx
import numpy as np
from numpy.linalg import inv
from scipy.spatial.distance import pdist, squareform
from prettytable import PrettyTable

DATA = 'projdata/1_networks/original'


def read_data(path):
    graphs = []
    i = 0
    for filename in sorted(os.listdir(path)):
        print(filename)
        start = time.time()
        i += 1
        g_type = nx.OrderedDiGraph if i == 3 else nx.OrderedGraph
        G = nx.read_weighted_edgelist(os.path.join(path, filename),
                                    create_using=g_type)
        graphs.append(G)
        end = time.time()
        print("Reading in network " + str(i) + " took " 
                + str(int(end-start)) + " seconds.\n\n")
    return graphs


def build_transition_matrix(adjacency_graph):
    degs = adjacency_graph.sum(axis=1)

    transition = adjacency_graph / degs[:,None]
    transition[degs==0,:] = 0
    transition[degs==0,degs==0] = 1
    return transition


def calc_hescotts(transition, iters, v=True, n=None):
    # params v and n are not used, left for compatibility with original dsdcore script
    nRw = iters
    p = transition
    n = p.shape[0]
    c = np.eye(n)
    c0 = np.eye(n)
    for i in xrange(nRw):
        c = np.dot(c, p) + c0
    return c


def calc_dsd(transition):
    return squareform(pdist(transition,metric='cityblock'))              


def add_self_edges(adjacency_graph, base_weight=1):    
    n = np.size(adjacency_graph[0])
    ident = np.identity(n)*base_weight
    return np.add(adjacency_graph,ident)


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
        if data['weight'] < threhold:
            G.remove_edge(u,v)
    return G


if __name__ == '__main__':
    graphs = read_data(DATA)
    show_eda_table(graphs)
