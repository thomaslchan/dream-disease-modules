import time, os
import numpy as np
import networkx as nx
from prettytable import PrettyTable

GRAPHS = '../data/networks/anonymized'
DSDs = '../data/networks/DSDs'
NODELISTS = '../data/nodelists'
NODES = 'nodes.pkl'
GENEFILE = '../data/ids/gene_ids.txt'


def read_graphs(path=GRAPHS):
    """
    Reads in graphs from .txt files that have weighted edges in 
    [node] [node] [weight] format, and returns a list of Networkx graphs.

    Args:
    -----------------------------------------------------------------
    - path: Directory with .txt files of graphs (default=GRAPHS)
    """
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


def get_nodelists(graphs):
    """
    Get node-lists from graphs in sorted order and return a list of 
    list of integers.

    Args:
    -----------------------------------------------------------------
    - graphs: List of Networkx graphs
    """
    sgs = [max(nx.connected_component_subgraphs(G), key=len) 
                if not nx.is_directed(G)
                else max(nx.weakly_connected_component_subgraphs(G), key=len)
                    for G in graphs]
    nodelists = [sorted(map(int, sg.nodes)) for sg in sgs]
    return nodelists


def load_DSDs(path=DSDs):
    """
    Loads distance matrices from a folder that contains .npy files.
    Returns list of distance matrices as Numpy arrays.

    Args:
    ---------------------------------------------------------------
    - path: Relative path to directory with .npy files
    """ 
    start = time.time()
    matrices = []
    for filename in sorted(os.listdir(path)):
        matrix = np.load(os.path.join(path, filename))
        matrices.append(matrix)
    end = time.time()
    print("Loaded networks in " + str(end-start) + " seconds")
    return matrices


def describe_dsds(matrices):
    """Prints distribution information about DSD matrices (e.g. quartiles)."""
    rows = [[i+1, np.mean(m), np.std(m),
                  np.percentile(m, 25), np.percentile(m, 75)] 
                        for (i, m) in enumerate(matrices)]
    t = PrettyTable()
    t.field_names = ['Network #', 'Mean', 'Std. Dev.', 
                    'First Quartile', 'Third Quartile']
    [t.add_row(row) for row in rows]
    print(t)
    

def read_nodelists(path=NODELISTS):
    """
    Loads in node-lists from a folder that contains .txt files with the
    indices of the nodes in the network's largest connected component.
    Returns a list of lists. 

    Args:
    ---------------------------------------------------------------
    - path: Relative path to directory with .txt files
    """
    nodelists = []
    for filename in sorted(os.listdir(path)):
        with open(os.path.join(path, filename), "r") as file:
            nodelist = []
            for line in file.readlines():
                nodelist.append(int(line.strip()))
            nodelists.append(nodelist)
    return nodelists


def load_nodelists(path=NODES):
    """
    Loads in nodelists from a pickle file, and returns a list of lists.
    
    Args:
    -------------------------------------------------------------
    - path: Path to .pkl file with nodelists
    """
    with open(path, 'rb') as f:
        nodelists = pickle.load(f)
    return nodelists


def gene_id_dict(filename=GENEFILE):
    """
    Loads gene IDs from a .txt file. Returns dictionary with node index
    as key, and the HGNC gene symbol as value.

    Args:
    ---------------------------------------------------------------
    - filename: Path to .txt file with gene IDs.
    """
    gene_ids = {}
    with open(filename, "r") as file:
        for line in file.readlines():
            d = line.strip().split("\t")
            gene_ids[int(d[1])] = d[0]
    return gene_ids


if __name__ == '__main__':
    dsds = load_DSDs()
    describe_dsds(dsds)