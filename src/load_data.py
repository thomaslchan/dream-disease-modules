import time, os
import numpy as np


DSDs = '../data/networks/DSDs'
NODELISTS = '../data/nodelists'
GENEFILE = '../data/ids/gene_ids.txt'


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
