import sys
import networkx as nx
import os

def read_data(path):
	graphs = []
	i = 1
	for filename in os.listdir(path):
		print (i)
		i += 1
		edges = []
		with open(os.path.join(path, filename)) as file:
                        for line in file.readlines():
			        edges.append(line.strip('\r\n').split('\t'))
		G = nx.MultiGraph()
		G.add_weighted_edges_from(edges)
		graphs.append(G)
	return graphs

def main(argv):
	graphs = read_data(argv[1])

if __name__ == '__main__':
        main(sys.argv)
