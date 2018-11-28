import sys
import networkx as nx

def main(argv):
	G = nx.MultiGraph()
        edges = []
        with open(argv[1], 'r') as f:
                for line in f.readlines():
                        edges.append(line.strip("\n\r").split("\t"))
	G.add_weighted_edges_from(edges)

	with open("output.txt", 'w') as f:
		for e in G.edges():
			weights = G[e[0]][e[1]]
			sum = 0.0
			num = 0.0
			for k, v in weights.items():
				sum += float(v['weight'])
				num += 1
			f.write("{}\t{}\t{}\n".format(e[0], e[1], sum/num))
if __name__ == '__main__':
        main(sys.argv)
