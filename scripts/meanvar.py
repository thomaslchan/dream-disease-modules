import sys
import numpy as np

def main(argv):
	weights = []
	with open(argv[1], 'r') as f:
		for line in f.readlines():
			w = float(line.strip("\n\r").split("\t")[2])
			weights.append(w)

	print("Mean: {}".format(np.mean(weights)))
	print("Var: {}".format(np.var(weights)))

if __name__ == '__main__':
	main(sys.argv)
