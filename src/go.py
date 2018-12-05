import sys
if __name__ == '__main__':
    go_dict = {}
    with open(sys.argv[1], "r") as f:
	for line in f.readlines():
	    l = line.strip().split('\t')
	    go_dict[l[0][3:]] = l[2].split(" ")
    gene_dict = {}
    m = 0
    for k,v in go_dict.items():
	length = len(v)
	for i in range(0, len(v)):
            temp = []
	    try:
		temp = gene_dict[v[i]]
	    except KeyError:
		gene_dict[v[i]] = []
	    temp.append([k, length])
	    gene_dict[v[i]] = temp
    with open("output.txt", "w") as f:
	for k,v in gene_dict.items():
	    f.write("{}\t".format(k))
	    for i in range(0, len(v)):
		f.write("{},{} ".format(v[i][0], v[i][1]))
	    f.write("\n")
