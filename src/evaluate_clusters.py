import sys, time, os, pickle
import gseapy as gp

def score(clusters, pval_cutoff=0.05):
	"""
	Calculates score for clustering of genes by statistical significance of
	given p value cutoff using gseapy and returns percentage of enriched
	clusters

	Args:
	-----------------------------------------------------------------
	- clusters: List of clusters
	- pval_cutoff: Integer for p value cutoff
	"""
	start = time.time()
	numerator = 0
	denominator = 0
	for cluster in clusters:
		if len(cluster) > 2:
			try:
				df = gp.enrichr(gene_list=cluster,
					gene_sets=['GO_Biological_Process_2018',
					'GO_Molecular_Function_2018'], cutoff=10).results

				for row in df[['Overlap', 'Adjusted P-value']].iterrows():
					if (row['Adjusted P-value'] < 0.05 and
					   int(row['Overlap'].split('/')[1]) < 100):
						numerator += 1
						break
			except:
				pass
			denominator += 1
	end = time.time()
	print("Scored clusters in " + str(end - start) + " seconds")
	
	if denominator == 0:
		return denominator
	else:
		return numerator / denominator

if __name__ == '__main__':
	path = sys.argv[1]
	for filename in sorted(os.listdir(path)):
		start = time.time()
		print("Testing")
		print(filename)
		clusters = pickle.load(open(os.path.join(path, filename), 'rb'))
		print ("Score: " + str(score(clusters)))
