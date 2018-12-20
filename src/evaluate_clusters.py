import sys, time, os, pickle
from gseapy import enrichr

PATH = '../converted_results/'


def extract_below_pval(df, cutoff):
	pval_col = 'Adjusted P-value'
	return df[df[pval_col] < 0.05].sort_values(by=pval_col)


def check_term_coverage(frac, max_val=1000):
	return int(frac.split('/')[1]) <= max_val


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
	enriched = total = 0

	for cluster in clusters:
		if len(cluster) > 2:
			total += 1
			df = enrichr(gene_list=cluster,
				gene_sets=['GO_Biological_Process_2018',
				'GO_Molecular_Function_2018']).results
			sorted_df = extract_below_pval(df, cutoff=pval_cutoff)
			if sorted_df['Overlap'].apply(check_term_coverage).any():
				enriched += 1
		print("Current score: " + str(enriched) + '/' + str(total))
	end = time.time()
	print("Scored clusters in " + str(end - start) + " seconds")
	
	return enriched, total, enriched/total

if __name__ == '__main__':
	for filename in sorted(os.listdir(PATH)):
		print(filename)
		clusters = pickle.load(open(os.path.join(PATH, filename), 'rb'))
		enriched, total, score = score(clusters)
		print ("Enriched Clusters, Score: ", enriched, score)
