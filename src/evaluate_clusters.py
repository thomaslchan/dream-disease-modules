import sys, time, os, pickle
from gseapy import enrichr

PATH = '../converted_results/'


def extract_below_pval(df, cutoff):
	pval_col = 'Adjusted P-value'
	return df[df[pval_col] < 0.05].sort_values(by=pval_col)


def check_term_coverage(frac, max_val=1000):
	return int(frac.split('/')[1]) <= max_val


def evaluate(clusters, pval_cutoff=0.05):
	"""
	Calculates enrichment score for modules of genes by using a p-value cutoff
	for statistical significance and using Enrichr for multiply hypothesis
	correction. This counts a module as enriched if it is 3 genes or larger,
	and if it 

	Args:
	-----------------------------------------------------------------
	- clusters: List of clusters, where each cluster is a list of gene names
	- pval_cutoff: Float for p-value cutoff
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
	for filename in sorted(os.listdir(PATH))[4:]:
		print(filename)
		clusters = pickle.load(open(os.path.join(PATH, filename), 'rb'))
		enriched, total, score = evaluate(clusters)
		print ("Enriched Clusters, Score: ", enriched, score)
