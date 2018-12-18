import time, os
from funcassociate.client import _fc
#no access to gene_ids!!

def score(clusters, pval_cutoff=0.05):
    """
    Calculates score for clustering of genes by statistical significance of
    given p value cutoff using FuncAssociate and returns percentage of enriched
    clusters

    Args:
    -----------------------------------------------------------------
    - clusters: List of clusters
    - pval_cutoff: Integer for p value cutoff
    """
    start = time.time()
    score = 0
    cluster_terms = []
    fc = _fc()
    for i in range(len(clusters)):
        try:
            cluster = clusters[i]
            ret = fc.functionate(query=cluster, species="Homo sapiens",
                        namespace="hgnc_symbol", mode="ordered",
                        cutoff=(pval_cutoff))["under"]

            # Result pval is the product of the p value and the number of
            # genes in the query that have the attribute
            # pval = ret[0][4] * (ret[0][1])
            pval = ret[0][4]

            # Cluster terms is appending the attribute ID
            cluster_terms.append(ret[0][7])
        except:
            pval = 1
            cluster_terms.append("")
        if (pval < pval_cutoff):
            score += 1
    end = time.time()
    print("Scored clusters in " + str(end - start) + " seconds")
    print(str(score/float(len(clusters))))

if __name__ == '__main__':

    path = sys.argv[1]
    for filename in sorted(os.listdir(path)):
        print(filename)
        start = time.time()
        clusters = pickle.load(open(os.path.join(path, filename), 'rb'))
        for cluster in clusters:
            score(cluster)
