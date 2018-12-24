# Comparing Imputation Techniques for Multigraph Clustering -- The DREAM Challenge.

## Scripts:

```
python aggregate.py
```

Finds the scored aggregated clustering of a collection of np matrices
located in data/networks/DSDs, using node lists in data/nodelists and
gene_ids in data/ids. Takes imputation method and number of clusters as
optional commandline arguments with flags "--impute" and "--n_clusters",
respectively.

```
python dsd.py
```

Computes the DSD matrices and nodelists of a collection of edgelists
located in data/networks/anonymized. Code sourced from Tufts University's
Diffusion State Distance research team, led by Prof. Lenore Cowen. 

```
python evaluate_clusters.py
```

Evaluates all clusterings for GO enrichment in converted_results folder.

```
bash launch_experiment.sh submit
```

Runs all the experiments on the Tufts High Performance Computing (HPC) Cluster. 

## Results

| Imputation Method | Enriched | Total | Score |
|:-----------------:|:--------:|:-----:|:------:
| Max               |  120    | 250   | 48%    |
| Min               | 143 | 250  | 57.2% |
| Mean              | 134 | 250  | 53.6% |
| Median            | 150 | 250  | 60% |
| Max Local         | 98 | 250  | 39.2% |
| Min Local         | 112 | 250  | 44.8% |
| Mean Local		| 98 | 250  | 39.2% |
| Median Local      | 98 | 250  | 39.2% |

## Acknowledgments

Thanks to Prof. Lenore Cowen from Tufts University for supervising us during
this work. We would also like to acknowledge Christian Zinck, who contributed
to earlier versions of this code base and helped produce some of the cleaned
data used in the project.
