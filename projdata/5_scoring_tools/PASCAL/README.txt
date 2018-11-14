=============================================================================================
These data have been shared confidentially with participants of the Disease Module
Identification DREAM Challenge.

ALL DATA AND RESULTS ARE EMBARGOED UNTIL PUBLICATION OF THE MAIN CHALLENGE PAPER.
See the Challenge website for further information and data access conditions.

Challenge website: https://www.synapse.org/modulechallenge
=============================================================================================


This is the version of the Pascal tool that was used for the Disease Module Identification
DREAM11 Challenge. Pascal is an efficient tool for computation of gene and pathway scores from
GWAS SNP p-values. See also:

- Website: http://www2.unil.ch/cbg/index.php?title=Pascal

- Reference: Lamparter D, Marbach D, Rueedi R, Kutalik Z, and Bergmann S. Fast and rigorous
computation of gene and pathway scores from SNP-based summary statistics. PLoS Computational
Biology 12, e1004714, 2016.

Here, we describe how to proceed to score module predictions using the Pascal tool as done in
the challenge. We assume that you have read the Pascal paper and user manual included here so
that you are familiar with the method and the tool.


0. Requirements
---------------

In order to reporduce the results of the challenge exactly, Java 8 is required. The enrichment
p-values can vary slightly when using previous Java versions due to different random ordering
of genes with identical scores.


1. Install Pascal
-----------------

Follow the instructions in the user manual to install Pascal. Note that if installation of the
native libraries fails, Pascal can still be run but it will be much slower. See Section 7 of 
the user manual.

If the native libraries were successfully installed, you will see this message when running
Pascal (see next step):
	INFO: successfully loaded /tmp/jniloader7434279127983874993netlib-native_system-linux-x86_64.so
	INFO: already loaded netlib-native_system-linux-x86_64.so

If the native libraries could not be loaded, you will see this warning when running Pascal:
	WARNING: Failed to load implementation from: com.github.fommil.netlib.NativeSystemBLAS
	WARNING: Failed to load implementation from: com.github.fommil.netlib.NativeRefBLAS
	WARNING: Failed to load implementation from: com.github.fommil.netlib.NativeSystemLAPACK
	WARNING: Failed to load implementation from: com.github.fommil.netlib.NativeRefLAPACK

The result is identical with or without native libraries, only the runtime is affected.


2. Run Pascal
-------------

We included an example module prediction and GWAS data in the directory dream11_example.
The modules are the prediction for network 1 from the best performer. The GWAS trait is
total cholesterol. Note that there is both a file for the SNP p-values and the pre-computed
gene scores. Pascal uses pre-computed gene scores for independent genes of a module. However,
genes of the same module that are proximal on the genome are not independent due to linkage
disequilibrium and must be treated as a unit to take the full correlation structure of the
locus into account (see the Pascal paper). For this reason, we must also specify the file with
the SNP p-values.

From the terminal, go to PASCAL directory and run the following command:

./Pascal --set=dream11_settings/settings_leaderboard_v3--1_ppi_anonym_v2.txt --runpathway=on --genescoring=sum\
 --pval=dream11_example/EUR.GLGC.jointGwasMc_TC.txt.gz\
 --genescorefile=dream11_example/EUR.GLGC.jointGwasMc_TC.sum.genescores.txt.gz\
 --genesetfile dream11_example/3344522.7320912.1_ppi_anonym_v2.txt.gz\
 --outdir .

The options that we use are:

--set=dream11_settings/settings_leaderboard_v3--1_ppi_anonym_v2.txt

This specifies the settings file. The settings file defines all the options of Pascal that are
not given on the command line. Note that we use a different settings file for each network
of the challenge. The only difference between these settings files is the gene annotation file.
(In retrospect, it would have been more practical to add a command-line option for this). For
network 1, we specify the corresponding annotation file, which includes only the genes of
network 1:
	ucscAnnotationFile = resources/annotation/ucsc/ucsc_known_genes_2013-09-03--1_ppi_string.txt
	backgroundAnnotationFile = resources/annotation/ucsc/ucsc_known_genes_2013-09-03--1_ppi_string.txt

This ensures that only the genes of the network are used as background to compute enrichment.
If your module predictions include all genes of the network, you may also use the full
annotation and set the flag onlyPathwayGenesAsBackground:
	ucscAnnotationFile = resources/annotation/ucsc/ucsc_known_genes_2013-09-03.txt
	backgroundAnnotationFile = 
	onlyPathwayGenesAsBackground = 0

--runpathway=on

Tells Pascal to run the pathway (module) enrichment analysis

 --genescoring=sum
 
Use the sum of chi-squared statistic for the gene scoring

--pval=dream11_example/EUR.GLGC.jointGwasMc_TC.txt.gz

The GWAS SNP p-values.

--genescorefile=dream11_example/EUR.GLGC.jointGwasMc_TC.sum.genescores.txt.gz

The precomputed GWAS gene scores. IMPORTANT: if you change any of the settings regarding the
gene scoring in the settings file, you have to recompute the gene scores. (Otherwise,
different settings would have been used for the scoring of the precomputed individual genes
and fused/meta genes that are computed dynamically).

--genesetfile dream11_example/3344522.7320912.1_ppi_anonym_v2.txt.gz

The module predictions / pathways given in the format used in the challenge (using gene
symbols, not the anonymized gene IDs of the challenge).

--outdir 

You can specify a custom output directory.


3. Output files
---------------

Pascal writes a couple of files with logging information that you can ignore. The relevant
file with the module enrichment p-values is:
- EUR.GLGC.jointGwasMc_TC.PathwaySet--3344522.7320912.1_ppi_anonym_v2--sum.txt

Your output file should be identical to the one provided in dream11_examples. Note that these
are the nominal p-values before any multiple testing correction.


4. Run Pascal using all GWAS datasets of the challenge
------------------------------------------------------

Download the file 4_gwas_datasets.zip, which includes the SNP p-values and pre-computed gene
scores for over 200 GWAS datasets. Note that not all GWASs were used in the challenge, some
were excluded, e.g. because they showed no signal in exploratory analyses. The GWAS sets used
for the leaderboard phase and the final evaluation are listed in two files that are included
in the same zip file. 

To score a module prediction for one network, you have to run Pascal for each of the GWAS
datasets (i.e., 180 jobs for the leaderboard and hold-out GWAS sets). This is only practical
on a computing cluster. Due to differences between queuing systmes, you will have to write
your own scripts to launch the jobs on your cluster. We recommend that you write all output
files for a given module prediction in the same output directory (using the option --outdir).


5. Compute final score
----------------------

See the README in the directory Scoring_scripts.


--
Daniel Marbach, Sarvenaz Choobdar, Sven Bergmann
November 3, 2016
