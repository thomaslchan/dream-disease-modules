=============================================================================================
These data have been shared confidentially with participants of the Disease Module
Identification DREAM Challenge.

ALL DATA AND RESULTS ARE EMBARGOED UNTIL PUBLICATION OF THE MAIN CHALLENGE PAPER.
See the Challenge website for further information and data access conditions.

Challenge website: https://www.synapse.org/modulechallenge
=============================================================================================


We provide here a simple python script demonstrating how the final scores can be computed from
the Pascal output files. The procedure is:

	significant_modules = empty set

	For every GWAS g {
	   Load the nominal module enrichment p-values from the Pascal output file
	   Remove modules with NA values (see below for comment)
	   Adjust p-values using the Benjamini-Hochberg procedure (after removing NAs!)
	   significant_modules_g = (modules with adjusted p-value < 0.05)
	   significant_modules = union(significant_modules, significant_modules_g)
	 }

     score = count(significant_modules)

The score should be computed for each of the six networks. The overall score of the challenge 
is the sum of the scores across the six networks.

Run the python script with the provided test data using the following command:

>> python multiple_testing.py example_pascal_output/3344522.7320912.1_ppi_anonym_v2/ 3344522.7320912

Output:
>> Multiple testing correction of PASCAL output for genesetfile: 3344522.7320912
>> Number of Significant modules: 16

Which corresponds to the score of the best performing team (3344522.7320912) on network 1 using
the hold-out GWAS set (you can change the GWAS set in the script).



--
Daniel Marbach, Sarvenaz Choobdar, Sven Bergmann
November 3, 2016
