# COMUNET
There are following tutorials available:

**clustering_patternSearch:**
1. 1_PrepareData.ipynb
2. 2_CellPhoneDBanalysis.ipynb
3. 3_Clustering_PatternSearch.ipynb

**comparativeAnalysis:**
1. 1_PrepareData.ipynb
2. 2_CellPhoneDBanalysis.ipynb
3. 3_ComparativeAnalysis.ipynb

## clustering_patternSearch tutorial
These tutorials use mouse embryo data.
### 1_PrepareData.ipynb notebook
This notebook uses R code, please make sure an R kernel is installed for your jupiter notebook or copy the code and run it in an R environment.

This notebook describes the steps for:
* downloading the data
* cleaning the data (subset for E6.5 stage)
* filtering of cells (cell types that are represented by less than 5 cells are filtered out)
* subsetting the dataset randomly for 1500 cells
* normalising data with scran
* filtering genes (genes with mean expression less than 1 in all expressing cells are filtered out)
* generate UMAP of the data

The **output** of this notebook is:
* counts_filtered_norm.txt a normalised (NOT log2 transformed) filtered count matrix (gene names in rows, cell IDs in columns)
* anno_filtered.txt an annotation matrix (containing columns "Cell" and "cell_type").

### 2_CellPhoneDBanalysis.ipynb notebook
This notebook uses CellPhoneDB venv, please make sure the venv kernel is installed for you jupiter notebook or run CellPhoneDB from the command line as described in the (CellPhoneDB tutorial)[https://github.com/Teichlab/cellphonedb].

This notebook describes the steps for:
* installation of CellPhoneDB
* setting parameters for CellPhoneDB (gene_input_type = "gene_name", subsample = False) 
* running CellPhoneDB

The **input** of this notebook is:
* counts_filtered_norm.txt a normalised (NOT log2 transformed) filtered count matrix (gene names in rows, cell IDs in columns)
* anno_filtered.txt an annotation matrix (containing columns "Cell" and "cell_type").

The **output** of this notebook is a "cpdb_output" folder that contains following files:
* deconvoluted.txt
*  means.txt
* pvalues.txt
* significant_means.txt

**PLEASE NOTE** that CellPhoneDB analysis might give slightly different results when rerun on the same data, which might affect the downstream analysis. If you wish to reproduce the exact figures from the manuscript, please run the downstream analysis on cpdb_output_embryo output.

### 3_Clustering_PatternSearch.ipynb notebook
This notebook uses R code, please make sure an R kernel is installed for your jupiter notebook or copy the code and run it in an R environment.

This notebook describes the steps for:
* defining parameters for COMUNET:
	* export_results = TRUE the results will be exported into "comunet_output" folder
	* minClusterSize = 6 not less than 6 ligand-receptor pairs are allowed per cluster
* reading in input flies
* converting CellPhoneDB output (significant_means.txt  file) into a weight array
* running Clusterins analysis of COMUNET
* running Pattern Search analysis of COMUNET

The **input** for this notebook is 
* significant_means.txt from "cpdb_output" folder 
* counts_filtered_norm.txt a normalised (NOT log2 transformed) filtered count matrix (gene names in rows, cell IDs in columns)
* anno_filtered.txt an annotation matrix (containing columns "Cell" and "cell_type").
* complex_input.csv obtained from CellPhoneDB site
* gene_input.csv obtained from CellPhoneDB site

The **output** of this notebook is comunet_output folder with following files:

*for clustering analysis:*
* result.RData an list that contains:
	* cellNumbers: number of cells per cell type
	* nodes: cell types
	* ligandReceptorPairs: dataframe of ligand-receptor pairs with columns: 
		* "pair" contains values in a form "ligand:receptor", i.e. ligand being at the first place, receptor being at the second place, e.g. "TNFSF13:TNFRSF17"
		*  "ligand" contains ligand names, e.g. "TNFSF13"
		* "ligand_complex_composition" if ligand is a complex (e.g. "aXb2_complex"), contains genes in the ligand complex separated with a comma, e.g. "ITGAX,ITGB2", else contains ""
		* "receptor" contains receptor names, e.g. "TNFRSF17"
		* "receptor_complex_composition" if receptor is a complex (e.g. "NKG2D_II_receptor"), contains genes in the receptor complex separated with a comma, e.g. "KLRK1,HCST", else contains ""
	* weight_array: array of weighted adjacency matrices (#nodes x #nodes x #ligand-receptor pairs)
	* degree_array: array of node degrees (#nodes x 3 (in, out , delta) x # ligand-receptor pairs)
	* dissim_matrix: dissimilarity matrix
	* clusters: cluster assignment for each legnd-receptor pair
	* weight_array_byCluster: array of weighted adjacency matrices of clusters (#nodes x #nodes x #clusters). Each weighted adjacency matrix contains average edge weights among ligand-receptor pairs in a cluster.
	* degree_array_byCluster: array of node degrees (#nodes x 3 (in, out , delta) x # clusters). Each each degree matrix contains average -in , -out, and -delta node degree values among ligand-receptor pairs in a cluster.
* cellNumbers.csv: number of cells per cell type
* cluster.plots.pdf: in a single pdf for each cluster a cluster communication graph and a graph of ligands and receptors that belong to this cluster
* clusterX.pdf: one .pdf per cluster (i.e. cluster1.pdf for cluster1). One .pdf contains for each ligand-receptor pair that belong to this cluster a ligand-receptor pair communication graph and log10 expression levels for ligand(s) and receptor(s) that belong to this ligand-receptor pair.
* clusters.csv: a table of ligand-receptor pairs and their corresponding cluster number.
* overview.pdf: one .pdf that contains a heat map of all ligand-receptor pairs sorted by cluster, a UMAP of all ligand-receptor pairs collared by cluster.
* unclustered_LRPs.pdf: one .pdf that contains for each ligand-receptor pair that doesn't belong to any cluster a ligand-receptor pair communication graph and log10 expression levels for ligand(s) and receptor(s) that belong to this ligand-receptor pair.
	
*for pattern search analysis:*
* pattern_results.csv: a table with all lingand-receptor pairs with dissimilarity < 1 to the pattern of interest sorted by increasing dissimilarity.
* pattern.pdf: one .pdf that contains for each ligand-receptor pair from pattern_results.csv (sorted by increasing dissimilarity) a ligand-receptor pair communication graph and log10 expression levels for ligand(s) and receptor(s) that belong to this ligand-receptor pair.

## comparativeAnalysis tutoarial
These tutorials use human AML data.
### 1_PrepareData.ipynb notebook
This notebook uses R code, please make sure an R kernel is installed for your jupiter notebook or copy the code and run it in an R environment.

This notebook describes the steps for:
* downloading the data
* filtering of cells (cell types that are represented by less than 5 cells are filtered out)
* filtering of cell types (cell types that are present only in one condition, but not in the other are filtered out)
* merge datasets for joint normalisation 
* normalising data with scran
* splitting data after normalisation
* filtering lowly expressed genes (genes with mean expression less than 1 in all expressing cells are filtered out)
* filtering unshared genes (genes that are present only in one condition, but not in the other are filtered out)
* generate UMAP of the data

The **output** of this notebook  are two folders:
* AML328_d0_forComparative
* AML328_d29_forComparative

each of them containing following files:
* counts_filtered_norm.txt a normalised (NOT log2 transformed) filtered count matrix (gene names in rows, cell IDs in columns)
* anno_filtered.txt an annotation matrix (containing columns "Cell" and "cell_type").

### 2_CellPhoneDBanalysis.ipynb notebook
THIS NOTEBOOK SHOULD BE RUN SEPARATELY FOR AML328_d0 AND AML328_d29 data!

This notebook uses CellPhoneDB venv, please make sure the venv kernel is installed for you jupiter notebook or run CellPhoneDB from the command line as described in the CellPhoneDB tutorial.

This notebook describes the steps for:
* installation of CellPhoneDB
* setting parameters for CellPhoneDB (gene_input_type = "gene_name", subsample = False) 
* running CellPhoneDB

The **input** of this notebook is:
* counts_filtered_norm.txt a normalised (NOT log2 transformed) filtered count matrix (gene names in rows, cell IDs in columns)
* anno_filtered.txt an annotation matrix (containing columns "Cell" and "cell_type").

The **output** of this notebook is a "cpdb_output" folder that contains following files:
* deconvoluted.txt
* means.txt
* pvalues.txt
* significant_means.txt

**PLEASE NOTE** that CellPhoneDB analysis might give slightly different results when rerun on the same data, which might affect the downstream analysis. If you wish to reproduce the exact figures from the manuscript, please run the downstream analysis on **cpdb_output_AML328_d0** and **cpdb_output_AML328_d29** outputs.

### 3_ComparativeAnalysis.ipynb notebook
This notebook uses R code, please make sure an R kernel is installed for your jupiter notebook or copy the code and run it in an R environment.

This notebook describes the steps for:
* defining parameters for COMUNET:
	* export_results = TRUE the results will be exported into "comunet_output" folder
	* minClusterSize = 6 not less than 6 ligand-receptor pairs are allowed per cluster
* reading in input flies (for both conditions)
* converting CellPhoneDB output (significant_means.txt  file) into a weight array (for both conditions)
* running Clusterins analysis of COMUNET (for both conditions)
* running Comparative analysis of COMUNET

The **input** for this notebook is 
* significant_means.txt from "cpdb_output" folder (for both conditions)
* counts_filtered_norm.txt a normalised (NOT log2 transformed) filtered count matrix (gene names in rows, cell IDs in columns) (for both conditions)
* anno_filtered.txt an annotation matrix (containing columns "Cell" and "cell_type") (for both conditions)
* complex_input.csv obtained from CellPhoneDB site
* gene_input.csv obtained from CellPhoneDB site

The **output** of this notebook is comunet_output folder with following files:
* comparative_dissimMatrix.csv: a pairwise dissimilarity matrix between ligand receptor pairs of two conditions
* comparative_overview.pdf: a .pdf file that contains an euler diagram of ligand-receptor pairs of both conditions, one clustered and one unclustered heatmap of pairwise dissimilarity of ligand-receptor pairs from both conditions 
* comparative_plots.pdf:  a .pdf that contains for each ligand-receptor of both conditions a ligand-receptor pair communication graph and log10 expression levels for ligand(s) and receptor(s) that belong to this ligand-receptor pair.
* comparative_sortedLRPlist.csv: a table of ligand-receptor pairs sorted by increasing dissimilarity between two conditions
