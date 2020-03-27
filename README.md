# COMUNET

Intercellular communication plays an essential role in multicellular organisms and several algorithms to analyse it from single-cell transcriptional data have been recently published, but the results are often hard to visualize and interpret. 
We developed COMUNET (Cell cOMmunication exploration with MUltiplex NETworks), a tool that streamlines the interpretation of the results from cell-cell communication analyses. 
COMUNET uses multiplex networks to represent and cluster all potential communication pathways between cell types. The algorithm also enables the search for specific patterns of communication and can perform comparative analysis between two biological conditions. To exemplify its use, here we apply COMUNET to investigate cell communication patterns in single-cell transcriptomic datasets from mouse embryos and from an acute myeloid leukemia patient at diagnosis and after treatment. 

There are following tutorials available:

**[Tutorial_embryo_lrpClustering](https://github.com/ScialdoneLab/COMUNET/tree/COUMNET-dev/Tutorials/Tutorial_embryo_lrpClustering):**
1. 1_PrepareData.ipynb
2. 2_CellPhoneDBanalysis.ipynb
3. 3_LRP_clustering.ipynb

These tutorials use mouse embryo data from [Nowotschin *et al.*, 2019](https://www.nature.com/articles/s41586-019-1127-1).

The 1_PrepareData.ipynb tutorial describes the stepst of data processing from the raw counts and annotation to normalized filtered counts and filtered annotation. If you prefere to skip this step, you can find the normalized filtered counts and filtered annotation in the tutorial subfolders AML328_d0_forComparative and AML328_d29_forComparative.

The 2_CellPhoneDBanalysis.ipynb tutorial descripbes the steps to run CellPhoneDB analysis using the normalized filtered counts and filtered annotation. If you prefere to skip this step, you can find the results in the cpdb_output subfolders of the AML328_d0_forComparative and AML328_d29_forComparative. Please note that COMUNET does not necessarily rely on the results of CellPhoneDB algorithm. You can use any algorithm of choise, which generates weight matrices, which represent how likely/strong is the interaction between two different cell types mediated by a given pair of proteins. 

The 3_LRP_clustering.ipynb tutorial describes the COMUNET ligand-receptor clustering procedure.

**[Tutorial_embryo_patternSearch](https://github.com/ScialdoneLab/COMUNET/tree/COUMNET-dev/Tutorials/Tutorial_embryo_patternSearch):**
1. 1_PrepareData.ipynb
2. 2_CellPhoneDBanalysis.ipynb
3. 3_PatternSearch.ipynb

These tutorials use mouse embryo data from [Nowotschin *et al.*, 2019](https://www.nature.com/articles/s41586-019-1127-1).

The tutorials 1_PrepareData.ipynb and 2_CellPhoneDBanalysis.ipynb are the same, as decribed above.

The tutorial 3_PatternSearch.ipynb describes the COMUNET cell-cell communication pattern search procedure.

**[Tutorial_AML_comparative](https://github.com/ScialdoneLab/COMUNET/tree/COUMNET-dev/Tutorials/Tutorial_AML_comparative):**
1. 1_PrepareData.ipynb
2. 2_CellPhoneDBanalysis.ipynb
3. 3_ComparativeAnalysis.ipynb

These tutorials use human AML data from [van Galen *et al.,* 2019](https://doi.org/10.1016/j.cell.2019.01.031).

The tutorials 1_PrepareData.ipynb and 2_CellPhoneDBanalysis.ipynb describe equivalent procedures, as described above.

The tutorial 3_ComparativeAnalysis.ipynb describes the COMUNET comparative analysis between communication patterns in a day 0 AML sample and the sample of the same patient at day 29 after therapy.
