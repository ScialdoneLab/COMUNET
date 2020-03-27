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
The 1_PrepareData.ipynb describes the stepst of data processing from the raw data to normalized filtered data. 

**[Tutorial_embryo_patternSearch](https://github.com/ScialdoneLab/COMUNET/tree/COUMNET-dev/Tutorials/Tutorial_embryo_patternSearch):**
1. 1_PrepareData.ipynb
2. 2_CellPhoneDBanalysis.ipynb
3. 3_PatternSearch.ipynb
These tutorials use mouse embryo data from [Nowotschin *et al.*, 2019](https://www.nature.com/articles/s41586-019-1127-1).

**[Tutorial_AML_comparative](https://github.com/ScialdoneLab/COMUNET/tree/COUMNET-dev/Tutorials/Tutorial_AML_comparative):**
1. 1_PrepareData.ipynb
2. 2_CellPhoneDBanalysis.ipynb
3. 3_ComparativeAnalysis.ipynb
These tutorials use human AML data from [van Galen *et al.,* 2019](https://doi.org/10.1016/j.cell.2019.01.031).
