#' Comparative analysis of AML328-d0 and AML328-d31 samples.
#'
#' Dataset contains a data frame of sorted ligand-receptor pairs and
#' a matrix of pairwise dissimilarity between all ligand-receptor pairs in the two conditions.
#'
#' @format A data list of three variables
#' \describe{
#' \item{sorted_LRP_df}{
#'
#'   Data frame with columns:
#'   \itemize{
#'     \item "pair" (character string vector): names of ligand-receptor pairs in the same form as they are in ligand_receptor_pair_df$pair.
#'     \item "presence" (character string vector): whether the ligand-receptor pair is present in both conditions ("shared") or only in one of them.
#'     \item "dissimilarity" (numeric vector): dissimilarity value between the topology of ligand-receptor pair graph in two conditions.
#'     The smaller the dissimilarity value, the more similar is the graph topology between the two conditions.
#'     If a ligand-receptor pair is present only in one of the conditions, the dissimilarity is equal to 1.
#'
#'   }
#'   }
#'   \item{dissim_cond1_cond2}{
#'
#'   Numeric matrix:
#'    pairwise dissimilarity between all ligand-receptor pairs in the two conditions (condition 1 in rows, condition 2 in columns).
#'   }
#' }
#'
"comparative_analysis"

# load AML328_d0_interactions
data("AML328_d0_interactions")

# load AML328_d31_interactions
data("AML328_d29_interactions")

AML_comparative_analysis <- comparative_analysis(cond1_weight_array = AML328_d0_interactions$weight_array
                                             ,cond2_weight_array = AML328_d29_interactions$weight_array
                                             ,cond1_ligand_receptor_pair_df = AML328_d0_interactions$ligand_receptor_pair_df
                                             ,cond2_ligand_receptor_pair_df = AML328_d29_interactions$ligand_receptor_pair_df
                                             ,cond1_nodes = AML328_d0_interactions$nodes
                                             ,cond2_nodes = AML328_d29_interactions$nodes
                                             ,cond1_name = "AML328_d0"
                                             ,cond2_name = "AML328_d29")

use_data(AML_comparative_analysis
         ,overwrite = TRUE)
