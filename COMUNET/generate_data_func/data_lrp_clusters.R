#' Clusters of ligand-receptor pairs.
#'
#' Dataset contains a matrix of pairwise dissimilarity between all ligand-receptor pairs,
#' cluster assignment for each ligand-receptor pair,
#' array of weight adjacency matrices (one weight adjacency matrix per cluster).
#'
#' @format A data list of three variables
#' \describe{
#' \itemize{
#'   \item{dissim_matrix}{
#'
#'   Numeric matrix: pairwise dissimilarity between all ligand-receptor pairs.
#'   }
#'   \item{clusters}{
#'
#'   Numeric vector: cluster assignment for each ligand-receptor pair.
#'   }
#'   \item{weight_array_by_cluster}{
#'
#'   Numeric array:   array of weight adjacency matrices (one weight adjacency matrix per cluster).
#'
#'   The array of weight adjacency has the following dimensions: number of nodes, number of nodes, number of ligand-receptor clusters.
#'
#'   Each edge weight in a cluster is calculated as the mean of the weights of this edge in all ligand-receptor pairs that belong to this cluster.
#'
#'   Each delta degree of a node in a cluster is calculated as the mean of the delta degrees of this node in all ligand-receptor pairs that belong to this cluster.
#'   }
#'   }
#' }
"lrp_clusters"

options(stringsAsFactors = F)

# load embryo_interactions data
data("embryo_interactions")

# calcualte clusters if ligand-receptor pairs
lrp_clusters <- lrp_clustering(weight_array = embryo_interactions$weight_array
                               ,ligand_receptor_pair_df = embryo_interactions$ligand_receptor_pair_df
                               ,nodes = embryo_interactions$nodes
)

# print(str(lrp_clusters))

use_data(lrp_clusters
         ,overwrite = TRUE)
