#' @rdname lrp_clustering
#' @export
#'
#' @title
#' Clusters interacting protein pairs
#'
#' @description
#' Clusters interacting protein pairs by performing hierarchical clustering of the dissimilarity matrix of layers
#' (by default, the \code{\link[hclust:stats]{hclust}} R function with the \code{average} agglomeration method is used) and the results can be visualized as a heatmap or a UMAP plot.
#'
#' The number of clusters is estimated using the  \code{\link[cutreeHybrid:dynamicTreeCut]{cutreeHybrid}} R function
#' with \code{deepSplit} equal to 0 and default \code{minClusterSize} equal to 6 by default (both can be adjusted).
#'
#' For each cluster, a graph represents the average communication pattern in the cluster and is built by averaging the edge weights and the delta degrees of the nodes of all the graphs in the cluster.
#'
#' Please nota that for simplicity, we address all interacting partners (including non-directional partners such as adhesion molecules) as "ligand-receptor pairs".
#'
#' @author
#' Maria Solovey \email{maria.solovey@helmholtz-muenchen.de}
#'
#' @param       weight_array Numeric array: array of weighted adjacency matrices with dimensions [number of nodes, number of nodes, number of ligand-receptor pairs].
#'
#' @param       ligand_receptor_pair_df Character string data frame: data frame with columns:
#'
#'  \itemize{
#'  \item "pair" contains values in a form "ligand:receptor", i.e. ligand being at the first place, receptor being at the second place, e.g. "TNFSF13:TNFRSF17".
#'  \item "ligand" contains ligand names, e.g. "TNFSF13".
#'  \item "ligand_complex_composition" if ligand is a complex (e.g. "aXb2_complex"),
#'  contains genes in the ligand complex separated with a comma, e.g. "ITGAX,ITGB2", else contains empty string "".
#'  \item "receptor" contains receptor names, e.g. "TNFRSF17".
#'  \item "receptor_complex_composition" if receptor is a complex (e.g. "NKG2D_II_receptor"),
#'  contains genes in the receptor complex separated with a comma, e.g. "KLRK1,HCST", else contains empty string "".
#'  }
#'
#' @param       nodes Character string vector: a vector with all cell types in the data.
#'
#' @param       dissimilarity Function: dissimilarity function. Default value is d_normWeightDiff.
#'
#' @param       deep_split Numeric: defines the \code{deepSplit} parameter of the \code{\link[cutreeHybrid:dynamicTreeCut]{cutreeHybrid}} function. Default value: 0.
#'
#' @param       min_cluster_size Numeric: defines the \code{minClusterSize} parameter of the \code{\link[cutreeHybrid:dynamicTreeCut]{cutreeHybrid}} function, i.e. the minimum number of ligand-receptor pairs per cluster. Default value: 6.
#'
#' @return A list of:
#'
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
#'
#' @examples
#' # load embryo_interactions data
#' data("embryo_interactions")
#'
#' # calculate clusters if ligand-receptor pairs
#' lrp_clusters <- lrp_clustering(weight_array = embryo_interactions$weight_array,
#'           ligand_receptor_pair_df = embryo_interactions$ligand_receptor_pair_df,
#'           nodes = embryo_interactions$nodes
#'           )
#'
#' # print(str(lrp_clusters))

lrp_clustering <- function(weight_array
                           ,ligand_receptor_pair_df
                           ,nodes
                           ,dissimilarity = d_normWeightDiff
                           ,deep_split = 0
                           ,min_cluster_size = 6

){
        # calculate degree array
        degree_array <- calc_degrees(weight_array = weight_array)
        hist(degree_array[,3,], main = "Histogram of degree delta")

        ## run analysis ####
        # check if there are more than 1 ligand receptor pairs
        if(nrow(ligand_receptor_pair_df) <= 1) {
                dissim_matrix <- NA
                my_clusters <- 0
                weight_array_by_cluster <- NA
                degree_array_for_clusters <- NA
        } else {
                # dissimilatriy ####
                dissim_matrix <- calc_dissim_matrix(weight_array1 = weight_array
                                                    ,weight_array2 = weight_array
                                                    ,dissimilarity = dissimilarity
                )
                hist(dissim_matrix)

                # cluster layers ####
                # convert dissimilarity matrix to distance matrix
                my_dist <- as.dist(dissim_matrix)
                # Make hierarchical clustering
                my_tree <- hclust(my_dist
                                  , method= "average"
                )

                # Cut the tree to identify clusters
                my_clusters <- dynamicTreeCut::cutreeHybrid(
                        my_tree
                        ,distM=dissim_matrix
                        ,deepSplit = deep_split
                        ,minClusterSize = min_cluster_size
                )$label
                names(my_clusters) <- dimnames(dissim_matrix)[[1]]

                # Any unassigned cells?
                if(sum(my_clusters == 0) != 0) print("Warning: some graphs are not assigned to any cluster")
                # print out number of clusters
                print(paste("We have", max(as.numeric(my_clusters)), "clusters"))

                # calculate adjacency matrices for each community (using mean weight for each edge)
                if(any(my_clusters != 0)){
                        weight_array_by_cluster <- array(NA
                                                         ,dim = c(length(nodes)
                                                                  ,length(nodes)
                                                                  ,max(my_clusters)
                                                         )
                        )
                        for(communityNr in 1:max(my_clusters)){
                                # calculate indices for the community
                                myCommunityIdx <- which(as.numeric(my_clusters) == communityNr)
                                # populate matrix
                                for(i in 1:length(nodes)){
                                        for(j in 1:length(nodes)){
                                                weight_array_by_cluster[i,j,communityNr] <- mean(weight_array[i,j,myCommunityIdx])
                                        }
                                }
                        }
                        dimnames(weight_array_by_cluster) <- list(nodes
                                                                  ,nodes
                                                                  ,paste("cluster"
                                                                         ,1:max(my_clusters))
                        )
                        check_NAs(weight_array_by_cluster)

                        # calculate degree arrays for each community (put it as well as an array)
                        degree_array_for_clusters <- calc_degrees(weight_array = weight_array_by_cluster)

                } else {
                        weight_array_by_cluster <- NA
                        degree_array_for_clusters <- NA
                }

                # end of definitions and data processing

        }


        ## produce report ####
        # make object ####
        result <- list(dissim_matrix = dissim_matrix
                       ,clusters = my_clusters
                       ,weight_array_by_cluster = weight_array_by_cluster
        )


        ## return result ####
        return(result)
}
