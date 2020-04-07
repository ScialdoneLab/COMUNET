#' @rdname plot_lig_rec
#' @export
#'
#' @title
#' Plots ligands and receptors in a particular cluster
#'
#' @description
#' \code{plot_lig_rec} plots interacting partners in a particular cluster.
#' If the pair of interacting partners is a ligand-receptor pair, then the ligand is colored red,
#' the receptor is colored blue, and the arrow goes from the ligand to the receptor.
#' In case no directionality is specified for the pair of interacting proteins A and B (i.e., as for adhesion molecules),
#' the arrows start at the node expressing partner A and point to the node expressing partner B.
#'
#' @author
#' Maria Solovey \email{maria.solovey@helmholtz-muenchen.de}
#'
#' @param       cluster_of_interest Integer: number of the cluster, for which the ligands and receptors should be plotted.
#'
#' @param       lrp_clusters Integer vector: cluster assignment for each ligand-receptor pair.
#'
#' @param       ligand_receptor_pair_df Character dataframe: data frame with columns "pair", "ligand", "ligand_complex_composition", "receptor", "receptor_complex_composition".
#'
#' Column "pair" contains values in a form "ligand:receptor", i.e. ligand being at the first place, receptor being at the second place, e.g. "TNFSF13:TNFRSF17".
#'
#' Column "ligand" contains ligand names, e.g. "TNFSF13".
#'
#' Column "ligand_complex_composition" if ligand is a complex (e.g. "aXb2_complex"),
#' contains genes in the ligand complex separated with a comma, e.g. "ITGAX,ITGB2", else contains empty string "".
#'
#' Column "receptor" contains receptor names, e.g. "TNFRSF17".
#'
#' Column "receptor_complex_composition" if receptor is a complex (e.g. "NKG2D_II_receptor"),
#' contains genes in the receptor complex separated with a comma, e.g. "KLRK1,HCST", else contains empty string "".
#'
#' @param       lig_rec_color Character string vector: colours for ligand and receptors: first colour for the ligand, second colour for the receptor. Default value: c("red", "blue").
#'
#' @param       edge_arrow_size Numeric: size of arrow edge. Default value: 0.25. See same argument of plot.igraph().
#'
#' @param       node_label_cex Numeric: size of node labels. Default value: 1.
#'
#' @param       layout Function: see layout from igraph package v0.2.1. Default value: layout_nicely.
#'
#' @param       legend_position Numeric vector: x, y coordinates of the legend position. Default value: c(-1, -1.1).
#'
#' @param       legend_pt_cex Numeric: expansion factor(s) for the points in the legend. See legend() from graphics v3.6.2. Default value: 2.
#'
#' @param       legend_cex Numeric: character expansion factor in the legend. See legend() from graphics v3.6.2. Default value: 0.8.
#'
#' @param       ... Any other plot.igraph parameters.
#'
#' @return     graph plot
plot_lig_rec <- function(cluster_of_interest
                         ,lrp_clusters
                         ,ligand_receptor_pair_df
                         ,lig_rec_color = c("red", "blue")
                         ,edge_arrow_size = 0.25
                         ,node_label_cex = 1
                         ,layout = layout_nicely
                         ,legend_position = c(-1, -1.1)
                         ,legend_pt_cex=2
                         ,legend_cex=0.8
                         ,...
){
        # make adjacency matrix:
        # rows contain ligands
        # columns contain receptors
        # 1 if ligand and receptor is a pair, 0 if ligand and receptor is no pair
        lrp <- ligand_receptor_pair_df[lrp_clusters == cluster_of_interest,]
        ligands <- unique(lrp$ligand)
        receptors <- unique(lrp$receptor)
        adj_matrix <- matrix(0
                             ,nrow = length(ligands) + length(receptors)
                             ,ncol = length(ligands) + length(receptors)
        )
        dimnames(adj_matrix) <- list(c(ligands,receptors)
                                     ,c(ligands,receptors)
        )
        for(idx in 1:nrow(lrp)){
                idx_row <- which(lrp$ligand[idx] == rownames(adj_matrix))
                idx_col <- which(lrp$receptor[idx] == colnames(adj_matrix))
                adj_matrix[idx_row, idx_col] <- 1
        }

        # create a directed weighted graph object
        lrp_graph <- igraph::graph.adjacency(adj_matrix
                                             ,mode = "directed"
        )
        # add color
        V(lrp_graph)$color <- {
                color <- c(rep(lig_rec_color[1]
                               ,length(ligands))
                           ,rep(lig_rec_color[2]
                                ,length(receptors)))
                names(color) <- c(ligands
                                  ,receptors)
                color <- color[colnames(adj_matrix)]
                color
        }
        # add label size
        V(lrp_graph)$label.cex = node_label_cex

        legend_lab <- c("ligand"
                        ,"receptor")
        legend_col <- lig_rec_color
        title <- paste("cluster"
                       ,cluster_of_interest)

        # set seed
        set.seed(1123)

        # create a plot
        igraph::plot.igraph(lrp_graph
                            ,edge.arrow.size = edge_arrow_size
                            ,vertex.label.color=V(lrp_graph)$color
                            ,vertex.shape="none"
                            ,edge.curved=F
                            ,main = title
                            ,layout=layout
                            ,... # other plot.igraph parameters
        )
        legend(x=legend_position[1]
               ,y=legend_position[2]
               ,legend = legend_lab
               ,pch=21
               ,pt.bg=legend_col
               ,pt.cex=legend_pt_cex
               ,cex=legend_cex
               ,bty="n"
               ,ncol=1
        )
}
