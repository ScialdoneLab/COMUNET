#' @rdname plot_communication_graph
#' @export
#'
#' @title
#' Plots communication graph
#'
#' @description
#' Plots a communication graph for a pair of interacting partners.
#'
#' If the pair of interacting partners is a ligand-receptor pair, a positive value of the node delta degree indicates that the node (cell type)
#' is mostly communicating by producing and secreting the ligand ("sending" node); conversely, a negative node delta degree marks nodes
#' that communicate mainly by receiving signals through the receptor ("receiving node").
#'
#' The arrows of the graph start at the sending nodes (which expresses the ligand) and point to the receiving nodes (which expresses the receptor),
#' while the thickness of an edge indicates the edge weight.
#' In case no directionality is specified for the pair of interacting proteins A and B (i.e., as for adhesion molecules),
#' the arrows start at the node expressing partner A and point to the node expressing partner B.
#'
#' Please nota that for simplicity, we address all interacting partners (including non-directional partners such as adhesion molecules) as "ligand-receptor pairs".
#'
#' @author
#' Maria Solovey \email{maria.solovey@helmholtz-muenchen.de}
#'
#' @param      LRP Character string: name of ligand-receptor pair for which the plot should be constructed.
#'
#' Note the the LRP name should be written in the same from as in ligand_receptor_pair_df$pair.
#' Alternatively, it can be a name of a cluster or a pattern.
#'
#' @param       weight_array Numeric array (3D) or numeric matrix (2D):
#' array of weighted adjacency matrices with dimensions [number of nodes, number of nodes, number of ligand-receptor pairs] or one adjacency matrix with dimensions [number of nodes, number of nodes].
#'
#' \itemize{
#'
#' \item first dimension: sending nodes (note: this dimension has all possible nodes, even if some of them are silent for a particular ligand-receptor pair).
#' \item second dimension: receiving nodes (note: this dimension has all possible nodes, even if some of them are silent for a particular ligand-receptor pair).
#' \item third dimension (for arrays): ligand-receptor pairs.
#'}
#'
#' Note that the weight_array should contain dimension names: dimnames = list(nodes, nodes, ligand-receptor pairs).
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
#' @param       nodes Character vector: a vector with all cell types in the data.
#'
#' @param       is_pattern Logical: should a graph of a pattern (and not a particular ligand-receptor pair) be plotted. Default value = FALSE.
#'
#' @param       is_cluster Logical: should a graph of a cluster (and not a particular ligand-receptor pair) be plotted. Default value = FALSE.
#'
#' @param       title Character string: title of the plot. Default value: LRP.
#'
#' @param       subtitle Character string:  subtitle of the plot. Default value: "".
#'
#' @param       node_color_palette Character string vector: vector of colours for nodes. Default values: c("blue", "gray40", "red").
#'
#' The colours will be used to construct a colour gradient for the node delta degree: blue representing a receiving node (negative delta degree),
#' grey representing a neutral node (zero delta degree), red representing a sending node (positive delta degree).
#'
#' @param       node_label_cex Numeric: size of node labels. Default value: 1.
#'
#' @param       vertex_shape Character string: shape of node. Default value: "none". See also \code{vertex.shape} argument of the \code{\link[plot.igraph:igraph]{plot.igraph}} function.
#'
#' @param       vertex_size Integer: size of node. Default value: 15.
#'
#' @param       edge_arrow_size  Numeric: size of arrow edge. Default value: 0.5. See also \code{edge.arrow.size} argument of the \code{\link[plot.igraph:igraph]{plot.igraph}} function.
#'
#' @param       edge_width_scale Numeric: scaling factor for edge width. Default value: 2.5.
#'
#' @param       legend_size Numeric: size of the legend text. Default value: 0.75.
#'
#' @param       legend_gradient_x_coord Numeric vector: x coordinates of the label position. Default value: c(1.15,1.25,1.25,1.15).
#'
#' @param       legend_gradient_y_coord Numeric numeric vector: y coordinates of the label position. Default value: c(-0.5,-0.5,-1,-1).
#'
#' @param       ... Any other parameters of the \code{\link[plot.igraph:plot.igraph]{plot.igraph}} function..
#'
#' @return      Graph plot.
#'
#' @examples
#'
#' # load embryo_interactions
#' data("embryo_interactions")
#'
#' # load lrp_clusters
#' data("embryo_lrp_clusters")
#'
#' # For a single ligand-receptor pair:
#' plot_communication_graph(LRP = "IGF2:IGF2R",
#'                         weight_array = embryo_interactions$weight_array,
#'                         ligand_receptor_pair_df = embryo_interactions$ligand_receptor_pair_df,
#'                         nodes = interactions$node)
#'
#' # For a cluster:
#' plot_communication_graph(LRP = "cluster 1",
#'                        weight_array = embryo_lrp_clusters$weight_array_by_cluster[,,"cluster 1"],
#'                        ligand_receptor_pair_df = embryo_interactions$ligand_receptor_pair_df,
#'                        nodes = embryo_interactions$nodes,
#'                        is_cluster = T)
#'
#' # for a patern:
#' test_pattern <- matrix(c(rep(0
#'                             ,dim(embryo_interactions$weight_array)[[1]]-1)
#'                             ,1)
#'                        ,nrow = dim(embryo_interactions$weight_array)[[1]]
#'                        ,ncol = dim(embryo_interactions$weight_array)[[2]]
#'                        ,byrow = F
#'                        )
#' dimnames(test_pattern) <- dimnames(embryo_interactions$weight_array)[c(1,2)]
#'
#' plot_communication_graph(LRP = "my pattern of interest",
#'                        weight_array = test_pattern,
#'                        ligand_receptor_pair_df = embryo_interactions$ligand_receptor_pair_df,
#'                        nodes = embryo_interactions$node,
#'                        is_pattern = T)
#'
plot_communication_graph <- function(LRP
                                     ,weight_array
                                     ,ligand_receptor_pair_df
                                     ,nodes
                                     ,is_pattern = FALSE
                                     ,is_cluster = FALSE
                                     ,title = LRP
                                     ,subtitle = ""
                                     ,node_color_palette = c("blue", "gray40", "red")
                                     ,node_label_cex = 1
                                     ,vertex_shape = "none"
                                     ,vertex_size = 15
                                     ,edge_arrow_size = 0.5
                                     ,edge_width_scale = 2.5
                                     ,legend_size = 0.75
                                     ,legend_gradient_x_coord = c(1.15,1.25,1.25,1.15)
                                     ,legend_gradient_y_coord = c(-0.5,-0.5,-1,-1)
                                     ,...
){
        # check that weight_array has dimnames
        if(length(dimnames(weight_array)) == 0) stop("ERROR: Dimension names of weight_array are not defined. Please define dimnames(weight_array)")

        # check whether the LRP is present in ligand_receptor_pair_df
        # In some cases of between sample comparison it might be absent,
        # so we will have to plot an empty plot
        present <- LRP %in% ligand_receptor_pair_df$pair

        # define weight matrix
        if(is_pattern | is_cluster) {
                adjacency_matrix <- weight_array
        } else if(present){ # LRP is present in the ligand_receptor_pair_df

                # subset array
                adjacency_matrix <- weight_array[,,LRP]
                # convert to matrix
                adjacency_matrix <- array2matrix(adjacency_matrix)

        } else { # LRP is NOT present in the ligand_receptor_pair_df

                # create an empty weight matrix
                adjacency_matrix <- create_empty_adj_matrix(nodes)
                # convert to matrix
                adjacency_matrix <- array2matrix(adjacency_matrix)
        }

        # calculate delta degree matrix
        degree_matrix <- calc_degrees(weight_array = adjacency_matrix)
        # convert to matrix
        degree_matrix <- array2matrix(degree_matrix)

        # create a directed weighted graph object
        layer_graph <- igraph::graph.adjacency(adjacency_matrix
                                               ,mode = "directed"
                                               ,weighted = T
        )

        # add colour scheme to the graph object
        colors <- delta_degree2color(degree_matrix[,3]
                                     ,nodes
                                     ,node_color_palette
        )
        V(layer_graph)$color <- colors

        # create colour map for the legend
        colormap <- colorRampPalette(c(min(colors)
                                       ,node_color_palette[2]
                                       ,max(colors)
        )
        )(20)
        leg_labels <- c(round(min(degree_matrix[,3]),2)
                        ,round(max(degree_matrix[,3]),1)
        )
        leg_title <- "delta\ndegree"

        # add size of node labels
        V(layer_graph)$label.cex = node_label_cex

        V(layer_graph)$size = vertex_size

        # define edge.width
        if(length(E(layer_graph)$weight) != 0){
                edge.width <- (E(layer_graph)$weight/max(E(layer_graph)$weight))*edge_width_scale
        }else edge.width <- 0

        # set seed
        set.seed(1123)

        # create a plot
        igraph::plot.igraph(
                layer_graph
                ,edge.width=edge.width
                ,edge.arrow.size = edge_arrow_size
                ,layout = layout.circle(layer_graph)
                ,edge.curved=TRUE
                ,vertex.label.color=V(layer_graph)$color
                ,vertex.shape=vertex_shape
                ,main = title
                ,sub = subtitle
                ,... # other plot.igraph parameters
        )

        #points for the gradient legend
        pnts = cbind(x = legend_gradient_x_coord
                     ,y = legend_gradient_y_coord
        )
        #create the gradient legend
        SDMTools::legend.gradient(
                pnts = pnts
                ,cols = colormap
                ,limits = leg_labels
                ,title = leg_title
                ,cex = legend_size
        )


}
