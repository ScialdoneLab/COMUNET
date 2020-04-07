#' @rdname plot_communication_graph
#' @export
#'
#' @title
#' Plots communication graph
#'
#' @description
#' \code{plot_communication_graph} plots a communication graph for a pair of interacting partners.
#' If the pair of interacting partners is a ligand-receptor pair, a positive value of the node delta degree indicates that the node (cell type)
#' is mostly communicating by producing and secreting the ligand (“sending” node); conversely, a negative node delta degree marks nodes
#' that communicate mainly by receiving signals through the receptor (“receiving node”).
#' The arrows of the graph start at the sending nodes (which expresses the ligand) and point to the receiving nodes (which expresses the receptor),
#' while the thickness of an edge indicates the edge weight.
#' In case no directionality is specified for the pair of interacting proteins A and B (i.e., as for adhesion molecules),
#' the arrows start at the node expressing partner A and point to the node expressing partner B.
#'
#' @author
#' Maria Solovey \email{maria.solovey@helmholtz-muenchen.de}
#'
#' @param      LRP character string: name of ligand-receptor pair for which the plot should be constructed.
#' Note the the LRP name should be written in the same from as in ligand_receptor_pair_df$pair.
#' Alternatively, it can be a name of a cluster or a pattern.
#'
#' @param       weight_array 3D numeric array or 2D numeric matrix:
#' array of weighted adjacency matrices with dimensions [number of nodes, number of nodes, number of ligand-receptor pairs] or one adjacency matrix with dimensions [number of nodes, number of nodes].
#'
#' First dimension: sending nodes (note: this dimension has all possible nodes, even if some of them are silent for a particular ligand-receptor pair).
#'
#' Second dimension: receiving nodes (note: this dimension has all possible nodes, even if some of them are silent for a particular ligand-receptor pair).
#'
#' Third dimension (for arrays): ligand-receptor pairs.
#'
#' Note that the weight_array should contain dimnames: dimnames = list(nodes, nodes, ligand-receptor pairs).
#'
#' @param       ligand_receptor_pair_df character dataframe: data frame with columns "pair", "ligand", "ligand_complex_composition", "receptor", "receptor_complex_composition".
#'
#' Column "pair" contains values in a form "ligand:receptor", i.e. ligand being at the first place, receptor being at the second place, e.g. "TNFSF13:TNFRSF17".
#'
#' Column "ligand" contains ligand names, e.g. "TNFSF13".
#'
#' Column ligand_complex_composition" if ligand is a complex (e.g. "aXb2_complex"),
#' contains genes in the ligand complex separated with a comma, e.g. "ITGAX,ITGB2", else contains empty string "".
#'
#' Column "receptor" contains receptor names, e.g. "TNFRSF17".
#'
#' Column "receptor_complex_composition" if receptor is a complex (e.g. "NKG2D_II_receptor"),
#' contains genes in the receptor complex separated with a comma, e.g. "KLRK1,HCST", else contains empty string "".
#'
#' @param       nodes character vector: a vector with all cell types in the data.
#'
#' @param       is_pattern logical: should a graph of a pattern (and not a particular ligand-receptor pair) be plotted. Default value = F.
#'
#' @param       is_cluster logical: should a graph of a cluster (and not a particular ligand-receptor pair) be plotted. Default value = F.
#'
#' @param       title character string: title of the plot. Default value: LRP.
#'
#' @param       subtitle character string:  subtitle of the plot. Default value: "".
#'
#' @param       node_color_palette character string vector: vector of colours for nodes. Default values: c("blue", "gray40", "red").
#' The colours will be used to construct a colour gradient for the node delta degree: blue representing a receiving node (negative delta degree),
#' grey representing a neutral node (zero delta degree), red representing a sending node (positive delta degree).
#'
#' @param       node_label_cex numeric: size of node labels. Default value: 1.
#'
#' @param       vertex_shape character string: shape of node. Default value: "none". See same argument of plot.igraph().
#'
#' @param       vertex_size integer: size of node. Default value: 15. See same argument of plot.igraph().
#'
#' @param       edge_arrow_size  numeric: size of arrow edge. Default value: 0.5. See same argument of plot.igraph().
#'
#' @param       edge_width_scale numeric: scaling factor for edge width. Default value: 2.5. See same argument of plot.igraph().
#'
#' @param       legend_size numeric: size of the legend text. Default value: 0.75.
#'
#' @param       legend_gradient_x_coord numeric vector: x coordinates of the label position. Default value: c(1.15,1.25,1.25,1.15).
#'
#' @param       legend_gradient_y_coord numeric numeric vector: y coordinates of the label position. Default value: c(-0.5,-0.5,-1,-1).
#'
#' @param       ...                             other plot.igraph parameters
#'
#' @return      graph plot
plot_communication_graph <- function(LRP
                                     ,weight_array
                                     ,ligand_receptor_pair_df
                                     ,nodes
                                     ,is_pattern = F
                                     ,is_cluster = F
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
