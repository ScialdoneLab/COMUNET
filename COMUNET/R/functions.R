# Functions ####

# make a function to check NAs
# INPUT: 
#       object_to_check: a vector or matrix containing any values
# OUPUT:
#       a warning message is case the object contains NAs
check_NAs <- function(object_to_check){
        if(!all(!is.na(object_to_check))) stop("WARNING: the object contains NAs.")
}

step_func <- function(x) {as.integer(x != 0)}

# prevent dropping dimensions
`[` <- function(...) base::`[`(...,drop=FALSE)

# filter out layers with no edges
filter_empty_layers <- function(weight_array
                                ,ligand_receptor_pair_df
){
        # check which layers have no edges
        no_edge_layers <- sapply(1:nrow(ligand_receptor_pair_df)
                                 ,function(layer){
                                         sum(weight_array[,,layer]) == 0
                                 })
        names(no_edge_layers) <- dimnames(weight_array)[3]
        # delete these layers form ligand_receptor_pair_df
        ligand_receptor_pair_df <- ligand_receptor_pair_df[!no_edge_layers,]
        # delete these layers from the array of the adjacency matrices
        weight_array <- weight_array[,,!no_edge_layers]
        # return
        return(list(ligand_receptor_pair_df = ligand_receptor_pair_df
                    ,weight_array = weight_array
        )
        )
}

# function to convert an array of dimensions [n x n x 1] to an [n x n] matrix
array2matrix <- function(my_2D_array){
        
        my_matrix <- matrix(my_2D_array
                            ,nrow = dim(my_2D_array)[1]
                            ,ncol = dim(my_2D_array)[2]
        )
        dimnames(my_matrix) <- dimnames(my_2D_array)[c(1,2)]
        my_matrix
}

# calculate an array of three dimensions: (nr of nodes) x (nr of degrees) x (nr of layers).
# nr of nodes is N
# nr of degrees in 3: in- out- and (out-in)-
# nr of layers is M
#function to calculate degree
calc_degrees <- function(weight_array
){
        N <- dim(weight_array)[1]
        if(length(dim(weight_array)) == 3){
                M <- dim(weight_array)[3]
        } else { M <- 1}
        
        # make a dummy array
        my_degree_array <- array(NA
                                 ,dim = c(N,3,M)
        )
        
        # add dimension names
        if(M == 1){
                dimnames(my_degree_array) <- list(dimnames(weight_array)[[1]]
                                                  ,c("in", "out", "delta")
                                                  ,1
                )        
        } else {
                dimnames(my_degree_array) <- list(dimnames(weight_array)[[1]]
                                                  ,c("in", "out", "delta")
                                                  ,dimnames(weight_array)[[3]]
                )
        }
        
        my_degree_array[,"in",] <- colSums(weight_array)
        if(M==1){
                my_degree_array[,"out",] <- rowSums(weight_array)
        }else{
                my_degree_array[,"out",] <- colSums(aperm(weight_array
                                                          ,c(2,1,3)))
        }
        my_degree_array[,"delta",] <- my_degree_array[,"out",] - my_degree_array[,"in",]
        # check NAs
        check_NAs(my_degree_array)
        # return my_degree_array
        my_degree_array
}

# dissimilarity function
# IN
#       adj1: numeric matrix [0,1], adjacency matrix for the layer 1
#       adj2: numeric matrix [0,1], adjacency matrix for the layer 2
# OUT
#       dissimCoefficient: numeric [0,1], the dissimilarity coefficient for the two layers. The identical matrices will have a coefficient of 0 and completely different matrices will have a coefficient of 1
d_normWeightDiff <- function(adj1, adj2){
        # make a union adjacency matrix: 
        # 1. apply theta function on each of the two matrices to convert weights into zeros and ones
        # 2. then perform bitwise or operation to get the union matrix
        union_of_edges <- apply(adj1, 1:2, step_func) | apply(adj2, 1:2, step_func)
        # calculate the number of edges in the union
        n_union_of_edges <- sum(union_of_edges)
        
        if(n_union_of_edges == 0) my_dissim <- 1 else {
                # calculate the sum of normalised weight differences
                # matrix of absolute differences of weights for each edge
                abs_weight_difference <- abs(adj1 - adj2)
                # matrix of sum of weights for each edge
                wight_sum <- adj1 + adj2
                # Convert zeros to -1 to avoid division by zero in the next step
                wight_sum[wight_sum == 0] <- -1
                # matrix of normalised differences 
                norm_weight_difference <- abs_weight_difference / wight_sum
                # 
                norm_weight_difference_sum <- sum(norm_weight_difference)
                my_dissim <- norm_weight_difference_sum / n_union_of_edges
        }
        my_dissim
}

# calculate dissimilarity
calc_dissim_matrix <- function(weight_array1
                               ,weight_array2
                               ,dissimilarity = d_normWeightDiff
){
        d <- dissimilarity
        
        N <- dim(weight_array1)[[3]]
        M <- dim(weight_array2)[[3]]
        # make a dummy matrix
        dissim_matrix <- matrix(NA
                                ,nrow = N
                                ,ncol = M
        )
        # calculate an N x M dissimilarity matrix for all layers
        if(identical(weight_array1,weight_array2)) { # if it is the same weight matrix, calculate only the upper triangular matrix
                # populate diagonal with 0
                diag(dissim_matrix) <- 0
                
                for (i in 1:(N-1)){
                        for(j in (i+1):M){
                                dissim_value <- d(weight_array1[,,i]
                                                  ,weight_array2[,,j]
                                )
                                dissim_matrix[i,j] <- dissim_value
                                dissim_matrix[j,i] <- dissim_value
                        }
                }
        } else {
                for (i in 1:N){
                        for(j in 1:M){
                                dissim_value <- d(weight_array1[,,i]
                                                  ,weight_array2[,,j]
                                )
                                dissim_matrix[i,j] <- dissim_value
                        }
                }
        }
        
        check_NAs(dissim_matrix)
        dimnames(dissim_matrix) <- list(dimnames(weight_array1)[[3]]
                                        ,dimnames(weight_array2)[[3]]
        )
        dissim_matrix
}

# assigns colour to each node according to its delta degree
# IN 
# delta_degree_vector
# nodes
# palette
# OUT
delta_degree2color <- function(delta_degree_vector
                               ,nodes
                               ,node_color_palette
){
        fine <- round(max(abs(delta_degree_vector)) 
                      ,2) * 100 * 2 + 2
        if(fine == 2) {
                colors <- rep("gray40"#FFFFFF"
                              ,length(delta_degree_vector))
                names(colors) <- nodes
                colors
        } else {
                palette <- colorRampPalette(node_color_palette)
                sapply(delta_degree_vector
                       ,function(i){
                               palette(fine)[ceiling(round(i,2)*100 
                                                     + ((fine - 1)/2)
                               )
                               ]
                       })       
        }
}

create_empty_adj_matrix <- function(nodes){
        adjacency_matrix <- matrix(0
                                   ,nrow=length(nodes)
                                   ,ncol=length(nodes)
        )
        dimnames(adjacency_matrix) <- list(nodes
                                           ,nodes)
        adjacency_matrix
}


# in case of a complex ligand / receptor, converts the complex name into a vector of individual ligands / receptors that comprise the complex
complex2genes <- function(complex
                          ,complex_input
                          ,gene_input
){
        if(complex %in% complex_input$complex_name){
                uniprots <- complex_input[complex_input$complex_name == complex
                                          ,2:5]
                # remove NA and ""
                uniprots <- uniprots[(!is.na(uniprots)) & (!(uniprots == "")) ]
                # translate uniprots to given gene names
                genes <- sapply(uniprots
                                ,function(uniprot) {
                                        unlist(gene_input[gene_input$uniprot == uniprot
                                                          ,"gene_name"])
                                }
                )
                genes
        } else {
                gene <- unique(unlist(gene_input[gene_input$gene_name == complex
                                                 ,"gene_name"]))
                gene
        }
}

#' @rdname plot_communication_graph
#' @export
#' 
#' @title 
#' Plots communication graph
#' 
#' @description 
#' \code{plot_lig_rec} plots a communication graph for a pair of interacting partners.
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
#' First dimension: sending nodes (note: this dimension has all possible nodes, even if some of them are silent for a particular ligand-receptor pair). 
#' Second dimension: receiving nodes (note: this dimension has all possible nodes, even if some of them are silent for a particular ligand-receptor pair). 
#' Third dimension (for arrays): ligand-receptor pairs. 
#' Note that the weight_array should contain dimnames: dimnames = list(nodes, nodes, ligand-receptor pairs).
#' 
#' @param       ligand_receptor_pair_df character dataframe: data frame with columns "pair", "ligand", "ligand_complex_composition", "receptor", "receptor_complex_composition". 
#' Column "pair" contains values in a form "ligand:receptor", i.e. ligand being at the first place, receptor being at the second place, e.g. "TNFSF13:TNFRSF17". 
#' Column "ligand" contains ligand names, e.g. "TNFSF13".  Column ligand_complex_composition" if ligand is a complex (e.g. "aXb2_complex"), 
#' contains genes in the ligand complex separated with a comma, e.g. "ITGAX,ITGB2", else contains empty string "". 
#' Column "receptor" contains receptor names, e.g. "TNFRSF17". 
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
#' @param       cluster_of_interest numeric: number of the cluster, for which the ligands and receptors should be plotted.
#' 
#' @param       lrp_clusters numeric vector: cluster assignment for each ligand-receptor pair.
#' 
#' @param       ligand_receptor_pair_df character dataframe: data frame with columns "pair", "ligand", "ligand_complex_composition", "receptor", "receptor_complex_composition". 
#' Column "pair" contains values in a form "ligand:receptor", i.e. ligand being at the first place, receptor being at the second place, e.g. "TNFSF13:TNFRSF17". 
#' Column "ligand" contains ligand names, e.g. "TNFSF13". 
#' Column "ligand_complex_composition" if ligand is a complex (e.g. "aXb2_complex"), 
#' contains genes in the ligand complex separated with a comma, e.g. "ITGAX,ITGB2", else contains empty string "". 
#' Column "receptor" contains receptor names, e.g. "TNFRSF17". 
#' Column "receptor_complex_composition" if receptor is a complex (e.g. "NKG2D_II_receptor"), 
#' contains genes in the receptor complex separated with a comma, e.g. "KLRK1,HCST", else contains empty string "".
#' 
#' @param       lig_rec_color character string vector: colours for ligand and receptors: first colour for the ligand, second colour for the receptor. Default value: c("red", "blue").
#' 
#' @param       edge_arrow_size numeric: size of arrow edge. Default value: 0.25. See same argument of plot.igraph().
#' 
#' @param       node_label_cex numeric: size of node labels. Default value: 1.
#' 
#' @param       layout function: see layout from igraph package v0.2.1. Default value: layout_nicely.
#' 
#' @param       legend_position numeric vector: x, y coordinates of the legend position. Default value: c(-1, -1.1).
#' 
#' @param       legend_pt_cex numeric: expansion factor(s) for the points in the legend. See legend() from graphics v3.6.2. Default value: 2.
#' 
#' @param       legend_cex numeric: character expansion factor in the legend. See legend() from graphics v3.6.2. Default value: 0.8.
#' 
#' @param       ... other plot.igraph parameters
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

#' @rdname plot_dissimilarity_heatmaps
#' @export
#' 
#' @title 
#' Plots heatmaps of dissimilarity between two conditions
#' 
#' @description 
#' \code{plot_dissimilarity_heatmaps} plots heatmaps of dissimilarity between two conditions. 
#' In each heatmap, interacting protein pairs from condition 1 are in the rows, and those from condition 2 are in the columns.
#' On the first heatmap, both rows and columns are clustered by the pairwise dissimilarity.
#' On the second heatmap, rows and columns are sorted by presence (shared or unshared) of an interacting protein pairs in both conditions
#' and by decreasing pairwise dissimilarity (in shared interacting pairs).
#' Labels of interacting pairs are colored by their presence (shared or unshared).
#' Each point of the heatmap is colored by the value of the pairwise dissimilarity.
#' 
#' @author
#' Maria Solovey \email{maria.solovey@helmholtz-muenchen.de}
#'
#' @param dissim_cond1_cond2 numeric matrix: pairwise dissimilarity between all ligand-receptor pairs in the two conditions (condition 1 in rows, condition 2 in columns).
#' 
#' @param sorted_LRP_df dataframe with columns: 
#' "pair" (character string): names of ligand-receptor pairs in the same form as they are in ligand_receptor_pair_df$pair; 
#' "presence" (character string): name of the condition in which the ligand-receptor pair is present or "shared" if present in both conditions; 
#' "dissimilarity" (numeric): dissimilarity value between the topology of ligand-receptor pair graph in two conditions. 
#' The smaller the dissimilarity value, the more similar is the graph topology between the two conditions. 
#' If a ligand-receptor pair is present only in one of the conditions, the dissimilarity is equal to 1.
#' 
#' @param cond1_name character string: sample name for condition 1.
#' 
#' @param cond2_name character string: sample name for condition 2.
#' 
#' @param colors_lrp character string vector: colours for ligand-receptor labels. Default: green for shared, black for unshared.
#' 
#' @param show_legend logical: parameter of the HeatmapAnnotation function. Default is T.
#' 
#' @param row_names_fontsize numeric: parameter of the Heatmap function. Default is 5.
#' 
#' @param colomn_names_fontsize numeric: parameter of the Heatmap function. Default is 5.
#' 
#' @param show_column_names logical: parameter of the Heatmap function. Default is T.
#' 
#' @param show_row_names logical: parameter of the Heatmap function. Default  is T.
#' 
#' @param width object of class "unit": parameter of the HeatmapAnnotation function. Default is unit(0.1, "cm").
#' 
#' @param legend_direction character string: parameter of the Heatmap function. Default is "horizontal".
#' 
#' @param legend_width object of class "unit": parameter of the Heatmap function. Default is unit(5, "cm").
#' 
#' @param title_position character string: parameter of the Heatmap function. Default is "lefttop".
#' 
#' @param row_dend_side character string: parameter of the Heatmap function. Default is "left"
#' 
#' @param column_dend_side character string: parameter of the Heatmap function. Default is "top".
#' 
#' @param column_title_side character string: parameter of the draw function. Default is "top".
#' 
#' @param heatmap_legend_side character string: parameter of the draw function. Default is "bottom".
#' 
#' @param ... other parameters of the Heatmap function.
#' 
#' @return     
#' two heatmaps:
#' one with ligand-receptor pairs clustered by their pairwise dissimilarity
#' one with ligand-receptor pairs sorted by decreasing dissimilarity of shared ligand-receptor pairs
plot_dissimilarity_heatmaps <- function(dissim_cond1_cond2
                                        ,sorted_LRP_df
                                        ,cond1_name
                                        ,cond2_name
                                        ,colors_lrp = c("green"
                                                        ,"black")
                                        ,show_legend = T
                                        ,row_names_fontsize = 5
                                        ,colomn_names_fontsize = 5
                                        ,show_column_names = T
                                        ,show_row_names = T
                                        ,width = unit(0.1, "cm")
                                        ,legend_direction = "horizontal"
                                        ,legend_width = unit(5, "cm")
                                        ,title_position = "lefttop"
                                        ,row_dend_side = "left"
                                        ,column_dend_side = "top"
                                        ,column_title_side = "top"
                                        ,heatmap_legend_side = "bottom"
                                        ,...
){
        # colour vector for ligand-receptor pairs in condition 1
        cond1_LRP_colors <- sapply(rownames(dissim_cond1_cond2)
                                   ,function(i){
                                           if(as.character(sorted_LRP_df[i,]$presence) == "shared") colors_lrp[1]
                                           else colors_lrp[2]
                                   })
        # colour vector for ligand-receptor pairs in condition 2
        cond2_LRP_colors <- sapply(colnames(dissim_cond1_cond2)
                                   ,function(i){
                                           if(as.character(sorted_LRP_df[i,]$presence) == "shared") colors_lrp[1]
                                           else colors_lrp[2]
                                   })
        # row annotations
        haRow = ComplexHeatmap::HeatmapAnnotation(
                df = data.frame(ligand_receptor_pair = c("shared"
                                                         ,"not shared")
                )
                ,col = list(ligand_receptor_pair = c("shared" = colors_lrp[1]
                                                     ,"not shared" = colors_lrp[2]))
                ,which = "row"
                ,width = width
                ,show_legend = show_legend
        )
        h_clustered <- ComplexHeatmap::Heatmap(
                dissim_cond1_cond2
                ,cluster_rows = T
                ,row_names_gp = gpar(fontsize = row_names_fontsize
                                     ,col = cond1_LRP_colors
                )
                ,cluster_columns = T
                ,column_names_gp = gpar(fontsize = colomn_names_fontsize
                                        ,col = cond2_LRP_colors
                )
                ,name = "dissimilarity coefficient"
                ,show_column_names = show_column_names
                ,show_row_names = show_row_names
                ,heatmap_legend_param = list(legend_direction = legend_direction
                                             ,legend_width = legend_width
                                             ,title_position = title_position
                )
                ,row_title = cond1_name
                ,end_side = row_dend_side
                ,column_title = cond2_name
                ,column_dend_side = column_dend_side
                ,...
        )
        heatmap_comparative_clustered <- ComplexHeatmap::draw(
                h_clustered + haRow
                ,column_title = "Comparative analysis: clustered"
                ,column_title_side = column_title_side
                ,heatmap_legend_side = heatmap_legend_side
        )
        
        h_sorted <- ComplexHeatmap::Heatmap(
                dissim_cond1_cond2
                ,cluster_rows = F
                ,row_order = {
                        sub_sortedList <- sorted_LRP_df[rownames(sorted_LRP_df) %in% rownames(dissim_cond1_cond2),]
                        rownames(sub_sortedList[order(sub_sortedList$presence),])
                }
                ,row_names_gp = gpar(fontsize = 5
                                     ,col = cond1_LRP_colors
                )
                ,cluster_columns = F
                ,column_order = {
                        sub_sortedList <- sorted_LRP_df[rownames(sorted_LRP_df) %in% colnames(dissim_cond1_cond2),]
                        rownames(sub_sortedList[order(sub_sortedList$presence),])
                }
                ,column_names_gp = gpar(fontsize = 5
                                        ,col = cond2_LRP_colors
                )
                ,name = "dissimilarity coefficient"
                ,show_column_names = show_column_names
                ,show_row_names = show_row_names
                ,heatmap_legend_param = list(legend_direction = legend_direction
                                             ,legend_width = legend_width
                                             , title_position = title_position
                )
                ,row_title = cond1_name
                ,row_dend_side = row_dend_side
                ,column_title = cond2_name
                ,column_dend_side = column_dend_side
                ,...
        )
        heatmap_comparative_sorted <- ComplexHeatmap::draw(
                h_sorted + haRow
                ,column_title = "Comparative analysis: sorted"
                ,column_title_side = column_title_side
                ,heatmap_legend_side = heatmap_legend_side
        )
}

#' @rdname plot_cluster_heatmap
#' @export
#' 
#' @title 
#' Plots heatmap for clustered ligand-receptor pairs
#' 
#' @description 
#' \code{plot_cluster_heatmap} plots a heatmap for clustered ligand-receptor pairs.
#' Rows and columns are sorted by cluster number, the color of the heatmap represents dissimilarity between two pairs.
#' 
#' @author
#' Maria Solovey \email{maria.solovey@helmholtz-muenchen.de}
#'
#' @param       dissim_matrix numeric matrix: pairwise dissimilarity matrix. Note that rownames and colnames should be defined.
#' 
#' @param       lrp_clusters numeric vector: cluster assignment for each ligand-receptor pair. 
#' Note that the the vector should be named and ordered in the same way as the rownames and colnames of the dissim_matrix.
#' 
#' @param       cluster_colors character string vector: vector with colours for each cluster. Default value: "".
#' 
#' @param       unclustered_LRP_color character string: colour for ligand-receptor pairs that remained unclustered (cluster assignment 0). Default value: "black".
#' 
#' @param       show_colomn_names logical: argument of the Heatmap function. Default value: F.
#' 
#' @param       show_row_names logical: argument of the Heatmap function. Default value: F.
#' 
#' @param       legend_direction character string: argument of the Heatmap function. Default value: "horizontal".
#' 
#' @param       legend_width object of class "unit": argument of the Heatmap function. Default is: unit(5, "cm").
#' 
#' @param       title_position character string: argument of the Heatmap function. Default value: "lefttop".
#' 
#' @param       column_title character string: title of the heatmap plot. Argument of the draw function. Default value: "Clustering of ligand-receptor pairs".
#' 
#' @param       column_title_side character string:  position of the title of the heatmap plot. Argument of the draw function. Default value: "top".
#' 
#' @param       heatmap_legend_side character string: position of the legend of the heatmap plot. Argument of the draw function. Default value: "bottom".
#' 
#' @param       ... any arguments "Heatmap" function takes.
#' 
#' @return      heatmap with ligand-receptor pairs ordered by their cluster assignment
plot_cluster_heatmap <- function(dissim_matrix
                                 ,lrp_clusters
                                 ,cluster_colors = ""
                                 ,unclustered_LRP_color = "black"
                                 ,show_column_names = F
                                 ,show_row_names = F
                                 ,legend_direction = "horizontal"
                                 ,legend_width = unit(5, "cm")
                                 ,title_position = "lefttop"
                                 ,column_title = "Clustering of ligand-receptor pairs"
                                 ,column_title_side = "top"
                                 ,heatmap_legend_side = "bottom"
                                 ,...
){
        # check that dissim_matrix  has rownames and colnames
        if(length(rownames(dissim_matrix )) == 0) stop("ERROR: dissim_matrix  has no rownames. Please define rownames as ligand-receptor pairs")
        if(length(colnames(dissim_matrix )) == 0) stop("ERROR: dissim_matrix  has no rownames. Please define colnames as ligand-receptor pairs")
        
        # check that rownames and colnames of the dissim_matrix  are identical
        if(!identical(rownames(dissim_matrix )
                      ,colnames(dissim_matrix ))) stop("ERROR: rownames and colnames of the dissim_matrix  are not identical. Please make sure they are identical")
        
        # check that the lrp_clusters is named
        if(length(names(lrp_clusters)) == 0) stop("ERROR: lrp_clusters is not a named vector. Please nake sure the names are ligand-receptor pairs")
        
        # check that rownames and colnames of the dissim_matrix and names of lrp_clusters are identical
        if(!identical(rownames(dissim_matrix)
                      ,names(lrp_clusters))) stop("ERROR: rownamesod dissim_matrix and names of lrp_clusters are not identical. Please make identical")
        if(!identical(colnames(dissim_matrix)
                      ,names(lrp_clusters))) stop("ERROR: colnamesod dissim_matrix and names of lrp_clusters are not identical. Please make identical")
        
        # assign colours to clusters
        if(cluster_colors == "") cluster_colors <- structure(rainbow(length(levels(as.factor(lrp_clusters)))))
        names(cluster_colors) <- levels(as.factor(lrp_clusters))
        
        # if there are unclustered ligand-receptor pairs, colour them with a defined colour
        if(as.factor(0) %in% names(cluster_colors)) cluster_colors[names(cluster_colors) == as.factor(0)] <- unclustered_LRP_color
        
        # mark the as "unclustered"
        if(max(lrp_clusters) != 0) {
                annotation_legend_param <-  list(labels = c("unclustered"
                                                            ,1:max(lrp_clusters))
                )} else { # all ligand-receptor pairs remained unclustered
                        annotation_legend_param <-  list(labels = "unclustered"
                        )
                }
        
        # define column annotation
        haCol <- ComplexHeatmap::HeatmapAnnotation(
                df = data.frame(cluster = as.vector(lrp_clusters))
                ,col = list(cluster = cluster_colors)
                ,annotation_legend_param = annotation_legend_param
        )
        
        # heatmap ordered by clusters
        h1 <- ComplexHeatmap::Heatmap(
                dissim_matrix
                ,cluster_rows = F
                ,row_order = rownames(dissim_matrix)[order(as.numeric(lrp_clusters))]
                ,cluster_columns = F
                ,name = "dissimilarity coefficient"
                ,column_order = colnames(dissim_matrix)[order(as.numeric(lrp_clusters))]
                ,top_annotation = haCol
                ,show_column_names = show_column_names
                ,show_row_names = show_row_names
                ,heatmap_legend_param = list(legend_direction = legend_direction
                                             ,legend_width = legend_width
                                             ,title_position = title_position
                )
                ,...
        )
        
        heatmap_clusters <- ComplexHeatmap::draw(
                h1
                ,column_title = column_title
                ,column_title_side = column_title_side
                ,heatmap_legend_side = heatmap_legend_side
        )
}

#' @rdname plot_cluster_UMAP
#' @export
#' 
#' @title 
#' Plots UMAP of ligand-receptor pairs
#' 
#' @description 
#' \code{plot_cluster_UMAP} plots UMAP of ligand-receptor pairs colored and shaped by cluster.
#' 
#' @author
#' Maria Solovey \email{maria.solovey@helmholtz-muenchen.de}
#'
#' @param       ligand_receptor_pair_df character dataframe: data frame with columns "pair", "ligand", "ligand_complex_composition", "receptor", "receptor_complex_composition". 
#' Column "pair" contains values in a form "ligand:receptor", i.e. ligand being at the first place, receptor being at the second place, e.g. "TNFSF13:TNFRSF17". 
#' Column "ligand" contains ligand names, e.g. "TNFSF13".  
#' Column ligand_complex_composition" if ligand is a complex (e.g. "aXb2_complex"), 
#' contains genes in the ligand complex separated with a comma, e.g. "ITGAX,ITGB2", else contains empty string "". 
#' Column "receptor" contains receptor names, e.g. "TNFRSF17". 
#' Column "receptor_complex_composition" if receptor is a complex (e.g. "NKG2D_II_receptor"), 
#' contains genes in the receptor complex separated with a comma, e.g. "KLRK1,HCST", else contains empty string "".
#' 
#' @param       dissim_matrix numeric matrix: pairwise dissimilarity matrix. Note that rownames and colnames should be defined.
#' 
#' @param       lrp_clusters numeric vector: cluster assignment for each ligand-recetor pair.
#'  Note that the the vector should be named and ordered in the same way as the rownames and colnames of the dissim_matrix.
#'  
#' @param       seed random seed for UMAP. Default value: 100.
#' 
#' @param       cluster_colors character string vector: vector with colours for each cluster. Default value: "".
#' 
#' @param       unclustered_LRP_color character string: colour for ligand-receptor pairs that remained unclustered (cluster assignment 0). Default value: "black".
#' 
#' @param       cluster_shape numeric vector: vector with shapes for each cluster. Default value: "".
#' 
#' @param       point_size numeric: size of points. Default value: 2.
#' 
#' @param       title character string: title of the heatmap plot. Argument of the draw function. Default value: "umap of ligand-receptor pairs".
#' 
#' @param       legend_position character string: argument of the theme function. Default value: "bottom" guide_legend.
#' 
#' @param       legend_direction character string: argument of the guide_legend function. Default value: "horizontal".
#' 
#' @param       legend_title_position character string: argument of the guide_legend function. Default value: "top".
#' 
#' @param       legend_label_position character string: argument of the guide_legend function. Default value: "bottom".
#' 
#' @param       legend_byrow logical: argument of the guide_legend function. Default value: TRUE.
#' 
#' @return      UMAP of ligand-receptor pairs coloured and shaped by cluster
plot_cluster_UMAP <- function(ligand_receptor_pair_df
                              ,dissim_matrix
                              ,lrp_clusters
                              ,seed = 100
                              ,cluster_colors = ""
                              ,unclustered_LRP_color = "black"
                              ,cluster_shape = ""
                              ,point_size = 2
                              ,title = "umap of ligand-receptor pairs"
                              ,legend_position="bottom"
                              ,legend_direction = "horizontal"
                              ,legend_title_position = "top"
                              ,legend_label_position= "bottom"
                              ,legend_byrow=TRUE
){
        
        # set seed
        set.seed(seed)
        
        # define n_neighbors parameter
        n_neighbors <- umap::umap.defaults$n_neighbors
        if(n_neighbors >= nrow(ligand_receptor_pair_df)) n_neighbors <- round(nrow(ligand_receptor_pair_df) /2)
        
        # check wether there are enough ligand-receptor pairs to plot a UMAP
        if(n_neighbors > 1){
                custom.settings <-  umap.defaults
                custom.settings$input <- "dist"
                custom.settings$n_neighbors <- n_neighbors
                
                # calculate the UMAP
                my.umap <- umap::umap(dissim_matrix
                                      ,config = custom.settings)
                
                # make data frame for UMAP plotting
                dfForUmap <- data.frame(x=my.umap$layout[,1]
                                        ,y=my.umap$layout[,2]
                                        ,cluster=lrp_clusters
                )
                
                # factorize clusters
                dfForUmap$cluster <- as.factor(dfForUmap$cluster)
                
                # assign colours to clusters
                if(cluster_colors == "") cluster_colors <- structure(rainbow(length(levels(as.factor(lrp_clusters)))))
                names(cluster_colors) <- levels(as.factor(lrp_clusters))
                
                # if there are unclustered ligand-receptor pairs, colour them with a defined colour
                if(as.factor(0) %in% names(cluster_colors)) cluster_colors[names(cluster_colors) == as.factor(0)] <- unclustered_LRP_color
                
                # define cluster shape
                if(cluster_shape == "") cluster_shape <- as.integer(names(cluster_colors)) %% 5 +21
                names(cluster_shape) <- names(cluster_colors)
                
                # plot UMAP
                umap.clusters <- ggplot2::ggplot(dfForUmap
                                                 ,aes(x=x
                                                      ,y=y
                                                      ,color=cluster
                                                      ,shape = cluster
                                                 )
                ) +
                        geom_point(size=point_size) +
                        scale_color_manual(name = "Cluster"
                                           ,labels = c("unclustered"
                                                       ,1:max(lrp_clusters))
                                           ,values = cluster_colors) +
                        scale_shape_manual(name = "Cluster"
                                           ,labels = c("unclustered"
                                                       ,1:max(lrp_clusters))
                                           ,values = cluster_shape) +
                        ylab("umap2") +
                        xlab("umap1") +
                        ggtitle(title) +
                        theme(legend.position=legend_position) +
                        guides(fill=guide_legend(nrow=2
                                                 ,direction = legend_direction
                                                 ,title.position = legend_title_position
                                                 ,label.position= legend_label_position
                                                 ,byrow=legend_byrow
                        )
                        )
                print(umap.clusters)
        } else print("WARNING: too few ligand-receptor pairs. UMAP can not be plotted.")
}

#' @rdname make_pattern_matrix
#' @export
#' 
#' @title 
#' Makes a weight matrix for pattern of interest
#'
#' @description 
#' \code{make_pattern_matrix} makes a weight matrix for pattern of interest.
#' The rows and columns of the matrix will contain all nodes.
#' The rows are considered as sending nodes and the columns are considered as receiving nodes. 
#' The value of 1 represents an edge (i.e. cell types are communicating) and 0 represents no edge (i.e. cell types are not communicating). 
#' 
#' @author
#' Maria Solovey \email{maria.solovey@helmholtz-muenchen.de}
#'
#' @param communicating_nodes character string vector: in form of c("cellType1_to_celltype2", "cell_type1_to_celltype3", ...), e.g. c("HSC_to_Tcells", "HSC_to_bcells", ...).
#' 
#' @param nodes character string vector: of length n with all cell types in the data.
#' 
#' @return     weight matrix: numeric n x n matrix with values 0 or 1. Rows are nodes, columns are nodes. 
#' Rows are regarded as sending cell type, columns are regarded as receiving cell types. 1 defines communication, 0 defines no communication
make_pattern_matrix <- function(communicating_nodes
                                ,nodes
){
        # make a dummy matrix
        pattern_w_matrix <- matrix(0
                                   ,nrow = length(nodes)
                                   ,ncol = length(nodes)
        )
        dimnames(pattern_w_matrix) <- list(nodes,nodes)
        
        # define sending nodes
        sending_nodes <- gsub( "_.*$", "", communicating_nodes)
        print("sending nodes are:")
        print(sending_nodes)
        
        # check correctness of the node names
        if(!all(sending_nodes %in% nodes)) stop("ERROR: cell type names in 'communicating_nodes' argument are not correct. Please check spelling.")
        
        # define receiving nodes
        receiving_nodes <- gsub('.*to_', "", communicating_nodes)
        print("recieving nodes are:")
        print(receiving_nodes)
        
        
        # check correctness of the node names
        if(!all(receiving_nodes %in% nodes)) stop("ERROR: cell type names in 'communicating_nodes' argument are not correct. Please check spelling.")
        
        # fill in 1 to communication nodes
        for(i in 1:length(communicating_nodes)){
                #print(i)
                pattern_w_matrix[sending_nodes[i],receiving_nodes[i]] <- 1
        }
        
        # return
        return(pattern_w_matrix)
}

#' @rdname convert_CellPhoneDB_output
#' @export
#' 
#' @title 
#' Converts CellPhoneDB output matrix into an array of weighted adjacency matrices
#'
#' @description 
#' \code{convert_CellPhoneDB_output} converts CellPhoneDB output matrix (significant_means.txt) into an array of weighted adjacency matrices.
#' The CellPhoneDB weight matrix (significant_means.txt) contains ligand-receptor pairs in the rows and all pairs of cell types in the columns. 
#' COMUNET transforms it into a stack of weight matrices (i.e. weight array), one weight matrix per ligand-receptor pair. 
#' If CellPhoneDB output contains \code{m} non-empty ligand-receptor pairs, then the number of matrices in the stack will be equal to \code{m}. 
#' Each such weight matrix has all cell types in the rows and all cell types in the columns. 
#' So if there are \code{n} cell types in the data, each weight matrix is an \code{n x n} matrix. 
#' The rows of each matrix represent sending cell types, i.e. cell types that express a ligand of a ligand-receptor pair, 
#' or just partner A of an undirected interacting pair (e.g. pair of adhesion molecules). 
#' The columns of each matrix represent receiving cell types, i.e. cell types that express a receptor of a ligand-receptor pair, 
#' or partner B of an undirected interacting pair (e.g. pair of adhesion molecules). 
#' The NA values of the CellPhoneDB output are substituted with 0.
#' 
#' @author
#' Maria Solovey \email{maria.solovey@helmholtz-muenchen.de}
#'
#' @param       CellPhoneDB_output significant_means.txt file in the output of CellPhoneDB.
#' 
#' @param       complex_input complex_input.csv file. Download from CellPhoneDB.
#' 
#' @param       gene_input gene_input.csv file. Download from CellPhoneDB.
#' 
#' @return list of:
#'               weight_array: 3D numeric array                       array of weighted adjacency matrices with dimensions [number of nodes, number of nodes, number of ligand-receptor pairs]. 
#' Note that the function filters out empty weight arrays and the corresponding ligand-receptor pairs
#'                       First dimension: sending nodes (note: this dimension has all possible nodes, even if some of them are silent for a particular ligand-receptor pair)
#'                       Second dimension: receiving nodes (note: this dimension has all possible nodes, even if some of them are silent for a particular ligand-receptor pair)
#'                       Third dimension: ligand-receptor pairs
#'                       Note that the weight_array should contain dimnames: dimnames = list(nodes, nodes, ligand-receptor pairs)
#'                       
#'               ligand_receptor_pair_df: string dataframe          data frame with columns "pair", "ligand", "ligand_complex_composition", "receptor", "receptor_complex_composition"
#'                        "pair" contains values in a form "ligand:receptor", i.e. ligand being at the first place, receptor being at the second place, e.g. "TNFSF13:TNFRSF17"
#'                        "ligand" contains ligand names, e.g. "TNFSF13"
#'                        "ligand_complex_composition" if ligand is a complex (e.g. "aXb2_complex"), contains genes in the ligand complex separated with a comma, e.g. "ITGAX,ITGB2", else contains ""
#'                        "receptor" contains receptor names, e.g. "TNFRSF17"
#'                        "receptor_complex_composition" if receptor is a complex (e.g. "NKG2D_II_receptor"), contains genes in the receptor complex separated with a comma, e.g. "KLRK1,HCST", else contains ""
#'               nodes: character string vector                            a vector with all cell types in the data
convert_CellPhoneDB_output <- function(CellPhoneDB_output
                                       ,complex_input
                                       ,gene_input
){
        # make ligand-receptor pair list
        # ligand_receptor_pair_df is data frame with columns "pair", "ligand", "ligand_complex_composition", "receptor", "receptor_complex_composition"
        # "pair" contains values in a form "ligand:receptor", i.e. ligand being at the first place, receptor being at the second place, e.g. "TNFSF13:TNFRSF17"
        # "ligand" contains ligand names, e.g. "TNFSF13"
        # "ligand_complex_composition" if ligand is a complex (e.g. "aXb2_complex"), contains genes in the ligand complex separated with a comma, e.g. "ITGAX,ITGB2", else contains ""
        # "receptor" contains receptor names, e.g. "TNFRSF17"
        # "receptor_complex_composition" if receptor is a complex (e.g. "NKG2D_II_receptor"), contains genes in the receptor complex separated with a comma, e.g. "KLRK1,HCST", else contains ""
        
        # ligands and receptors: each pair will have a form "ligand:receptor"
        
        # indices of receptors
        idx_receptor <- as.matrix(data.frame(a = CellPhoneDB_output$receptor_a
                                             ,b = CellPhoneDB_output$receptor_b))
        rownames(idx_receptor) <- CellPhoneDB_output$interacting_pair
        
        # which pairs are directed
        idx_directed <- rowSums(idx_receptor) == 1
        
        # make empty ligand_receptor_pair_df data.frame
        ligand_receptor_pair_df <- as.data.frame(matrix(NA
                                                        ,nrow = nrow(CellPhoneDB_output)
                                                        ,ncol = 5
        ))
        colnames(ligand_receptor_pair_df) <- c("pair"
                                               ,"ligand"
                                               ,"ligand_complex_composition"
                                               ,"receptor"
                                               ,"receptor_complex_composition"
        )
        
        # put in ligand and receptors
        split_LRP <- {
                # check if any ligands or receptor complexes are written with "_"
                # we will split pairs by "_", so any "_" in a complex name will interfere and we want to replace it with " "
                which_to_replace <- function(idx_vector # all positions of "_" in the pair name
                                             ,partner # which partner is a complex
                                             ,how_many # how many "_" does this complex contain (i.e. how many do we want to replace with " ")
                ){
                        if(partner == "partner_a") {
                                idx_vector <- idx_vector[1:how_many] # take first "how_many" values
                        } else {
                                idx_vector <- idx_vector[(length(idx_vector) - how_many + 1
                                ):length(idx_vector)] # take last "how_many
                        }
                } # the function returns a vector with exact positions of "_" in the pair name that should be replaced with " "
                
                # replace "_" with " " in the complex names in the CellPhoneDB_output$interacting_pair
                for(partner in c("partner_a"
                                 ,"partner_b"
                )
                ){
                        # which complexes in partern_a or in pertner_b colomn contain "_"
                        idx_complex <- grep("_"
                                            ,unlist(CellPhoneDB_output[,partner])
                        )
                        # replace "_" for space in each position
                        if(length(idx_complex) != 0){
                                for(lrp in idx_complex) {
                                        # how many "_" does a particular partner_a or partner_b complex contain
                                        how_many <- length(unlist(gregexpr(pattern = "_"
                                                                           ,unlist(CellPhoneDB_output[lrp,partner]
                                                                           ))
                                        )
                                        )
                                        # find all positions of "_" in the interacting_pair
                                        idx_char <- unlist(gregexpr(pattern = "_"
                                                                    ,CellPhoneDB_output$interacting_pair[lrp]
                                        )
                                        )
                                        # take only those that correspond to the prtner_a or prtner_b
                                        idx_char <- which_to_replace(idx_vector = idx_char
                                                                     ,partner = partner
                                                                     ,how_many = how_many)
                                        # replace them for " "
                                        for(idx in idx_char){
                                                substr(CellPhoneDB_output$interacting_pair[lrp]
                                                       ,idx
                                                       ,idx) <- " "
                                        }
                                }
                        }
                }
                
                # update rownames
                rownames(CellPhoneDB_output) <- CellPhoneDB_output$interacting_pair
                rownames(idx_receptor) <- CellPhoneDB_output$interacting_pair
                
                # split by "_"
                split_LRP <- strsplit(CellPhoneDB_output$interacting_pair
                                      ,"_")
                
                # check splitting correctness
                nr_partners <- sapply(1:nrow(CellPhoneDB_output)
                                      ,function(i) length(split_LRP[[i]])
                )
                if(any(nr_partners!=2)) stop("ERROR: Not all pairs contain 2 partners.")
                # translate to data frame
                split_LRP <- matrix(unlist(split_LRP)
                                    ,ncol = 2
                                    ,byrow = T)
        }
        ligand_receptor_pair_df$ligand <- as.vector(split_LRP[,1])
        ligand_receptor_pair_df$receptor <- as.vector(split_LRP[,2])
        
        # check if any ligand and receptor pairs are swapped
        idx_swapped <- !idx_receptor[,2] & idx_directed
        ligands <- ligand_receptor_pair_df$receptor[idx_swapped]
        ligand_receptor_pair_df$receptor[idx_swapped] <- ligand_receptor_pair_df$ligand[idx_swapped]
        ligand_receptor_pair_df$ligand[idx_swapped] <- ligands
        
        # put in pairs
        ligand_receptor_pair_df$pair <- paste0(ligand_receptor_pair_df$ligand
                                               ,":"
                                               ,ligand_receptor_pair_df$receptor)
        CellPhoneDB_output$sorted_pair <- ligand_receptor_pair_df$pair
        
        # put in complex composition
        ligand_receptor_pair_df$ligand_complex_composition <- suppressWarnings( # suppresses warning: In `[.data.frame`(..., drop = FALSE) : 'drop' argument will be ignored
                sapply(ligand_receptor_pair_df$ligand
                       ,function(g) {
                               genes <- complex2genes(g
                                                      ,complex_input = complex_input
                                                      ,gene_input = gene_input)
                               paste(unlist(genes)
                                     ,collapse = ",")
                       }
                ))
        ligand_receptor_pair_df$receptor_complex_composition <- suppressWarnings( # suppresses warning: In `[.data.frame`(..., drop = FALSE) : 'drop' argument will be ignored
                sapply(ligand_receptor_pair_df$receptor
                       ,function(g) {
                               genes <- complex2genes(g
                                                      ,complex_input = complex_input
                                                      ,gene_input = gene_input)
                               paste(unlist(genes)
                                     ,collapse = ",")
                       }
                )
        )
        
        # assign rownames as in CellPhoneDB_output$interacting_pair
        rownames(ligand_receptor_pair_df) <- CellPhoneDB_output$interacting_pair
        
        # check NAs
        check_NAs(ligand_receptor_pair_df)
        
        # make a key of nodes (in order to convert a vector of pairwise nodes to an adjacency matrix)
        node_combinations <- colnames(CellPhoneDB_output)[grepl("[|]"
                                                                ,colnames(CellPhoneDB_output))]
        node_key <- matrix(unlist(strsplit(node_combinations
                                           ,"[|]")
        )
        , ncol = 2
        , byrow = T)
        dimnames(node_key) <- list(node_combinations
                                   ,c("a", "b"))
        
        nodes <- unique(c(node_key[,"a"]
                          ,node_key[,"b"])
        )
        # make weight array
        # make a dummy array of weights
        weight_array <- array(0
                              ,dim = c(length(nodes) # sending nodes
                                       ,length(nodes) # recieving nodes
                                       ,nrow(ligand_receptor_pair_df) # ligand-receptor pairs
                              )
        )
        dimnames(weight_array) <- list(unique(c(node_key[,"a"]
                                                ,node_key[,"b"])
        )
        ,unique(c(node_key[,"a"]
                  ,node_key[,"b"])
        )
        ,ligand_receptor_pair_df$pair
        )
        # populate array of weights
        idx_swapped <- which(idx_swapped)
        for(layer in CellPhoneDB_output$interacting_pair){
                for(nodePair in node_combinations){
                        if(!is.na(CellPhoneDB_output[layer
                                                     ,nodePair])) {
                                idxLayerCPDB <- which(rownames(ligand_receptor_pair_df) == layer)
                                # control swapped pairs
                                if(idxLayerCPDB %in% idx_swapped) {
                                        ligand <- 2
                                        receptor <- 1
                                } else {
                                        ligand <- 1
                                        receptor <- 2
                                }
                                weight_array[node_key[nodePair,][ligand] # partner a (ligand)
                                             ,node_key[nodePair,][receptor] # partner b (receptor)
                                             ,as.character(ligand_receptor_pair_df[layer,"pair"])
                                             ] <- as.numeric(CellPhoneDB_output[idxLayerCPDB 
                                                                                ,nodePair])
                        }
                }
        }
        
        # filter empty weight matrices
        filtered_data <- filter_empty_layers(weight_array = weight_array
                                             ,ligand_receptor_pair_df = ligand_receptor_pair_df
        )
        ligand_receptor_pair_df <- filtered_data$ligand_receptor_pair_df
        weight_array_filtered <- filtered_data$weight_array
        rm(filtered_data)
        
        
        # print a warning of number of ligand-receptor pairs is lower than 10
        if(nrow(ligand_receptor_pair_df) < 10) print("WARNING: there are very few ligand-receptor pairs!")
        
        # return
        result <- list(weight_array = weight_array_filtered
                       ,ligand_receptor_pair_df = ligand_receptor_pair_df
                       ,nodes = nodes
        )
}

#' @rdname lrp_clustering
#' @export
#' 
#' @title 
#' Clusters interacting protein pairs
#'
#' @description 
#' \code{lrp_clustering} clusters interacting protein pairs by performing hierarchical clustering of the dissimilarity matrix of layers 
#' (by default, the \code{hclust} R function with the \code{average} agglomeration method is used) and the results can be visualized as a heatmap or a UMAP plot (McInnes et al., 2018). 
#' The number of clusters is estimated using the  \code{cutreeHybrid} R function (package dynamicTreeCut, version 1.63-1(Langfelder et al., 2008)) 
#' with \code{deepSplit} equal to 0 and default \code{minClusterSize} equal to 6 by default (both can be adjusted). 
#' For each cluster, a graph that represents the “average” pattern is built by averaging the adjacency matrices and the (see section 2.1) of the nodes of all the graphs in the cluster.
#' 
#' @author
#' Maria Solovey \email{maria.solovey@helmholtz-muenchen.de}
#'
#' @param       weight_array numeric array: array of weighted adjacency matrices with dimensions [number of nodes, number of nodes, number of ligand-receptor pairs].
#' 
#' @param       ligand_receptor_pair_df character dataframe: data frame with columns "pair", "ligand", "ligand_complex_composition", "receptor", "receptor_complex_composition". 
#' Column "pair" contains values in a form "ligand:receptor", i.e. ligand being at the first place, receptor being at the second place, e.g. "TNFSF13:TNFRSF17". 
#' Column "ligand" contains ligand names, e.g. "TNFSF13".  
#' Column ligand_complex_composition" if ligand is a complex (e.g. "aXb2_complex"), 
#' contains genes in the ligand complex separated with a comma, e.g. "ITGAX,ITGB2", else contains empty string "". 
#' Column "receptor" contains receptor names, e.g. "TNFRSF17". Column "receptor_complex_composition" if receptor is a complex (e.g. "NKG2D_II_receptor"), 
#' contains genes in the receptor complex separated with a comma, e.g. "KLRK1,HCST", else contains empty string "".
#' 
#' @param       nodes character string vector: a vector with all cell types in the data.
#' 
#' @param       dissimilarity function: dissimilarity function. Default value is d_normWeightDiff.
#' 
#' @param       deep_split numeric: defines the deepSplit parameter of the \code{cutreeHybrid} function. Default value: 0.
#' 
#' @param       min_cluster_size numeric: minimum number of ligand-receptor pairs per cluster. Default value: 6.
#' 
#' @return a list of:
#'               dissim_matrix: numeric matrix           pairwise dissimilarity between all ligand-receptor pairs
#'               clusters: numeric vector                cluster assignment for each ligand-receptor pair
#'               weight_array_by_cluster: numeric array   array of weighted adjacency matrices with dimensions [number of nodes, number of nodes, number of ligand-receptor clusters]. 
#' Each edge weight in a cluster is calculated as the mean of the weights of this edge in all ligand-receptor pairs that belong to this cluster.
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
                        
                        # calcualte degree arrays for each community (put it as well as an array)
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

#' @rdname pattern_search
#' @export
#' 
#' @title 
#' Performs pattern search
#'
#' @description 
#' \code{pattern_search} performs pattern search.
#' A pattern of interest is defined with a binary adjacency matrix, 1 indicating the presence of an edge and 0 the absence. 
#' The adjacency matrices of the layers are also binarized based on the presence/absence of the edges and then the dissimilarities with the user-specified adjacency matrix are calculated. 
#' The output is a list of pairs of interacting partners sorted by increasing dissimilarity with the user-specified pattern.
#' 
#' @author
#' Maria Solovey \email{maria.solovey@helmholtz-muenchen.de}
#'
#' @param       pattern_adj_matrix numeric matrix: an adjacency matrix for the pattern of interest.
#' 
#' @param       weight_array numeric array: array of weighted adjacency matrices with dimensions [number of nodes, number of nodes, number of ligand-receptor pairs].
#' 
#' @param       ligand_receptor_pair_df character dataframe: data frame with columns "pair", "ligand", "ligand_complex_composition", "receptor", "receptor_complex_composition". 
#' Column "pair" contains values in a form "ligand:receptor", i.e. ligand being at the first place, receptor being at the second place, e.g. "TNFSF13:TNFRSF17". 
#' Column "ligand" contains ligand names, e.g. "TNFSF13".  
#' Column "ligand_complex_composition" if ligand is a complex (e.g. "aXb2_complex"), 
#' contains genes in the ligand complex separated with a comma, e.g. "ITGAX,ITGB2", else contains empty string "". 
#' Column "receptor" contains receptor names, e.g. "TNFRSF17". 
#' Column "receptor_complex_composition" if receptor is a complex (e.g. "NKG2D_II_receptor"), 
#' contains genes in the receptor complex separated with a comma, e.g. "KLRK1,HCST", else contains empty string "".
#' 
#' @param       nodes character string vector: a vector with all cell types in the data.
#' 
#' @param       dissimilarity function: dissimilarity measure. Default value: d_normWeightDiff.
#' 
#' @return dataframe:                              The data frame is sorted by increasing dissimilarity (i.e. similar patterns a the top
#'               pair: character string vector             ligand-receptor pair names
#'               dissimilarity: numeric vector  dissimilarity of ligand-receptor pairs to the pattern of interest. Note here that for the calculation of the dissimilarity of a weight matrix to the pattern of interest, we don't take into consideration the actual weight of the edges in this weight matrix, but just check wether the edge is present (i.e. has a non-zero value) or absent (i.e. has a zero value).
pattern_search <- function(pattern_adj_matrix # in rows are sending nodes, in columns are receiving nodes, 0 for no communication, 1 for communication
                           ,weight_array
                           ,ligand_receptor_pair_df
                           ,nodes
                           ,dissimilarity = d_normWeightDiff
){
        # check whether pattern_adj_matrix contains all nodes
        if((!identical(rownames(pattern_adj_matrix)
                       ,nodes
        ) ) |
        (!identical(colnames(pattern_adj_matrix)
                    ,nodes)
        )
        ) stop("ERROR: Row names or column names of the pattern_adj_matrix are not identical to the nodes.")
        
        # convert pattern_adj_matrix to array, We need it to pass two arrays to the dissimilarity function
        pattern_adj_array <- array(pattern_adj_matrix
                                   ,dim = list(dim(pattern_adj_matrix)[[1]]
                                               ,dim(pattern_adj_matrix)[[2]]
                                               ,1)
        )
        
        # binarize weight_array
        binarized_weight_array <- (weight_array != 0) *1
        
        # calculate dissimilarity to the each LRP
        dissim <- sapply(1:dim(binarized_weight_array)[[3]]
                         ,function(k){
                                 dissimilarity(pattern_adj_array
                                               ,binarized_weight_array[,,k]
                                 )
                         })
        
        # result
        result <- data.frame(pair = dimnames(binarized_weight_array)[[3]]
                             ,dissimilarity = dissim
        )
        rownames(result) <- result$pair
        
        # sort by increasing dissimilarity
        result <- result[order(result$dissimilarity, decreasing = F),]
        
        # return result
        return(pattern_result = result)
}

#' @rdname comparative_analysis
#' @export
#' 
#' @title 
#' Performs comparative analysis
#'
#' @description 
#' \code{comparative_analysis} performs comparative analysis between two conditions (e.g. two samples).
#' Given the results of a cell-cell communication analysis of two datasets including the same cell types, 
#' COMUNET estimates the differences in the intercellular communication patterns between corresponding interacting partners in the two datasets.
#' 
#' @author
#' Maria Solovey \email{maria.solovey@helmholtz-muenchen.de}
#'
#' @param cond1_weight_array numeric array: array of weighted adjacency matrices for condition 1. 
#' Note that the function filters out empty weight arrays and the corresponding ligand-receptor pairs.
#' 
#' @param cond2_weight_array numeric array: array of weighted adjacency matrices for condition 2. 
#' Note that the function filters out empty weight arrays and the corresponding ligand-receptor pairs.
#' 
#' @param cond1_ligand_receptor_pair_df character dataframe: data frame with columns "pair", "ligand", "ligand_complex_composition", "receptor", "receptor_complex_composition" for condition 1.
#' 
#' @param cond2_ligand_receptor_pair_df character dataframe: data frame with columns "pair", "ligand", "ligand_complex_composition", "receptor", "receptor_complex_composition" for condition 2.
#' 
#' @param cond1_nodes character string vector: a vector with all cell types in the data. 
#' Please note that comparative_analysis expects nodes (i.e. cell types) to be the same in both conditions.
#' 
#' @param cond2_nodes character string vector: a vector with all cell types in the data. Default value: cond1_nodes. 
#' 
#' Please note that comparative_analysis expects nodes (i.e. cell types) to be the same in both conditions.
#' @param cond1_name character string: sample name for condition 1.
#' 
#' @param cond2_name character string: sample name for condition 2.
#' 
#' @param dissimilarity function: dissimilarity function. Default value is d_normWeightDiff.
#' 
#' @return a list of:
#' sorted_LRP_df: data.frame with columns:
#' 	pair: character string                          names of ligand-receptor pairs in the same form as they are in ligand_receptor_pair_df$pair
#' 	presence: character string                      whether the ligand-receptor pair is present in both conditions ("shared") or only in one of them.
#' 	dissimilarity: numeric                 dissimilarity value between the topology of ligand-receptor pair graph in two conditions. The smaller the dissimilarity value, the more similar is the graph topology between the two conditions. If a ligand-receptor pair is present only in one of the conditions, the dissimilarity is equal to 1.
#' dissim_cond1_cond2: numeric matrix            pairwise dissimilarity between all ligand-receptor pairs in the two conditions (condition 1 in rows, condition 2 in columns)
comparative_analysis <- function(cond1_weight_array
                                 ,cond2_weight_array
                                 ,cond1_ligand_receptor_pair_df
                                 ,cond2_ligand_receptor_pair_df
                                 ,cond1_nodes
                                 ,cond2_nodes = cond1_nodes
                                 ,cond1_name
                                 ,cond2_name
                                 ,dissimilarity = d_normWeightDiff 
){
        
        # check if the nodes in both conditions are the same
        if(!identical(cond1_nodes, cond2_nodes)) stop("ERROR: one of the conditions might contain cell populations that are missing in the other. Please make sure both conditions contain same cell populations.")
        
        # stratify LRPs by shared or not shared
        all_LRPs <- union(cond1_ligand_receptor_pair_df$pair
                          ,cond2_ligand_receptor_pair_df$pair
        )
        
        shared_LRP <- intersect(cond1_ligand_receptor_pair_df$pair
                                ,cond2_ligand_receptor_pair_df$pair
        )
        
        only_cond1_LRPs <- setdiff(cond1_ligand_receptor_pair_df$pair
                                   ,cond2_ligand_receptor_pair_df$pair
        )
        
        only_cond2_LRPs <- setdiff(cond2_ligand_receptor_pair_df$pair
                                   ,cond1_ligand_receptor_pair_df$pair
        )
        
        # calculate dissimilarity
        dissim_cond1_cond2 <- calc_dissim_matrix(weight_array1 = cond1_weight_array
                                                 ,weight_array2 = cond2_weight_array
                                                 ,dissimilarity = dissimilarity
        )
        
        # make sorted LRP list
        sorted_LRP_df <- as.data.frame(t(cbind(sapply(all_LRPs
                                                      ,function(lrp){
                                                              if((length(shared_LRP) != 0) & (lrp %in% shared_LRP)){
                                                                      c(lrp
                                                                        ,"shared"
                                                                        ,unlist(dissim_cond1_cond2[lrp,lrp]))
                                                              } else if((length(only_cond1_LRPs) != 0) & (lrp %in% only_cond1_LRPs)) {
                                                                      c(lrp
                                                                        ,paste("only"
                                                                               ,cond1_name)
                                                                        ,1)
                                                              } else if((length(only_cond2_LRPs) != 0) & (lrp %in% only_cond2_LRPs)) {
                                                                      c(lrp
                                                                        ,paste("only"
                                                                               ,cond2_name)
                                                                        ,1)
                                                              }
                                                      })
        )))
        colnames(sorted_LRP_df) <- c("pair"
                                     ,"presence"
                                     ,"dissimilarity")
        rownames(sorted_LRP_df) <- sorted_LRP_df$pair
        
        # sort by accending dissimilarity
        sorted_LRP_df$dissimilarity <- as.numeric(sorted_LRP_df$dissimilarity)
        sorted_LRP_df <- sorted_LRP_df[order(sorted_LRP_df$dissimilarity
                                             ,decreasing = T),]
        
        # sort by presence
        sorted_LRP_df$presence <- factor(sorted_LRP_df$presence
                                         ,levels = c("shared"
                                                     ,paste("only"
                                                            ,cond1_name)
                                                     ,paste("only"
                                                            ,cond2_name)
                                         )
                                         ,labels = c("shared"
                                                     ,paste("only"
                                                            ,cond1_name)
                                                     ,paste("only"
                                                            ,cond2_name)
                                         )
                                         ,ordered = T)
        sorted_LRP_df <- sorted_LRP_df[order(sorted_LRP_df$presence
                                             ,decreasing = F),]
        
        return(list(sorted_LRP_df = sorted_LRP_df
                    ,dissim_cond1_cond2 = dissim_cond1_cond2)
        )
        
}