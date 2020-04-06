# Functions ####

# make a function to check NAs
# INPUT: 
#       objectToCheck: a vector or matrix containing any values
# OUPUT:
#       a warning message is case the object contains NAs
check_NAs <- function(objectToCheck){
        if(!all(!is.na(objectToCheck))) stop("WARNING: the object contains NAs.")
}

step_func <- function(x) {as.integer(x != 0)}

# prevent dropping dimensions
`[` <- function(...) base::`[`(...,drop=FALSE)

# filter out layers with no edges
filter_empty_layers <- function(weight_array
                              ,ligandReceptorPairDF
){
        # check which layers have no edges
        noEdgeLayers <- sapply(1:nrow(ligandReceptorPairDF)
                               ,function(layer){
                                       sum(weight_array[,,layer]) == 0
                               })
        names(noEdgeLayers) <- dimnames(weight_array)[3]
        # delete these layers form ligandReceptorPairDF
        ligandReceptorPairDF <- ligandReceptorPairDF[!noEdgeLayers,]
        # delete these layers from the array of the adjacency matrices
        weight_array <- weight_array[,,!noEdgeLayers]
        # return
        return(list(ligandReceptorPairDF = ligandReceptorPairDF
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
        myDegreeArray <- array(NA
                               ,dim = c(N,3,M)
        )
        
        # add dimension names
        if(M == 1){
                dimnames(myDegreeArray) <- list(dimnames(weight_array)[[1]]
                                                ,c("in", "out", "delta")
                                                ,1
                )        
        } else {
                dimnames(myDegreeArray) <- list(dimnames(weight_array)[[1]]
                                                ,c("in", "out", "delta")
                                                ,dimnames(weight_array)[[3]]
                )
        }
        
        myDegreeArray[,"in",] <- colSums(weight_array)
        if(M==1){
                myDegreeArray[,"out",] <- rowSums(weight_array)
        }else{
                myDegreeArray[,"out",] <- colSums(aperm(weight_array
                                                        ,c(2,1,3)))
        }
        myDegreeArray[,"delta",] <- myDegreeArray[,"out",] - myDegreeArray[,"in",]
        # check NAs
        check_NAs(myDegreeArray)
        # return myDegreeArray
        myDegreeArray
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
        unionOfEdges <- apply(adj1, 1:2, step_func) | apply(adj2, 1:2, step_func)
        # calculate the number of edges in the union
        n_unionOfEdges <- sum(unionOfEdges)
        
        if(n_unionOfEdges == 0) myDissim <- 1 else {
                # calculate the sum of normalised weight differences
                # matrix of absolute differences of weights for each edge
                absWeightDifference <- abs(adj1 - adj2)
                # matrix of sum of weights for each edge
                wightSum <- adj1 + adj2
                # Convert zeros to -1 to avoid division by zero in the next step
                wightSum[wightSum == 0] <- -1
                # matrix of normalised differences 
                normWeightDifference <- absWeightDifference / wightSum
                # 
                normWeightDifferenceSum <- sum(normWeightDifference)
                myDissim <- normWeightDifferenceSum / n_unionOfEdges
        }
        myDissim
}

# dissimilarity function
# IN 
#       adj1: numeric matrix [0,1], adjacency matrix for the layer 1 
#       adj2: numeric matrix [0,1], adjacency matrix for the layer 2
# OUT
#       dissimCoefficient: numeric [0,Inf), the dissimilarity coefficient for the two layers. The identical matrices will have a coefficient of 0 
d_deltaDegree <- function(adj1, adj2){
        N <- dim(adj1)[1]
        myArray <- array(c(adj1
                           ,adj2)
                         ,dim = c(N,N,2)
        )
        deltaDegrees <- calc_degrees(myArray)
        # calculate distance
        myDist <- dist(x = t(matrix(deltaDegrees[,3,], nrow = N))
                       ,method = "euclidean"
        )[1]
        myDist
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
                                dissimValue <- d(weight_array1[,,i]
                                                 ,weight_array2[,,j]
                                                 #,...
                                )
                                dissim_matrix[i,j] <- dissimValue
                                dissim_matrix[j,i] <- dissimValue
                        }
                }
        } else {
                for (i in 1:N){
                        for(j in 1:M){
                                dissimValue <- d(weight_array1[,,i]
                                                 ,weight_array2[,,j]
                                                 #,...
                                )
                                dissim_matrix[i,j] <- dissimValue
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
# deltaDegreeVector
# nodes
# palette
# OUT
delta_degree2color <- function(deltaDegreeVector
                               ,nodes
                               ,node_color_palette
){
        fine <- round(max(abs(deltaDegreeVector)) 
                      ,2) * 100 * 2 + 2
        if(fine == 2) {
                colors <- rep("gray40"#FFFFFF"
                              ,length(deltaDegreeVector))
                names(colors) <- nodes
                colors
        } else {
                palette <- colorRampPalette(node_color_palette)
                sapply(deltaDegreeVector
                       ,function(i){
                               palette(fine)[ceiling(round(i,2)*100 
                                                     + ((fine - 1)/2)
                               )
                               ]
                       })       
        }
}

create_empty_adj_matrix <- function(nodes){
        adjacencyMatrix <- matrix(0
                                  ,nrow=length(nodes)
                                  ,ncol=length(nodes)
        )
        dimnames(adjacencyMatrix) <- list(nodes
                                          ,nodes)
        adjacencyMatrix
}

# define a function to plot a community
# function to make a plot of communication graph for a given ligand-receptor pair
# IN:
#       LRP: string                             name of ligand-receptor pair for which the plot should be constructed. Note the the LRP name should be written in the same from as in ligandReceptorPairDF$pair. Alternatively, it can be a name of a cluster or a pattern
#       weight_array: 3D numeric array or 2D matrix  array of weighted adjacency matrices or one adjacency matrix with dimensions [number of nodes, number of nodes, number of ligand-receptor pairs]
#               First dimension: sending nodes (note: this dimension has all possible nodes, even if some of them are silent for a particular ligand-receptor pair)
#               Second dimension: receiving nodes (note: this dimension has all possible nodes, even if some of them are silent for a particular ligand-receptor pair)
#               Third dimension: ligand-receptor pairs
#               Note that the weight_array should contain dimnames: dimnames = list(nodes, nodes, ligand-receptor pairs)
#       ligandReceptorPairDF: string dataframe  data frame with columns "pair", "ligand", "ligand_complex_composition", "receptor", "receptor_complex_composition"
#               "pair" contains values in a form "ligand:receptor", i.e. ligand being at the first place, receptor being at the second place, e.g. "TNFSF13:TNFRSF17"
#               "ligand" contains ligand names, e.g. "TNFSF13"
#               "ligand_complex_composition" if ligand is a complex (e.g. "aXb2_complex"), contains genes in the ligand complex separated with a comma, e.g. "ITGAX,ITGB2", else contains ""
#               "receptor" contains receptor names, e.g. "TNFRSF17"
#               "receptor_complex_composition" if receptor is a complex (e.g. "NKG2D_II_receptor"), contains genes in the receptor complex separated with a comma, e.g. "KLRK1,HCST", else contains ""
#       nodes: string vector                    a vector with all cell types in the data
#       isPattern: logical                      should a graph of a pattern (and not a particular ligand-receptor pair) be plotted. Default value = F
#       isCluster: logical                      should a graph of a cluster (and not a particular ligand-receptor pair) be plotted. Default value = F
#       title: string                           title of the plot. Default value: LRP
#       subtitle: string                        subtitle of the plot. Default value: ""
#       node_color_palette: string vector       vector of colours for nodes. Default values: c("blue", "gray40", "red"). The colours will be used to construct a colour gradient for the node delta degree: blue representing a receiving node (negative delta degree), grey representing a neutral node (zero delta degree), red representing a sending node (positive delta degree)
#       node.label.cex: numeric                 size of node labels. Default value: 1
#       vertex.shape: string                    shape of node. Default value: "none". See same argument of plot.igraph()
#       vertex.size: integer                    size of node. Default value: 15. See same argument of plot.igraph()
#       edge.arrow.size: numeric                size of arrow edge. Default value: 0.5. See same argument of plot.igraph()
#       edge.width.scale: numeric               scaling factor for edge width. Default value: 2.5. See same argument of plot.igraph()
#       legend.size: numeric                    size of the legend text. Default value: 0.75
#       legend.gradient.x.coord: numeric vector x coordinates of the label position. Default value: c(1.15,1.25,1.25,1.15)
#       legend.gradient.y.coord: numeric vector y coordinates of the label position. Default value: c(-0.5,-0.5,-1,-1)
#       ...                             other plot.igraph parameters
# OUT:
#       graph plot
plot_communication_graph <- function(LRP
                                   ,weight_array
                                   ,ligandReceptorPairDF
                                   ,nodes
                                   ,isPattern = F
                                   ,isCluster = F
                                   ,title = LRP
                                   ,subtitle = ""
                                   ,node_color_palette = c("blue", "gray40", "red")
                                   ,node.label.cex = 1
                                   ,vertex.shape = "none"
                                   ,vertex.size = 15
                                   ,edge.arrow.size = 0.5
                                   ,edge.width.scale = 2.5
                                   ,legend.size = 0.75
                                   ,legend.gradient.x.coord = c(1.15,1.25,1.25,1.15)
                                   ,legend.gradient.y.coord = c(-0.5,-0.5,-1,-1)
                                   ,...
){
        # check that weight_array has dimnames
        if(length(dimnames(weight_array)) == 0) stop("ERROR: Dimension names of weight_array are not defined. Please define dimnames(weight_array)")
        
        # check whether the LRP is present in ligandReceptorPairDF. 
        # In some cases of between sample comparison it might be absent, 
        # so we will have to plot an empty plot
        present <- LRP %in% ligandReceptorPairDF$pair
        
        # define weight matrix
        if(isPattern | isCluster) {
                adjacencyMatrix <- weight_array
        } else if(present){ # LRP is present in the ligandReceptorPairDF
                
                # subset array
                adjacencyMatrix <- weight_array[,,LRP]
                # convert to matrix
                adjacencyMatrix <- array2matrix(adjacencyMatrix)
                
        } else { # LRP is NOT present in the ligandReceptorPairDF
                
                # create an empty weight matrix
                adjacencyMatrix <- create_empty_adj_matrix(nodes)
                # convert to matrix
                adjacencyMatrix <- array2matrix(adjacencyMatrix)
        }
        
        # calculate delta degree matrix
        degreeMatrix <- calc_degrees(weight_array = adjacencyMatrix)
        # convert to matrix
        degreeMatrix <- array2matrix(degreeMatrix)
        
        # create a directed weighted graph object
        layer_graph <- igraph::graph.adjacency(adjacencyMatrix
                                               ,mode = "directed"
                                               ,weighted = T
        )
        
        # add colour scheme to the graph object
        colors <- delta_degree2color(degreeMatrix[,3]
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
        leg.labels <- c(round(min(degreeMatrix[,3]),2)
                        ,round(max(degreeMatrix[,3]),1)
        )
        leg.title <- "delta\ndegree"
        
        # add size of node labels
        V(layer_graph)$label.cex = node.label.cex 
        
        V(layer_graph)$size = vertex.size
        
        # define edge.width
        if(length(E(layer_graph)$weight) != 0){
                edge.width <- (E(layer_graph)$weight/max(E(layer_graph)$weight))*edge.width.scale
        }else edge.width <- 0
        
        # set seed
        set.seed(1123)
        
        # create a plot
        plot.igraph(layer_graph
                    ,edge.width=edge.width
                    ,edge.arrow.size = edge.arrow.size
                    ,layout = layout.circle(layer_graph)
                    ,edge.curved=TRUE
                    ,vertex.label.color=V(layer_graph)$color
                    ,vertex.shape=vertex.shape
                    ,main = title
                    ,sub = subtitle
                    ,... # other plot.igraph parameters
        )
        #points for the gradient legend
        pnts = cbind(x = legend.gradient.x.coord
                     , y = legend.gradient.y.coord
        )
        #create the gradient legend
        legend.gradient(pnts = pnts
                        ,cols = colormap
                        ,limits = leg.labels
                        ,title = leg.title
                        ,cex = legend.size
        )
        
        
}

# function to plot a graph with ligands (red) and receptors (blue) in a particular cluster
# IN:
#       clusterOfInterest: numeric              number of the cluster, for which the ligands and receptors should be plotted
#       lrp_clusters: numeric vector            cluster assignment for each ligand-receptor pair
#       ligandReceptorPairDF: string dataframe  data frame with columns "pair", "ligand", "ligand_complex_composition", "receptor", "receptor_complex_composition"
#               "pair" contains values in a form "ligand:receptor", i.e. ligand being at the first place, receptor being at the second place, e.g. "TNFSF13:TNFRSF17"
#               "ligand" contains ligand names, e.g. "TNFSF13"
#               "ligand_complex_composition" if ligand is a complex (e.g. "aXb2_complex"), contains genes in the ligand complex separated with a comma, e.g. "ITGAX,ITGB2", else contains ""
#               "receptor" contains receptor names, e.g. "TNFRSF17"
#               "receptor_complex_composition" if receptor is a complex (e.g. "NKG2D_II_receptor"), contains genes in the receptor complex separated with a comma, e.g. "KLRK1,HCST", else contains ""
#       ligRecColor: string vector              colours for ligand and receptors: first colour for the ligand, second colour for the receptor. Default value: c("red", "blue") 
#       edge.arrow.size: numeric                size of arrow edge. Default value: 0.25. See same argument of plot.igraph()
#       node.label.cex: numeric                 size of node labels. Default value: 1
#       layout: func                            see layout from igraph package v0.2.1. Default value: layout_nicely 
#       legend.position: numeric vector         x, y coordinates of the legend position. Default value: c(-1, -1.1) 
#       legend.pt.cex: numeric                  expansion factor(s) for the points in the legend. See legend() from graphics v3.6.2. Default value: 2 
#       legend.cex: numeric                     character expansion factor in the legend. See legend() from graphics v3.6.2. Default value: 0.8
#       ...                                     other plot.igraph parameters
# OUT:
#       graph plot
plot_lig_rec <- function(clusterOfInterest
                       ,lrp_clusters
                       ,ligandReceptorPairDF
                       ,ligRecColor = c("red", "blue")
                       ,edge.arrow.size = 0.25
                       ,node.label.cex = 1
                       ,layout = layout_nicely
                       ,legend.position = c(-1, -1.1)
                       ,legend.pt.cex=2
                       ,legend.cex=0.8
                       ,...
){
        # make adjacency matrix: 
        # rows contain ligands
        # columns contain receptors
        # 1 if ligand and receptor is a pair, 0 if ligand and receptor is no pair
        lrp <- ligandReceptorPairDF[lrp_clusters == clusterOfInterest,]
        ligands <- unique(lrp$ligand)
        receptors <- unique(lrp$receptor)
        adjMat <- matrix(0
                         ,nrow = length(ligands) + length(receptors)
                         ,ncol = length(ligands) + length(receptors)
        )
        dimnames(adjMat) <- list(c(ligands,receptors)
                                 ,c(ligands,receptors)
        )
        for(idx in 1:nrow(lrp)){
                idxRow <- which(lrp$ligand[idx] == rownames(adjMat))
                idxCol <- which(lrp$receptor[idx] == colnames(adjMat))
                adjMat[idxRow, idxCol] <- 1
        }
        
        # create a directed weighted graph object
        lrp_graph <- igraph::graph.adjacency(adjMat
                                             ,mode = "directed"
        )
        # add color
        V(lrp_graph)$color <- {
                color <- c(rep(ligRecColor[1]
                               ,length(ligands))
                           ,rep(ligRecColor[2]
                                ,length(receptors)))
                names(color) <- c(ligands
                                  ,receptors)
                color <- color[colnames(adjMat)]
                color
        }
        # add label size
        V(lrp_graph)$label.cex = node.label.cex
        
        legend.lab <- c("ligand"
                        ,"receptor")
        legend.col <- ligRecColor
        title <- paste("cluster"
                       ,clusterOfInterest)
        
        edge.arrow.size <- edge.arrow.size
        
        # set seed
        set.seed(1123)
        
        # create a plot
        plot.igraph(lrp_graph
                    ,edge.arrow.size = edge.arrow.size
                    ,vertex.label.color=V(lrp_graph)$color
                    ,vertex.shape="none"
                    ,edge.curved=F
                    ,main = title
                    ,layout=layout
                    ,... # other plot.igraph parameters
        )
        legend(x=legend.position[1]
               ,y=legend.position[2]
               ,legend = legend.lab
               ,pch=21
               ,pt.bg=legend.col
               ,pt.cex=legend.pt.cex
               ,cex=legend.cex
               ,bty="n"
               ,ncol=1
        )
}


# function to plot heatmap of dissimilarity between two conditions
# IN:
# dissim_cond1_cond2: numeric matrix            pairwise dissimilarity between all ligand-receptor pairs in the two conditions (condition 1 in rows, condition 2 in columns)
# sortedLRP_DF: data.frame with columns:
# pair: string                          names of ligand-receptor pairs in the same form as they are in ligandReceptorPairDF$pair
# presence: string                      whether the ligand-receptor pair is present in both conditions ("shared") or only in one of them.
# dissimilarity: numeric                 dissimilarity value between the topology of ligand-receptor pair graph in two conditions. The smaller the dissimilarity value, the more similar is the graph topology between the two conditions. If a ligand-receptor pair is present only in one of the conditions, the dissimilarity is equal to 1.
# cond1_name: character string                  sample name for condition 1
# cond2_name: character string                  sample name for condition 2
# colors_lrp: string vector                     colours for ligand-receptor labels. Default: green for shared, black for unshared
# show_legend: logical                          parameter of the HeatmapAnnotation function. Default is T
# row_names_fontsize: numeric                   parameter of the Heatmap function. Default is 5
# colomn_names_fontsize: numeric                parameter of the Heatmap function. Default is 5
# show_column_names: logical                    parameter of the Heatmap function. Default is T
# show_row_names: logical                       parameter of the Heatmap function. Default  is T
# width: object of class "unit"                 parameter of the HeatmapAnnotation function. Default is unit(0.1, "cm")
# legend_direction: string                      parameter of the Heatmap function. Default is "horizontal"
# legend_width: object of class "unit"          parameter of the Heatmap function. Default is unit(5, "cm") 
# title_position: string                        parameter of the Heatmap function. Default is "lefttop"
# row_dend_side: string                         parameter of the Heatmap function. Default is "left"
# column_dend_side: string                      parameter of the Heatmap function. Default is "top"
# column_title_side: string                     parameter of the draw function. Default is "top"
# heatmap_legend_side: string                   parameter of the draw function. Default is "bottom" 
# ...                                           other parameters of the Heatmap function.
# OUT:
# two heatmaps:
# - one with ligand-receptor pairs clustered by their pairwise dissimilarity
# - one with ligand-receptor pairs sorted by decreasing dissimilarity of shared ligand-receptor pairs
plot_dissimilarity_heatmaps <- function(dissim_cond1_cond2
                                        ,sortedLRP_DF
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
        cond1_LRPcolors <- sapply(rownames(dissim_cond1_cond2)
                                  ,function(i){
                                          if(as.character(sortedLRP_DF[i,]$presence) == "shared") colors_lrp[1]
                                          else colors_lrp[2]
                                  })
        # colour vector for ligand-receptor pairs in condition 2
        cond2_LRPcolors <- sapply(colnames(dissim_cond1_cond2)
                                  ,function(i){
                                          if(as.character(sortedLRP_DF[i,]$presence) == "shared") colors_lrp[1]
                                          else colors_lrp[2]
                                  })
        # row annotations
        haRow = HeatmapAnnotation(df = data.frame(
                ligand_receptor_pair = c("shared"
                                         ,"not shared"))
                , col = list(ligand_receptor_pair = c("shared" = colors_lrp[1]
                                                      ,"not shared" = colors_lrp[2]))
                , which = "row"
                , width = width
                , show_legend = show_legend
        )
        h_clustered <- Heatmap(dissim_cond1_cond2
                               ,cluster_rows = T
                               ,row_names_gp = gpar(fontsize = row_names_fontsize
                                                    ,col = cond1_LRPcolors
                               )
                               ,cluster_columns = T
                               ,column_names_gp = gpar(fontsize = colomn_names_fontsize
                                                       ,col = cond2_LRPcolors
                               )
                               ,name = "dissimilarity coefficient"
                               ,show_column_names = show_column_names
                               ,show_row_names = show_row_names
                               ,heatmap_legend_param = list(legend_direction = legend_direction
                                                            ,legend_width = legend_width
                                                            ,title_position = title_position
                               )
                               ,row_title = cond1_name
                               ,row_dend_side = row_dend_side
                               ,column_title = cond2_name
                               ,column_dend_side = column_dend_side
                               ,...
        )
        heatmap.comparative.clustered <- draw(h_clustered + haRow
                                              ,column_title = "Comparative analysis: clustered"
                                              ,column_title_side = column_title_side
                                              ,heatmap_legend_side = heatmap_legend_side
        )
        
        h_sorted <- Heatmap(dissim_cond1_cond2
                            ,cluster_rows = F
                            ,row_order = {
                                    sub_sortedList <- sortedLRP_DF[rownames(sortedLRP_DF) %in% rownames(dissim_cond1_cond2),]
                                    rownames(sub_sortedList[order(
                                            sub_sortedList$presence),])
                            }
                            ,row_names_gp = gpar(fontsize = 5
                                                 ,col = cond1_LRPcolors
                            )
                            ,cluster_columns = F
                            ,column_order = {
                                    sub_sortedList <- sortedLRP_DF[rownames(sortedLRP_DF) %in% colnames(dissim_cond1_cond2),]
                                    rownames(sub_sortedList[order(
                                            sub_sortedList$presence),])
                            }
                            ,column_names_gp = gpar(fontsize = 5
                                                    ,col = cond2_LRPcolors
                            )
                            ,name = "dissimilarity coefficient"
                            , show_column_names = show_column_names
                            , show_row_names = show_row_names
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
        heatmap.comparative.sorted <- draw(h_sorted + haRow
                                           ,column_title = "Comparative analysis: sorted"
                                           ,column_title_side = column_title_side
                                           ,heatmap_legend_side = heatmap_legend_side
        )
}

# function to plot heatmap for clustered ligand-receptor pairs
# IN:
#       dissim_matrix: numeric matrix           pairwise dissimilarity matrix. Note that rownames and colnames should be defined
#       lrp_clusters: numeric vector  cluster assignment for each ligand-receptor pair. Note that the the vector should be named and ordered in the same way as the rownames and colnames of the dissim_matrix
#       cluster_colors: string vector           vector with colours for each cluster. Default value: ""
#       unclustered_LRP_color: string           colour for ligand-receptor pairs that remained unclustered (cluster assignment 0). Default value: "black"
#       show_colomn_names: logical              argument of the Heatmap function. Default value: F
#       show_row_names: logical                 argument of the Heatmap function. Default value: F
#       legend_direction: string                argument of the Heatmap function. Default value: "horizontal"
#       legend_width: object of class "unit"    argument of the Heatmap function. Default is: unit(5, "cm")
#       title_position: string                  argument of the Heatmap function. Default value: "lefttop"
#       column_title: string                    title of the heatmap plot. Argument of the draw function. Default value: "Clustering of ligand-receptor pairs"
#       column_title_side: string               position of the title of the heatmap plot. Argument of the draw function. Default value: "top"
#       heatmap_legend_side: string             position of the legend of the heatmap plot. Argument of the draw function. Default value: "bottom"
#       ...                                     any arguments "Heatmap" function takes
# OUT:
#       heatmap with ligand-receptor pairs ordered by their cluster assignment
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
        haCol <- HeatmapAnnotation(df = data.frame(cluster = as.vector(lrp_clusters))
                                   ,col = list(cluster = cluster_colors)
                                   ,annotation_legend_param = annotation_legend_param
        )
        
        # heatmap ordered by clusters
        h1 <- Heatmap(dissim_matrix
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
                                                   ,title_position = title_position#"leftcenter" #"topcenter" # "topleft" #"lefttop"
                      )
                      ,...
        )
        
        heatmap.clusters <- draw(h1
                                 ,column_title = column_title
                                 ,column_title_side = column_title_side
                                 ,heatmap_legend_side = heatmap_legend_side
        )
        # print(heatmap.clusters)
}

# function to plot UMAP of ligand-receptor pairs
# IN:
#       ligandReceptorPairDF: string dataframe  data frame with columns "pair", "ligand", "ligand_complex_composition", "receptor", "receptor_complex_composition"
#       dissim_matrix: numeric matrix           pairwise dissimilarity matrix. Note that rownames and colnames should be defined
#       lrp_clusters: numeric vector            cluster assignment for each ligand-recetor pair. Note that the the vector should be named and ordered in the same way as the rownames and colnames of the dissim_matrix
#       seed                                    random seed for UMAP. Default value: 100
#       cluster_colors: string vector           vector with colours for each cluster. Default value: ""
#       unclustered_LRP_color: string           colour for ligand-receptor pairs that remained unclustered (cluster assignment 0). Default value: "black"
#       cluster_shape: numeric vector           vector with shapes for each cluster. Default value: ""
#       point_size: numeric                     size of points. Default value: 2
#       title: string                           title of the heatmap plot. Argument of the draw function. Default value: "umap of ligand-receptor pairs"
#       legend.position: string                 argument of the theme function. Default value: "bottom" guide_legend
#       legend_direction: string                argument of the guide_legend function. Default value: "horizontal"
#       legend_title.position: string           argument of the guide_legend function. Default value: "top"
#       legend_label.position: string           argument of the guide_legend function. Default value: "bottom"
#       legend_byrow: logical                   argument of the guide_legend function. Default value: TRUE
# OUT:
#       UMAP of ligand-receptor pairs coloured and shaped by cluster
plot_cluster_UMAP <- function(ligandReceptorPairDF
                              ,dissim_matrix
                              ,lrp_clusters
                              ,seed = 100
                              ,cluster_colors = ""
                              ,unclustered_LRP_color = "black"
                              ,cluster_shape = ""
                              ,point_size = 2
                              ,title = "umap of ligand-receptor pairs"
                              ,legend.position="bottom"
                              ,legend_direction = "horizontal"
                              ,legend_title.position = "top"
                              ,legend_label.position= "bottom"
                              ,legend_byrow=TRUE
){
        
        # set seed
        set.seed(seed)
        
        # define n_neighbors parameter
        n_neighbors <- umap.defaults$n_neighbors
        if(n_neighbors >= nrow(ligandReceptorPairDF)) n_neighbors <- round(nrow(ligandReceptorPairDF) /2)
        
        # check wether there are enough ligand-receptor pairs to plot a UMAP
        if(n_neighbors > 1){
                custom.settings <-  umap.defaults
                custom.settings$input <- "dist"
                custom.settings$n_neighbors <- n_neighbors
                
                # calculate the UMAP
                my.umap <- umap(dissim_matrix
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
                umap.clusters <- ggplot(dfForUmap
                                        , aes(x=x
                                              ,y=y
                                              ,color=cluster
                                              ,shape = cluster
                                        )
                )+
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
                        theme(legend.position=legend.position) +
                        guides(fill=guide_legend(nrow=2
                                                 ,direction = legend_direction
                                                 ,title.position = legend_title.position
                                                 ,label.position= legend_label.position
                                                 ,byrow=legend_byrow
                        )
                        )
                print(umap.clusters)
        } else print("WARNING: too few ligand-receptor pairs. UMAP can not be plotted.")
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

# function to make a weight matrix for pattern of interest
# IN:
# communicatingNodes    character string vector in form of c("cellType1_to_celltype2", "cell_type1_to_celltype3",...), e.g. c("HSC_to_Tcells", "HSC_to_bcells",...)
# nodes                 character string vector of length n with all cell types in the data
# OUT:
# weight matrix         numeric n x n matrix [0,1]. Rows are nodes, columns are nodes. Rows are regarded as sending cell type, columns are regarded as receiving cell types. 1 defines communication, 0 defines no communication
make_pattern_matrix <- function(communicatingNodes
                              ,nodes
){
        # make a dummy matrix
        pattern_w_matrix <- matrix(0
                                   ,nrow = length(nodes)
                                   ,ncol = length(nodes)
        )
        dimnames(pattern_w_matrix) <- list(nodes,nodes)
        
        # define sending nodes
        sendingNodes <- gsub( "_.*$", "", communicatingNodes)
        print("sending nodes are:")
        print(sendingNodes)
        
        # check correctness of the node names
        if(!all(sendingNodes %in% nodes)) stop("ERROR: cell type names in 'communicatingNodes' argument are not correct. Please check spelling.")
        
        # define receiving nodes
        receivingNodes <- gsub('.*to_', "", communicatingNodes)
        print("recieving nodes are:")
        print(receivingNodes)
        
        
        # check correctness of the node names
        if(!all(receivingNodes %in% nodes)) stop("ERROR: cell type names in 'communicatingNodes' argument are not correct. Please check spelling.")
        
        # fill in 1 to communication nodes
        for(i in 1:length(communicatingNodes)){
                #print(i)
                pattern_w_matrix[sendingNodes[i],receivingNodes[i]] <- 1
        }
        
        # return
        return(pattern_w_matrix)
}

# converts CellPhoneDB output matrix (significant_means.txt) into an array of weighted adjacency matrices
# IN:
#       CellPhoneDB_output      significant_means.txt file in the output of CellPhoneDB
#       complex_input           complex_input.csv file. Download from CellPhoneDB
#       gene_input              gene_input.csv file. Download from CellPhoneDB
# OUT:
#       list of:
#               weight_array: 3D numeric array                       array of weighted adjacency matrices with dimensions [number of nodes, number of nodes, number of ligand-receptor pairs]. Note that the function filters out empty weight arrays and the corresponding ligand-receptor pairs
#                       First dimension: sending nodes (note: this dimension has all possible nodes, even if some of them are silent for a particular ligand-receptor pair)
#                       Second dimension: receiving nodes (note: this dimension has all possible nodes, even if some of them are silent for a particular ligand-receptor pair)
#                       Third dimension: ligand-receptor pairs
#                       Note that the weight_array should contain dimnames: dimnames = list(nodes, nodes, ligand-receptor pairs)
#                       
#               ligandReceptorPairDF: string dataframe          data frame with columns "pair", "ligand", "ligand_complex_composition", "receptor", "receptor_complex_composition"
#                        "pair" contains values in a form "ligand:receptor", i.e. ligand being at the first place, receptor being at the second place, e.g. "TNFSF13:TNFRSF17"
#                        "ligand" contains ligand names, e.g. "TNFSF13"
#                        "ligand_complex_composition" if ligand is a complex (e.g. "aXb2_complex"), contains genes in the ligand complex separated with a comma, e.g. "ITGAX,ITGB2", else contains ""
#                        "receptor" contains receptor names, e.g. "TNFRSF17"
#                        "receptor_complex_composition" if receptor is a complex (e.g. "NKG2D_II_receptor"), contains genes in the receptor complex separated with a comma, e.g. "KLRK1,HCST", else contains ""
#               nodes: character string vector                            a vector with all cell types in the data
convert_CellPhoneDB_output <- function(CellPhoneDB_output
                                       ,complex_input
                                       ,gene_input
){
        # make ligand-receptor pair list
        # ligandReceptorPairDF is data frame with columns "pair", "ligand", "ligand_complex_composition", "receptor", "receptor_complex_composition"
        # "pair" contains values in a form "ligand:receptor", i.e. ligand being at the first place, receptor being at the second place, e.g. "TNFSF13:TNFRSF17"
        # "ligand" contains ligand names, e.g. "TNFSF13"
        # "ligand_complex_composition" if ligand is a complex (e.g. "aXb2_complex"), contains genes in the ligand complex separated with a comma, e.g. "ITGAX,ITGB2", else contains ""
        # "receptor" contains receptor names, e.g. "TNFRSF17"
        # "receptor_complex_composition" if receptor is a complex (e.g. "NKG2D_II_receptor"), contains genes in the receptor complex separated with a comma, e.g. "KLRK1,HCST", else contains ""
        
        # ligands and receptors: each pair will have a form "ligand:receptor"
        
        # indices of receptors
        idxReceptor <- as.matrix(data.frame(a = CellPhoneDB_output$receptor_a
                                            ,b = CellPhoneDB_output$receptor_b))
        rownames(idxReceptor) <- CellPhoneDB_output$interacting_pair
        
        # which pairs are directed
        idxDirected <- rowSums(idxReceptor) == 1
        
        # make empty ligandReceptorPairDF data.frame
        ligandReceptorPairDF <- as.data.frame(matrix(NA
                                                     ,nrow = nrow(CellPhoneDB_output)
                                                     ,ncol = 5
        ))
        colnames(ligandReceptorPairDF) <- c("pair"
                                            ,"ligand"
                                            ,"ligand_complex_composition"
                                            ,"receptor"
                                            ,"receptor_complex_composition"
        )
        
        # put in ligand and receptors
        splitLRP <- {
                # check if any ligands or receptor complexes are written with "_"
                # we will split pairs by "_", so any "_" in a complex name will interfere and we want to replace it with " "
                whichToReplace <- function(idxVector # all positions of "_" in the pair name
                                           ,partner # which partner is a complex
                                           ,howMany # how many "_" does this complex contain (i.e. how many do we want to replace with " ")
                ){
                        if(partner == "partner_a") {
                                idxVector <- idxVector[1:howMany] # take first "howMany" values
                        } else {
                                idxVector <- idxVector[(length(idxVector) - howMany + 1
                                ):length(idxVector)] # take last "howMany
                        }
                } # the function returns a vector with exact positions of "_" in the pair name that should be replaced with " "
                
                # replace "_" with " " in the complex names in the CellPhoneDB_output$interacting_pair
                for(partner in c("partner_a"
                                 ,"partner_b"
                )
                ){
                        # which complexes in partern_a or in pertner_b colomn contain "_"
                        idxComplex <- grep("_"
                                           ,unlist(CellPhoneDB_output[,partner])
                        )
                        # replace "_" for space in each position
                        if(length(idxComplex) != 0){
                                for(lrp in idxComplex) {
                                        # how many "_" does a particular partner_a or partner_b complex contain
                                        howMany <- length(unlist(gregexpr(pattern = "_"
                                                                          ,unlist(CellPhoneDB_output[lrp,partner]
                                                                          ))
                                        )
                                        )
                                        # find all positions of "_" in the interacting_pair
                                        idxChar <- unlist(gregexpr(pattern = "_"
                                                                   ,CellPhoneDB_output$interacting_pair[lrp]
                                        )
                                        )
                                        # take only those that correspond to the prtner_a or prtner_b
                                        idxChar <- whichToReplace(idxVector = idxChar
                                                                  ,partner = partner
                                                                  ,howMany = howMany)
                                        # replace them for " "
                                        for(idx in idxChar){
                                                substr(CellPhoneDB_output$interacting_pair[lrp]
                                                       ,idx
                                                       ,idx) <- " "
                                        }
                                }
                        }
                }
                
                # update rownames
                rownames(CellPhoneDB_output) <- CellPhoneDB_output$interacting_pair
                rownames(idxReceptor) <- CellPhoneDB_output$interacting_pair
                
                # split by "_"
                splitLRP <- strsplit(CellPhoneDB_output$interacting_pair
                                     ,"_")
                
                # check splitting correctness
                nrPartners <- sapply(1:nrow(CellPhoneDB_output)
                                     ,function(i) length(splitLRP[[i]])
                )
                if(any(nrPartners!=2)) stop("ERROR: Not all pairs contain 2 partners.")
                # translate to data frame
                splitLRP <- matrix(unlist(splitLRP)
                                   ,ncol = 2
                                   ,byrow = T)
        }
        ligandReceptorPairDF$ligand <- as.vector(splitLRP[,1])
        ligandReceptorPairDF$receptor <- as.vector(splitLRP[,2])
        
        # check if any ligand and receptor pairs are swapped
        idxSwapped <- !idxReceptor[,2] & idxDirected
        ligands <- ligandReceptorPairDF$receptor[idxSwapped]
        ligandReceptorPairDF$receptor[idxSwapped] <- ligandReceptorPairDF$ligand[idxSwapped]
        ligandReceptorPairDF$ligand[idxSwapped] <- ligands
        
        # put in pairs
        ligandReceptorPairDF$pair <- paste0(ligandReceptorPairDF$ligand
                                            ,":"
                                            ,ligandReceptorPairDF$receptor)
        CellPhoneDB_output$sorted_pair <- ligandReceptorPairDF$pair
        
        # put in complex composition
        ligandReceptorPairDF$ligand_complex_composition <- suppressWarnings( # suppresses warning: In `[.data.frame`(..., drop = FALSE) : 'drop' argument will be ignored
                sapply(ligandReceptorPairDF$ligand
                       ,function(g) {
                               genes <- complex2genes(g
                                                      ,complex_input = complex_input
                                                      ,gene_input = gene_input)
                               paste(unlist(genes)
                                     ,collapse = ",")
                       }
                ))
        ligandReceptorPairDF$receptor_complex_composition <- suppressWarnings( # suppresses warning: In `[.data.frame`(..., drop = FALSE) : 'drop' argument will be ignored
                sapply(ligandReceptorPairDF$receptor
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
        rownames(ligandReceptorPairDF) <- CellPhoneDB_output$interacting_pair
        
        # check NAs
        check_NAs(ligandReceptorPairDF)
        
        # make a key of nodes (in order to convert a vector of pairwise nodes to an adjacency matrix)
        nodeCombinations <- colnames(CellPhoneDB_output)[grepl("[|]"
                                                               ,colnames(CellPhoneDB_output))]
        node_key <- matrix(unlist(strsplit(nodeCombinations
                                           ,"[|]")
        )
        , ncol = 2
        , byrow = T)
        dimnames(node_key) <- list(nodeCombinations
                                   ,c("a", "b"))
        
        nodes <- unique(c(node_key[,"a"]
                          ,node_key[,"b"])
        )
        # make weight array
        # make a dummy array of weights
        weight_array <- array(0
                         ,dim = c(length(nodes) # sending nodes
                                  ,length(nodes) # recieving nodes
                                  ,nrow(ligandReceptorPairDF) # ligand-receptor pairs
                         )
        )
        dimnames(weight_array) <- list(unique(c(node_key[,"a"]
                                           ,node_key[,"b"])
        )
        ,unique(c(node_key[,"a"]
                  ,node_key[,"b"])
        )
        ,ligandReceptorPairDF$pair
        )
        # populate array of weights
        idxSwapped <- which(idxSwapped)
        for(layer in CellPhoneDB_output$interacting_pair){
                for(nodePair in nodeCombinations){
                        if(!is.na(CellPhoneDB_output[layer
                                                     ,nodePair])) {
                                idxLayerCPDB <- which(rownames(ligandReceptorPairDF) == layer)
                                # control swapped pairs
                                if(idxLayerCPDB %in% idxSwapped) {
                                        ligand <- 2
                                        receptor <- 1
                                } else {
                                        ligand <- 1
                                        receptor <- 2
                                }
                                weight_array[node_key[nodePair,][ligand] # partner a (ligand)
                                        ,node_key[nodePair,][receptor] # partner b (receptor)
                                        ,as.character(ligandReceptorPairDF[layer,"pair"])
                                        ] <- as.numeric(CellPhoneDB_output[idxLayerCPDB 
                                                                           ,nodePair])
                        }
                }
        }
        
        # filter empty weight matrices
        filteredData <- filter_empty_layers(weight_array = weight_array
                                          ,ligandReceptorPairDF = ligandReceptorPairDF
        )
        ligandReceptorPairDF <- filteredData$ligandReceptorPairDF
        weight_array_filtered <- filteredData$weight_array
        rm(filteredData)
        
        
        # print a warning of number of ligand-receptor pairs is lower than 10
        if(nrow(ligandReceptorPairDF) < 10) print("WARNING: there are very few ligand-receptor pairs!")
        
        # return
        result <- list(weight_array = weight_array_filtered
                       ,ligandReceptorPairDF = ligandReceptorPairDF
                       ,nodes = nodes
        )
}

# function to cluster ligand-receptor pairs
# IN:
#       weight_array: numeric array                     array of weighted adjacency matrices with dimensions [number of nodes, number of nodes, number of ligand-receptor pairs]
#       ligandReceptorPairDF: character dataframe          data frame with columns "pair", "ligand", "ligand_complex_composition", "receptor", "receptor_complex_composition"
#                                                       "pair" contains values in a form "ligand:receptor", i.e. ligand being at the first place, receptor being at the second place, e.g. "TNFSF13:TNFRSF17"
#                                                       "ligand" contains ligand names, e.g. "TNFSF13"
#                                                       "ligand_complex_composition" if ligand is a complex (e.g. "aXb2_complex"), contains genes in the ligand complex separated with a comma, e.g. "ITGAX,ITGB2", else contains ""
#                                                       "receptor" contains receptor names, e.g. "TNFRSF17"
#                                                       "receptor_complex_composition" if receptor is a complex (e.g. "NKG2D_II_receptor"), contains genes in the receptor complex separated with a comma, e.g. "KLRK1,HCST", else contains ""
#       nodes: character string vector                            a vector with all cell types in the data
#       dissimilarity: func                             dissimilarity function. Default value is d_normWeightDiff
#       minClusterSize: numeric                         minimum number of ligand-receptor pairs per cluster. Default value: 6
# OUT:
#       a list of:
#               dissim_matrix: numeric matrix           pairwise dissimilarity between all ligand-receptor pairs
#               clusters: numeric vector                cluster assignment for each ligand-receptor pair
#               weight_array_byCluster: numeric array   array of weighted adjacency matrices with dimensions [number of nodes, number of nodes, number of ligand-receptor clusters]. Each edge weight in a cluster is calculated as the mean of the weights of this edge in all ligand-receptor pairs that belong to this cluster.
lrp_clustering <- function(weight_array
                           ,ligandReceptorPairDF
                           ,nodes
                           ,dissimilarity = d_normWeightDiff
                           ,minClusterSize = 6
                           
){
        # calculate degree array
        degree_array <- calc_degrees(weight_array = weight_array)
        hist(degree_array[,3,], main = "Histogram of degree delta")
        
        ## run analysis ####
        # check if there are more than 1 ligand receptor pairs
        if(nrow(ligandReceptorPairDF) <= 1) {
                dissim_matrix <- NA
                my.clusters <- 0
                weight_array_byCluster <- NA
                degreeArrayForClusters <- NA
        } else {
                # dissimilatriy ####
                dissim_matrix <- calc_dissim_matrix(weight_array1 = weight_array
                                                    ,weight_array2 = weight_array
                                                    ,dissimilarity = dissimilarity
                )
                hist(dissim_matrix)
                
                # cluster layers ####
                # convert dissimilarity matrix to distance matrix
                my.dist <- as.dist(dissim_matrix)
                # Make hierarchical clustering
                my.tree <- hclust(my.dist
                                  , method= "average" 
                )
                
                # Cut the treee to identify clusters
                my.clusters <- cutreeHybrid(my.tree
                                            ,distM=dissim_matrix
                                            ,deepSplit = 0
                                            ,minClusterSize = minClusterSize 
                )$label
                names(my.clusters) <- dimnames(dissim_matrix)[[1]]
                
                # Any unassigned cells?
                if(sum(my.clusters == 0) != 0) print("Warning: some graphs are not assigned to any cluster")
                # print out number of clusters
                print(paste("We have", max(as.numeric(my.clusters)), "clusters"))
                
                # calculate adjacency matrices for each community (using mean weight for each edge)
                if(any(my.clusters != 0)){
                        weight_array_byCluster <- array(NA
                                                   ,dim = c(length(nodes)
                                                            ,length(nodes)
                                                            ,max(my.clusters)
                                                   )
                        )
                        for(communityNr in 1:max(my.clusters)){
                                # calculate indices for the community
                                myCommunityIdx <- which(as.numeric(my.clusters) == communityNr)
                                # populate matrix
                                for(i in 1:length(nodes)){
                                        for(j in 1:length(nodes)){
                                                weight_array_byCluster[i,j,communityNr] <- mean(weight_array[i,j,myCommunityIdx]) 
                                        }
                                }   
                        }
                        dimnames(weight_array_byCluster) <- list(nodes
                                                            ,nodes
                                                            ,paste("cluster"
                                                                   ,1:max(my.clusters))
                        )
                        check_NAs(weight_array_byCluster)
                        
                        # calcualte degree arrays for each community (put it as well as an array)
                        degreeArrayForClusters <- calc_degrees(weight_array = weight_array_byCluster)
                        
                } else {
                        weight_array_byCluster <- NA
                        degreeArrayForClusters <- NA
                }
                
                # end of definitions and data processing
                
        }
        
        
        ## produce report ####
        # make object ####
        result <- list(dissim_matrix = dissim_matrix
                       ,clusters = my.clusters
                       ,weight_array_byCluster = weight_array_byCluster
        )
        
        
        ## return result ####
        return(result)
}

# function to perform pattern search
# IN:
#       patternAdjMatrix: numeric matrix        an adjacency matrix for the pattern of interest
#       weight_array: numeric array             array of weighted adjacency matrices with dimensions [number of nodes, number of nodes, number of ligand-receptor pairs]
#       ligandReceptorPairDF: character dataframe  data frame with columns "pair", "ligand", "ligand_complex_composition", "receptor", "receptor_complex_composition"
#       nodes: character string vector                    a vector with all cell types in the data
#       dissimilarity: function                 dissimilarity measure. Default value: d_normWeightDiff
# OUT:
#       dataframe:                              The data frame is sorted by increasing dissimilarity (i.e. similar patterns a the top
#               pair: character string vector             ligand-receptor pair names
#               dissimilarity: numeric vector  dissimilarity of ligand-receptor pairs to the pattern of interest. Note here that for the calculation of the dissimilarity of a weight matrix to the pattern of interest, we don't take into consideration the actual weight of the edges in this weight matrix, but just check wether the edge is present (i.e. has a non-zero value) or absent (i.e. has a zero value).
pattern_search <- function(patternAdjMatrix # in rows are sending nodes, in columns are receiving nodes, 0 for no communication, 1 for communication
                          ,weight_array
                          ,ligandReceptorPairDF
                          ,nodes
                          ,dissimilarity = d_normWeightDiff
){
        # check whether patternAdjMatrix contains all nodes
        if((!identical(rownames(patternAdjMatrix)
                       ,nodes
        ) ) |
        (!identical(colnames(patternAdjMatrix)
                    ,nodes)
        )
        ) stop("ERROR: Row names or column names of the patternAdjMatrix are not identical to the nodes.")
        
        # convert patternAdjMatrix to array, We need it to pass two arrays to the dissimilarity function
        patternAdjArray <- array(patternAdjMatrix
                                 ,dim = list(dim(patternAdjMatrix)[[1]]
                                             ,dim(patternAdjMatrix)[[2]]
                                             ,1)
        )
        
        # binarize weight_array
        binarized_weight_array <- (weight_array != 0) *1
        
        # calculate dissimilarity to the each LRP
        dissim <- sapply(1:dim(binarized_weight_array)[[3]]
                         ,function(k){
                                 dissimilarity(patternAdjArray
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

# function to perform comparative analysis
# IN:
# cond1_weight_array: numeric array                     array of weighted adjacency matrices for condition 1. Note that the function filters out empty weight arrays and the corresponding ligand-receptor pairs
# cond2_weight_array: numeric array                     array of weighted adjacency matrices for condition 2. Note that the function filters out empty weight arrays and the corresponding ligand-receptor pairs
# cond1_ligandReceptorPairDF: character dataframe          data frame with columns "pair", "ligand", "ligand_complex_composition", "receptor", "receptor_complex_composition" for condition 1
# cond2_ligandReceptorPairDF: character dataframe          data frame with columns "pair", "ligand", "ligand_complex_composition", "receptor", "receptor_complex_composition" for condition 2
# cond1_nodes: character string vector                            a vector with all cell types in the data
# Please note that comparativeAnalysis expects nodes (i.e. cell types) to be the same in both conditions.
# cond2_nodes = cond1_nodes: character string vector              a vector with all cell types in the data
# Please note that comparativeAnalysis expects nodes (i.e. cell types) to be the same in both conditions.
# cond1_name: character string                          sample name for condition 1
# cond2_name: character string                          sample name for condition 2
# dissimilarity: func                                   dissimilarity function. Default value is d_normWeightDiff
# OUT:
# list of:
# sortedLRP_DF: data.frame with columns:
# 	pair: character string                          names of ligand-receptor pairs in the same form as they are in ligandReceptorPairDF$pair
# 	presence: character string                      whether the ligand-receptor pair is present in both conditions ("shared") or only in one of them.
# 	dissimilarity: numeric                 dissimilarity value between the topology of ligand-receptor pair graph in two conditions. The smaller the dissimilarity value, the more similar is the graph topology between the two conditions. If a ligand-receptor pair is present only in one of the conditions, the dissimilarity is equal to 1.
# dissim_cond1_cond2: numeric matrix            pairwise dissimilarity between all ligand-receptor pairs in the two conditions (condition 1 in rows, condition 2 in columns)

comparative_analysis <- function(cond1_weight_array
                                ,cond2_weight_array
                                ,cond1_ligandReceptorPairDF
                                ,cond2_ligandReceptorPairDF
                                ,cond1_nodes
                                ,cond2_nodes = cond1_nodes
                                ,cond1_name
                                ,cond2_name
                                ,dissimilarity = d_normWeightDiff 
){
        
        # check if the nodes in both conditions are the same
        if(!identical(cond1_nodes, cond2_nodes)) stop("ERROR: one of the conditions might contain cell populations that are missing in the other. Please make sure both conditions contain same cell populations.")
        
        # stratify LRPs by shared or not shared
        allLRPs <- union(cond1_ligandReceptorPairDF$pair
                         ,cond2_ligandReceptorPairDF$pair
        )
        
        sharedLRP <- intersect(cond1_ligandReceptorPairDF$pair
                               ,cond2_ligandReceptorPairDF$pair
        )
        
        onlyCond1LRPs <- setdiff(cond1_ligandReceptorPairDF$pair
                                 ,cond2_ligandReceptorPairDF$pair
        )
        
        onlyCond2LRPs <- setdiff(cond2_ligandReceptorPairDF$pair
                                 ,cond1_ligandReceptorPairDF$pair
        )
        
        # calculate dissimilarity
        dissim_cond1_cond2 <- calc_dissim_matrix(weight_array1 = cond1_weight_array
                                                 ,weight_array2 = cond2_weight_array
                                                 ,dissimilarity = dissimilarity
        )
        
        # make sorted LRP list
        sortedLRP_DF <- as.data.frame(t(cbind(sapply(allLRPs
                                                      ,function(lrp){
                                                              if((length(sharedLRP) != 0) & (lrp %in% sharedLRP)){
                                                                      c(lrp
                                                                        ,"shared"
                                                                        ,unlist(dissim_cond1_cond2[lrp,lrp]))
                                                              } else if((length(onlyCond1LRPs) != 0) & (lrp %in% onlyCond1LRPs)) {
                                                                      c(lrp
                                                                        ,paste("only"
                                                                               ,cond1_name)
                                                                        ,1)
                                                              } else if((length(onlyCond2LRPs) != 0) & (lrp %in% onlyCond2LRPs)) {
                                                                      c(lrp
                                                                        ,paste("only"
                                                                               ,cond2_name)
                                                                        ,1)
                                                              }
                                                      })
        )))
        colnames(sortedLRP_DF) <- c("pair"
                                     ,"presence"
                                     ,"dissimilarity")
        rownames(sortedLRP_DF) <- sortedLRP_DF$pair
        
        # sort by accending dissimilarity
        sortedLRP_DF$dissimilarity <- as.numeric(sortedLRP_DF$dissimilarity)
        sortedLRP_DF <- sortedLRP_DF[order(sortedLRP_DF$dissimilarity
                                             ,decreasing = T),]
        
        # sort by presence
        sortedLRP_DF$presence <- factor(sortedLRP_DF$presence
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
        sortedLRP_DF <- sortedLRP_DF[order(sortedLRP_DF$presence
                                             ,decreasing = F),]
        
        return(list(sortedLRP_DF = sortedLRP_DF
                    ,dissim_cond1_cond2 = dissim_cond1_cond2)
        )
        
}