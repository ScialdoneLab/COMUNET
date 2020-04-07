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
