#' @rdname make_pattern_matrix
#' @export
#'
#' @title
#' Makes a weight matrix for pattern of interest
#'
#' @description
#' \code{make_pattern_matrix} makes a weight matrix for pattern of interest.
#' The rows and columns of the matrix contain all nodes.
#' The rows are considered as sending nodes and the columns are considered as receiving nodes.
#' The value of 1 represents an edge (i.e. cell types are communicating) and 0 represents no edge (i.e. cell types are not communicating).
#'
#' @author
#' Maria Solovey \email{maria.solovey@helmholtz-muenchen.de}
#'
#' @param communicating_nodes Character string vector: in form of c("cellType1_to_celltype2", "cell_type1_to_celltype3", ...), e.g. c("HSC_to_Tcells", "HSC_to_bcells", ...).
#'
#' @param nodes Character string vector: of length n with all cell types in the data.
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
