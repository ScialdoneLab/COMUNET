#' @rdname make_pattern_matrix
#' @export
#'
#' @title
#' Makes a weight matrix for pattern of interest
#'
#' @description
#' Makes a weight matrix for pattern of interest.
#'
#' The rows and columns of the matrix contain all nodes.
#' The rows are considered as sending nodes and the columns are considered as receiving nodes.
#' The value of 1 represents an edge (i.e. cell types are communicating) and 0 represents no edge (i.e. cell types are not communicating).
#'
#' Please nota that for simplicity, we address all interacting partners (including non-directional partners such as adhesion molecules) as "ligand-receptor pairs".
#'
#' @author
#' Maria Solovey \email{maria.solovey@helmholtz-muenchen.de}
#'
#' @param communicating_nodes Character string vector: in form of c("cellType1_to_cellType2", "cellType1_to_cellType3", ...), e.g. c("HSC_to_Tcells", "HSC_to_Bcells", ...).
#'
#' @param nodes Character string vector: of length n with all cell types in the data.
#'
#' @return
#' \item{weight_matrix}{
#'
#'   Numeric \eqn{n} x \eqn{n} matrix with \eqn{n} being the number of nodes and the values of the matrix being 0 or 1.
#'   Rows are regarded as sending cell type, columns are regarded as receiving cell types. A value of 1 defines communication, a value of 0 defines no communication.
#' }
#'
#' @examples
#'
#' # load embryo_interactions
#' data("embryo_interactions")
#'
#' communicating_nodes <- c("exVE_to_EPI" ,"exVE_to_Mes","exVE_to_TE" ,"exVE_to_emVE" ,"exVE_to_exVE")
#'
#' test_pattern <- make_pattern_matrix(communicating_nodes = communicating_nodes,
#'                    nodes = embryo_interactions$nodes)
#' print(test_pattern)
#'
#'
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
        print("receiving nodes are:")
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
