#' @rdname pattern_search
#' @export
#'
#' @title
#' Performs pattern search
#'
#' @description
#' \code{pattern_search} performs a search search for ligand-receptor pairs with a pattern of communication similar to the pattern of interest.
#' A pattern of interest is defined with a binary adjacency matrix with values 0 or 1, the value 1 indicating the presence of an edge and 0 the absence.
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
#' @param       nodes character string vector: a vector with all cell types in the data.
#'
#' @param       dissimilarity function: dissimilarity measure. Default value: d_normWeightDiff.
#'
#' @return dataframe:                              The data frame is sorted by increasing dissimilarity (i.e. similar patterns a the top
#'               "pair" (character string vector): ligand-receptor pair names
#'               "dissimilarity" (numeric vector): dissimilarity of ligand-receptor pairs to the pattern of interest. Note here that for the calculation of the dissimilarity of a weight matrix to the pattern of interest, we don't take into consideration the actual weight of the edges in this weight matrix, but just check wether the edge is present (i.e. has a non-zero value) or absent (i.e. has a zero value).
pattern_search <- function(pattern_adj_matrix
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
