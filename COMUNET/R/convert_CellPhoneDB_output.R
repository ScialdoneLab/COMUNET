#' @rdname convert_CellPhoneDB_output
#' @export
#'
#' @title
#' Converts CellPhoneDB output matrix into an array of weighted adjacency matrices
#'
#' @description
#' Converts CellPhoneDB output matrix (significant_means.txt) into an array of weighted adjacency matrices.
#'
#' The CellPhoneDB weight matrix (significant_means.txt) contains ligand-receptor pairs in the rows and all pairs of cell types in the columns.
#' COMUNET transforms it into a stack of weight matrices (i.e. weight array), one weight matrix per ligand-receptor pair.
#' If CellPhoneDB output contains \eqn{m} non-empty ligand-receptor pairs, then the number of matrices in the stack will be equal to \eqn{m}.
#' Each such weight matrix has all cell types in the rows and all cell types in the columns.
#' By this, if there are \eqn{n} cell types in the data, each weight matrix is an \eqn{n} x \eqn{n} matrix.
#'
#' The rows of each matrix represent sending cell types, i.e. cell types that express a ligand of a ligand-receptor pair,
#' or just partner A of an undirected interacting pair (e.g. pair of adhesion molecules).
#' The columns of each matrix represent receiving cell types, i.e. cell types that express a receptor of a ligand-receptor pair,
#' or partner B of an undirected interacting pair (e.g. pair of adhesion molecules).
#' The NA values of the CellPhoneDB output are substituted with 0.
#'
#' @author
#' Maria Solovey \email{maria.solovey@helmholtz-muenchen.de}
#'
#' @param       CellPhoneDB_output Character string: significant_means.txt file in the output of CellPhoneDB.
#'
#' @param       complex_input Character string: complex_input.csv file. Download from CellPhoneDB.
#'
#' @param       gene_input Character string: gene_input.csv file. Download from CellPhoneDB.
#'
#' @return A list of:
#'
#' \itemize{
#'   \item{weight_array}{
#'
#'   Numeric array (3D): array of weighted adjacency matrices with dimensions [number of nodes, number of nodes, number of ligand-receptor pairs].
#'   \itemize{
#'     \item First dimension: sending nodes (note: this dimension has all possible nodes, even if some of them are silent for a particular ligand-receptor pair).
#'     \item Second dimension: receiving nodes (note: this dimension has all possible nodes, even if some of them are silent for a particular ligand-receptor pair).
#'     \item Third dimension: ligand-receptor pairs.
#'   }
#'   Note that the function filters out empty weight arrays and the corresponding ligand-receptor pairs.
#'
#'   Note that the weight_array should contain dimnames: dimnames = list(nodes, nodes, ligand-receptor pairs).
#' }
#'   \item{ligand_receptor_pair_df}{
#'
#'   Character string data frame: data frame with columns "pair", "ligand", "ligand_complex_composition", "receptor", "receptor_complex_composition".
#'   \itemize{
#'     \item "pair" contains values in a form "ligand:receptor", i.e. ligand being at the first place, receptor being at the second place, e.g. "TNFSF13:TNFRSF17".
#'     \item "ligand" contains ligand names, e.g. "TNFSF13".
#'     \item "ligand" contains ligand names, e.g. "TNFSF13".
#'     \item "ligand_complex_composition" if ligand is a complex (e.g. "aXb2_complex"), contains genes in the ligand complex separated with a comma, e.g. "ITGAX,ITGB2", else contains empty string "".
#'     \item "receptor" contains receptor names, e.g. "TNFRSF17".
#'     \item "receptor_complex_composition" if receptor is a complex (e.g. "NKG2D_II_receptor"), contains genes in the receptor complex separated with a comma, e.g. "KLRK1,HCST", else contains empty string "".
#'
#'   }
#'   }
#'   \item{nodes}{
#'
#'   Character string vector: a vector with all cell types in the data.
#'   }
#' }
#'
#' @references
#' CellPhoneDB:
#'
#' \href{https://www.nature.com/articles/s41596-020-0292-x?proof=trueMay}{Efremova et al., \emph{Nature Protocols} 2019}
#'
#' \href{http://dx.doi.org/10.1038/s41586-018-0698-6}{Vento-Tormo et al., \emph{Nature} 2018}
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
                        # which complexes in partern_a or in partner_b colomn contain "_"
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
