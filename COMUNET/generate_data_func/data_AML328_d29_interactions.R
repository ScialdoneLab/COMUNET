#' COMUNET format of CellPhoneDB communication output.
#'
#' Dataset contains a weight array for all ligand-receptor interactions,
#' a data frame with all active ligand-receptor pairs,
#' and a vector of all cell types (nodes) present in the data.
#'
#' @format A data list of three variables
#' \describe{
#' \item{weight_array}{
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
"AML328_d29_interactions"

options(stringsAsFactors = F)

# load CellPhoneDB output for embryo data
data("AML328_d29_CellPhoneDB_output")

# load complex_input table
data("complex_input")

# load gene_input table
data("gene_input")


# transform CellPhoneDB output
AML328_d29_interactions <- convert_CellPhoneDB_output(CellPhoneDB_output = AML328_d29_CellPhoneDB_output
                                                     ,complex_input = complex_input
                                                     ,gene_input = gene_input)


use_data(AML328_d29_interactions
         ,overwrite = TRUE)
