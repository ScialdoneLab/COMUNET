#' CellPhoneDB output of significant interactions.
#'
#' The dataset contains significant interactions called by CellPhoneDB algorithm on the AML328-d29 sample dataset.
#' Original CellPhoneDB file name is significant_means.txt.
#'
#' @format A data frame with 315 rows and 208 variables:
#'
#' Rows contain ligand-receptor pairs.
#'
#' Columns contain same columns as the significant_means.txt file.
#'
"AML328_d29_CellPhoneDB_output"

# define input.path
input.path <- getwd()

options(stringsAsFactors = F)

# read in CellPhoneDB output
AML328_d29_CellPhoneDB_output <- read.csv(paste0(input.path ,"/AML328_d29_significant_means.txt")
                                         ,sep = "\t" ,check.names = F)
# delete duplicates
AML328_d29_CellPhoneDB_output <-AML328_d29_CellPhoneDB_output[!duplicated(AML328_d29_CellPhoneDB_output$interacting_pair),]
# add row names (interacting pairs)
rownames(AML328_d29_CellPhoneDB_output) <- AML328_d29_CellPhoneDB_output$interacting_pair
# transform 'receptor_a' colomn into boolean
AML328_d29_CellPhoneDB_output$receptor_a <- sapply(AML328_d29_CellPhoneDB_output$receptor_a
                                                  ,function(i){
                                                          if(i == "True") T else F }
)
# transform 'receptor_b' colomn into boolean
AML328_d29_CellPhoneDB_output$receptor_b <- sapply(AML328_d29_CellPhoneDB_output$receptor_b
                                                  ,function(i){
                                                          if(i == "True") T else F }
)
use_data(AML328_d29_CellPhoneDB_output
         ,overwrite = TRUE)
