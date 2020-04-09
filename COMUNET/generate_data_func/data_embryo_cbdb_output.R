#' CellPhoneDB output of significant interactions.
#'
#' The dataset contains significant interactions called by CellPhoneDB algorithm on the embryo dataset.
#' Original CellPhoneDB file name is significant_means.txt.
#'
#' @format A data frame with 530 rows and 37 variables:
#'
#' Rows contain ligand-receptor pairs.
#'
#' Columns contain same columns as the significant_means.txt file.
#'
"embryo_CellPhoneDB_output"

# define input.path
input.path <- getwd()

options(stringsAsFactors = F)

# read in CellPhoneDB output
embryo_CellPhoneDB_output <- read.csv(paste0(input.path ,"/embryo_significant_means.txt")
                               ,sep = "\t" ,check.names = F)
# delete duplicates
embryo_CellPhoneDB_output <-embryo_CellPhoneDB_output[!duplicated(embryo_CellPhoneDB_output$interacting_pair),]
# add row names (interacting pairs)
rownames(embryo_CellPhoneDB_output) <- embryo_CellPhoneDB_output$interacting_pair
# transform 'receptor_a' colomn into boolean
embryo_CellPhoneDB_output$receptor_a <- sapply(embryo_CellPhoneDB_output$receptor_a
                                        ,function(i){
                                                if(i == "True") T else F }
                                        )
# transform 'receptor_b' colomn into boolean
embryo_CellPhoneDB_output$receptor_b <- sapply(embryo_CellPhoneDB_output$receptor_b
                                        ,function(i){
                                                if(i == "True") T else F }
                                        )
use_data(embryo_CellPhoneDB_output
         ,overwrite = TRUE)
