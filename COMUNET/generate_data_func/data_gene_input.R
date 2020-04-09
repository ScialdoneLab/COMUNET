#' CellPhoneDB gene_input.csv file.
#'
#' @format A data frame with 1252 rows and 4 variables.
#'
#'  @source \url{https://github.com/Teichlab/cellphonedb/blob/master/cellphonedb/src/core/data/gene_input.csv}
#'
"gene_input"

# define input.path
input.path <- getwd()

options(stringsAsFactors = F)

# read in cgene_input.csv files
gene_input <- read.csv(paste0(input.path
                                 ,"/gene_input.csv"))
gene_input$complex_name <- gsub("_"
                                   ," "
                                   ,gene_input$complex_name
)

use_data(gene_input
         ,overwrite = TRUE)
