#' CellPhoneDB complex_input.csv file.
#'
#' @format A data frame with 112 rows and 19 variables.
#'
#'  @source \url{https://github.com/Teichlab/cellphonedb/blob/master/cellphonedb/src/core/data/complex_input.csv}
#'
"complex_input"

# define input.path
input.path <- getwd()

options(stringsAsFactors = F)

# read in complex_input.csv files
complex_input <- read.csv(paste0(input.path
                                 ,"/complex_input.csv"))
complex_input$complex_name <- gsub("_"
                                   ," "
                                   ,complex_input$complex_name
                                   )

use_data(complex_input
         ,overwrite = TRUE)
