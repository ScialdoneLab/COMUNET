\name{convert_CellPhoneDB_output}
\alias{convert_CellPhoneDB_output}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to convert CellPhoneDB output matrix into an array of weighted adjacency matrices
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
convert_CellPhoneDB_output(CellPhoneDB_output, complex_input, gene_input)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CellPhoneDB_output}{
significant_means.txt file in the output of CellPhoneDB
}
  \item{complex_input}{
complex_input.csv file. Download from CellPhoneDB
}
  \item{gene_input}{
gene_input.csv file. Download from CellPhoneDB
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a list of:

\item{weight_array}{
3D numeric array: array of weighted adjacency matrices with dimensions [number of nodes, number of nodes, number of ligand-receptor pairs]. Note that the function filters out empty weight arrays and the corresponding ligand-receptor pairs. First dimension: sending nodes (note: this dimension has all possible nodes, even if some of them are silent for a particular ligand-receptor pair). Second dimension: receiving nodes (note: this dimension has all possible nodes, even if some of them are silent for a particular ligand-receptor pair). Third dimension: ligand-receptor pairs. Note that the weight_array should contain dimnames: dimnames = list(nodes, nodes, ligand-receptor pairs).
}
\item{ligandReceptorPairDF}{
character dataframe: data frame with columns "pair", "ligand", "ligand_complex_composition", "receptor", "receptor_complex_composition". Column "pair" contains values in a form "ligand:receptor", i.e. ligand being at the first place, receptor being at the second place, e.g. "TNFSF13:TNFRSF17". Column "ligand" contains ligand names, e.g. "TNFSF13".  Column ligand_complex_composition" if ligand is a complex (e.g. "aXb2_complex"), contains genes in the ligand complex separated with a comma, e.g. "ITGAX,ITGB2", else contains empty string "". Column "receptor" contains receptor names, e.g. "TNFRSF17". Column "receptor_complex_composition" if receptor is a complex (e.g. "NKG2D_II_receptor"), contains genes in the receptor complex separated with a comma, e.g. "KLRK1,HCST", else contains empty string "".
}
\item{nodes}{
character string vector: a vector with all cell types in the data.
}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Maria Solovey
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
