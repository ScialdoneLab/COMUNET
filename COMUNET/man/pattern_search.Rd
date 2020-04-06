\name{pattern_search}
\alias{pattern_search}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to perform pattern search
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pattern_search(pattern_adj_matrix, weight_array, ligand_receptor_pair_df, nodes, dissimilarity = d_normWeightDiff)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pattern_adj_matrix}{
numeric matrix: an adjacency matrix for the pattern of interest.
}
  \item{weight_array}{
numeric array: array of weighted adjacency matrices with dimensions [number of nodes, number of nodes, number of ligand-receptor pairs].
}
  \item{ligand_receptor_pair_df}{
character dataframe: data frame with columns "pair", "ligand", "ligand_complex_composition", "receptor", "receptor_complex_composition". Column "pair" contains values in a form "ligand:receptor", i.e. ligand being at the first place, receptor being at the second place, e.g. "TNFSF13:TNFRSF17". Column "ligand" contains ligand names, e.g. "TNFSF13".  Column ligand_complex_composition" if ligand is a complex (e.g. "aXb2_complex"), contains genes in the ligand complex separated with a comma, e.g. "ITGAX,ITGB2", else contains empty string "". Column "receptor" contains receptor names, e.g. "TNFRSF17". Column "receptor_complex_composition" if receptor is a complex (e.g. "NKG2D_II_receptor"), contains genes in the receptor complex separated with a comma, e.g. "KLRK1,HCST", else contains empty string "".
}
  \item{nodes}{
character string vector: a vector with all cell types in the data
}
  \item{dissimilarity}{
function: dissimilarity measure. Default value: d_normWeightDiff
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a dataframe with two columns:

 \item{pair}{character string vector: ligand-receptor pair names.}
\item{dissimilarity}{numeric vector: dissimilarity of ligand-receptor pairs to the pattern of interest. Note here that for the calculation of the dissimilarity of a weight matrix to the pattern of interest, we don't take into consideration the actual weight of the edges in this weight matrix, but just check wether the edge is present (i.e. has a non-zero value) or absent (i.e. has a zero value).}

The data frame is sorted by increasing dissimilarity (i.e. similar patterns a the top
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
