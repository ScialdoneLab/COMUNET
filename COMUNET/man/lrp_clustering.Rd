\name{lrp_clustering}
\alias{lrp_clustering}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to cluster ligand-receptor pairs
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
lrp_clustering(weight_array, ligand_receptor_pair_df, nodes, dissimilarity = d_normWeightDiff, min_cluster_size = 6)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{weight_array}{
numeric array: array of weighted adjacency matrices with dimensions [number of nodes, number of nodes, number of ligand-receptor pairs].
}
  \item{ligand_receptor_pair_df}{
character dataframe: data frame with columns "pair", "ligand", "ligand_complex_composition", "receptor", "receptor_complex_composition". Column "pair" contains values in a form "ligand:receptor", i.e. ligand being at the first place, receptor being at the second place, e.g. "TNFSF13:TNFRSF17". Column "ligand" contains ligand names, e.g. "TNFSF13".  Column ligand_complex_composition" if ligand is a complex (e.g. "aXb2_complex"), contains genes in the ligand complex separated with a comma, e.g. "ITGAX,ITGB2", else contains empty string "". Column "receptor" contains receptor names, e.g. "TNFRSF17". Column "receptor_complex_composition" if receptor is a complex (e.g. "NKG2D_II_receptor"), contains genes in the receptor complex separated with a comma, e.g. "KLRK1,HCST", else contains empty string "".
}
  \item{nodes}{
character string vector: a vector with all cell types in the data.
}
  \item{dissimilarity}{
function: dissimilarity function. Default value is d_normWeightDiff
}
  \item{min_cluster_size}{
numeric: minimum number of ligand-receptor pairs per cluster. Default value: 6.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a list of:

\item{dissim_matrix}{
numeric matrix: pairwise dissimilarity between all ligand-receptor pairs.
}
\item{clusters}{
numeric vector: cluster assignment for each ligand-receptor pair.
}
\item{weight_array_by_cluster}{
numeric array: array of weighted adjacency matrices with dimensions [number of nodes, number of nodes, number of ligand-receptor clusters]. Each edge weight in a cluster is calculated as the mean of the weights of this edge in all ligand-receptor pairs that belong to this cluster.
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
