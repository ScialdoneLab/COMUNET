\name{comparativeAnalysis}
\alias{comparativeAnalysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to perform comparative analysis
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
comparativeAnalysis(cond1_weight_array, cond2_weight_array, cond1_ligandReceptorPairDF, cond2_ligandReceptorPairDF, cond1_nodes, cond2_nodes = cond1_nodes, cond1_name, cond2_name, dissimilarity = d_normWeightDiff)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cond1_weight_array}{
numeric array: array of weighted adjacency matrices for condition 1. Note that the function filters out empty weight arrays and the corresponding ligand-receptor pairs.
}
  \item{cond2_weight_array}{
numeric array: array of weighted adjacency matrices for condition 2. Note that the function filters out empty weight arrays and the corresponding ligand-receptor pairs
}
  \item{cond1_ligandReceptorPairDF}{
character dataframe          data frame with columns "pair", "ligand", "ligand_complex_composition", "receptor", "receptor_complex_composition" for condition 1.
}
  \item{cond2_ligandReceptorPairDF}{
character dataframe: data frame with columns "pair", "ligand", "ligand_complex_composition", "receptor", "receptor_complex_composition" for condition 2
}
  \item{cond1_nodes}{
character string vector: a vector with all cell types in the data. Please note that comparativeAnalysis expects nodes (i.e. cell types) to be the same in both conditions.
}
  \item{cond2_nodes}{
character string vector: a vector with all cell types in the data. Default value: cond1_nodes. Please note that comparativeAnalysis expects nodes (i.e. cell types) to be the same in both conditions.
}
  \item{cond1_name}{
character string: sample name for condition 1
}
  \item{cond2_name}{
character string: sample name for condition 2
}
  \item{dissimilarity}{
function: dissimilarity function. Default value is d_normWeightDiff
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a list of:

\item{sortedLRP_DF}{
dataframe with columns: "pair" (character string): names of ligand-receptor pairs in the same form as they are in ligandReceptorPairDF$pair; "presence" (character string): name of the condition in which the ligand-receptor pair is present or "shared" if present in both conditions; dissimilarity (numeric): dissimilarity value between the topology of ligand-receptor pair graph in two conditions. The smaller the dissimilarity value, the more similar is the graph topology between the two conditions. If a ligand-receptor pair is present only in one of the conditions, the dissimilarity is equal to 1.
}
\item{dissim_cond1_cond2}{
numeric matrix: pairwise dissimilarity between all ligand-receptor pairs in the two conditions (condition 1 in rows, condition 2 in columns).
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
