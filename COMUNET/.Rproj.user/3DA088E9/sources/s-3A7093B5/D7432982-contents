\name{plot_cluster_UMAP}
\alias{plot_cluster_UMAP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to plot UMAP of ligand-receptor pairs
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot_cluster_UMAP(ligandReceptorPairDF, dissim_matrix, lrp_clusters, seed = 100, cluster_colors = "", unclustered_LRP_color = "black", cluster_shape = "", point_size = 2, title = "umap of ligand-receptor pairs", legend.position="bottom", legend_direction = "horizontal", legend_title.position = "top", legend_label.position= "bottom", legend_byrow=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ligandReceptorPairDF}{
character dataframe: data frame with columns "pair", "ligand", "ligand_complex_composition", "receptor", "receptor_complex_composition". Column "pair" contains values in a form "ligand:receptor", i.e. ligand being at the first place, receptor being at the second place, e.g. "TNFSF13:TNFRSF17". Column "ligand" contains ligand names, e.g. "TNFSF13".  Column ligand_complex_composition" if ligand is a complex (e.g. "aXb2_complex"), contains genes in the ligand complex separated with a comma, e.g. "ITGAX,ITGB2", else contains empty string "". Column "receptor" contains receptor names, e.g. "TNFRSF17". Column "receptor_complex_composition" if receptor is a complex (e.g. "NKG2D_II_receptor"), contains genes in the receptor complex separated with a comma, e.g. "KLRK1,HCST", else contains empty string "".
}
  \item{dissim_matrix}{
numeric matrix: pairwise dissimilarity matrix. Note that rownames and colnames should be defined.
}
  \item{lrp_clusters}{
numeric vector: cluster assignment for each ligand-recetor pair. Note that the the vector should be named and ordered in the same way as the rownames and colnames of the dissim_matrix.
}
  \item{seed}{
random seed for UMAP. Default value: 100
}
  \item{cluster_colors}{
character string vector: vector with colours for each cluster. Default value: "".
}
  \item{unclustered_LRP_color}{
character string: colour for ligand-receptor pairs that remained unclustered (cluster assignment 0). Default value: "black"
}
  \item{cluster_shape}{
numeric vector: vector with shapes for each cluster. Default value: "".
}
  \item{point_size}{
numeric: size of points. Default value: 2
}
  \item{title}{
character string: title of the heatmap plot. Argument of the draw function. Default value: "umap of ligand-receptor pairs".
}
  \item{legend.position}{
character string: argument of the theme function. Default value: "bottom" guide_legend.
}
  \item{legend_direction}{
character string: argument of the guide_legend function. Default value: "horizontal".
}
  \item{legend_title.position}{
character string: argument of the guide_legend function. Default value: "top".
}
  \item{legend_label.position}{
character string: argument of the guide_legend function. Default value: "bottom".
}
  \item{legend_byrow}{
logical: argument of the guide_legend function. Default value: TRUE.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
UMAP of ligand-receptor pairs coloured and shaped by cluster
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
