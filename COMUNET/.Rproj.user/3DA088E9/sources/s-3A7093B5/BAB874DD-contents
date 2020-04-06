\name{plotCommunicationGraph}
\alias{plotCommunicationGraph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Function to make a plot of communication graph for a given ligand-receptor pair.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotCommunicationGraph(LRP, weight_array, ligandReceptorPairDF, nodes, isPattern = F, isCluster = F, title = LRP, subtitle = "", node_color_palette = c("blue", "gray40", "red"), node.label.cex = 1, vertex.shape = "none", vertex.size = 15, edge.arrow.size = 0.5, edge.width.scale = 2.5, legend.size = 0.75, legend.gradient.x.coord = c(1.15,1.25,1.25,1.15), legend.gradient.y.coord = c(-0.5,-0.5,-1,-1), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{LRP}{
     character string: name of ligand-receptor pair for which the plot should be constructed. Note the the LRP name should be written in the same from as in ligandReceptorPairDF$pair. Alternatively, it can be a name of a cluster or a pattern.

}
  \item{weight_array}{
     3D numeric array or 2D numeric matrix: array of weighted adjacency matrices with dimensions [number of nodes, number of nodes, number of ligand-receptor pairs] or one adjacency matrix with dimensions [number of nodes, number of nodes]. First dimension: sending nodes (note: this dimension has all possible nodes, even if some of them are silent for a particular ligand-receptor pair). Second dimension: receiving nodes (note: this dimension has all possible nodes, even if some of them are silent for a particular ligand-receptor pair). Third dimension (for arrays): ligand-receptor pairs. Note that the weight_array should contain dimnames: dimnames = list(nodes, nodes, ligand-receptor pairs).

}
  \item{ligandReceptorPairDF}{
     character dataframe: data frame with columns "pair", "ligand", "ligand_complex_composition", "receptor", "receptor_complex_composition". Column "pair" contains values in a form "ligand:receptor", i.e. ligand being at the first place, receptor being at the second place, e.g. "TNFSF13:TNFRSF17". Column "ligand" contains ligand names, e.g. "TNFSF13".  Column ligand_complex_composition" if ligand is a complex (e.g. "aXb2_complex"), contains genes in the ligand complex separated with a comma, e.g. "ITGAX,ITGB2", else contains empty string "". Column "receptor" contains receptor names, e.g. "TNFRSF17". Column "receptor_complex_composition" if receptor is a complex (e.g. "NKG2D_II_receptor"), contains genes in the receptor complex separated with a comma, e.g. "KLRK1,HCST", else contains empty string "".
  }
  \item{nodes}{
     character vector: a vector with all cell types in the data.
  }
  \item{isPattern}{
     logical: should a graph of a pattern (and not a particular ligand-receptor pair) be plotted. Default value = F.
  }
  \item{isCluster}{
     logical: should a graph of a cluster (and not a particular ligand-receptor pair) be plotted. Default value = F.
  }
  \item{title}{
     character string: title of the plot. Default value: LRP.
  }
  \item{subtitle}{
     character string:  subtitle of the plot. Default value: ""
  }
  \item{node_color_palette}{
     character string vector: vector of colours for nodes. Default values: c("blue", "gray40", "red"). The colours will be used to construct a colour gradient for the node delta degree: blue representing a receiving node (negative delta degree), grey representing a neutral node (zero delta degree), red representing a sending node (positive delta degree).
  }
  \item{node.label.cex}{
     numeric: size of node labels. Default value: 1
  }
  \item{vertex.shape}{
     charater string: shape of node. Default value: "none". See same argument of plot.igraph().
  }
  \item{vertex.size}{
     integer: size of node. Default value: 15. See same argument of plot.igraph().
  }
  \item{edge.arrow.size}{
     numeric: size of arrow edge. Default value: 0.5. See same argument of plot.igraph().
  }
  \item{edge.width.scale}{
     numeric: scaling factor for edge width. Default value: 2.5. See same argument of plot.igraph().
  }
  \item{legend.size}{
     numeric: size of the legend text. Default value: 0.75
  }
  \item{legend.gradient.x.coord}{
     numeric vector: x coordinates of the label position. Default value: c(1.15,1.25,1.25,1.15).
  }
  \item{legend.gradient.y.coord}{
     numeric vector: y coordinates of the label position. Default value: c(-0.5,-0.5,-1,-1).
  }
  \item{...}{
     other plot.igraph parameters
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
graph plot
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
