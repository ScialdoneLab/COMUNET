#' @rdname plot_cluster_heatmap
#' @export
#'
#' @title
#' Plots heatmap for clustered ligand-receptor pairs
#'
#' @description
#' Plots a heatmap for clustered ligand-receptor pairs.
#'
#' Rows and columns are sorted by cluster number, the color of the heatmap represents dissimilarity between two ligand-receptor pairs.
#'
#' @author
#' Maria Solovey \email{maria.solovey@helmholtz-muenchen.de}
#'
#' @param       dissim_matrix Numeric matrix: pairwise dissimilarity matrix.
#'
#' Note that row names and column names should be defined.
#'
#' @param       lrp_clusters Numeric vector: cluster assignment for each ligand-receptor pair.
#'
#' Note that the the vector should be named and ordered in the same way as the row names and column names of the dissim_matrix.
#'
#' @param       cluster_colors Character string vector: vector with colours for each cluster. Default value: "".
#'
#' @param       unclustered_LRP_color Character string: colour for ligand-receptor pairs that remained unclustered (cluster assignment 0). Default value: "black".
#'
#' @param       show_colomn_names Logical: argument of the \code{\link[Heatmap:ComplexHeatmap]{Heatmap}} function. Default value: FALSE.
#'
#' @param       show_row_names Logical: argument of the \code{\link[Heatmap:ComplexHeatmap]{Heatmap}} function. Default value: FALSE.
#'
#' @param       legend_direction Character string: argument of the \code{\link[Heatmap:ComplexHeatmap]{Heatmap}} function. Default value: "horizontal".
#'
#' @param       legend_width Object of class "unit": argument of the \code{\link[Heatmap:ComplexHeatmap]{Heatmap}} function. Default is: unit(5, "cm").
#'
#' @param       title_position Character string: argument of the \code{\link[Heatmap:ComplexHeatmap]{Heatmap}} function. Default value: "lefttop".
#'
#' @param       column_title Character string: title of the heatmap plot. Argument of the \code{\link[draw:ComplexHeatmap]{draw}} function. Default value: "Clustering of ligand-receptor pairs".
#'
#' @param       column_title_side Character string:  position of the title of the heatmap plot. Argument of the \code{\link[draw:ComplexHeatmap]{draw}} function. Default value: "top".
#'
#' @param       heatmap_legend_side Character string: position of the legend of the heatmap plot. Argument of the \code{\link[draw:ComplexHeatmap]{draw}} function. Default value: "bottom".
#'
#' @param       ... Any arguments \code{\link[Heatmap:ComplexHeatmap]{Heatmap}} function takes.
#'
#' @return      Heatmap with ligand-receptor pairs ordered by their cluster assignment.
#'
#' @examples
#'
#' # load lrp_clusters
#' data("lrp_clusters")
#'
#' # plot heatmap
#' plot_cluster_heatmap(dissim_matrix = lrp_clusters$dissim_matrix, lrp_clusters = lrp_clusters$clusters)
#'
plot_cluster_heatmap <- function(dissim_matrix
                                 ,lrp_clusters
                                 ,cluster_colors = ""
                                 ,unclustered_LRP_color = "black"
                                 ,show_column_names = FALSE
                                 ,show_row_names = FALSE
                                 ,legend_direction = "horizontal"
                                 ,legend_width = unit(5, "cm")
                                 ,title_position = "lefttop"
                                 ,column_title = "Clustering of ligand-receptor pairs"
                                 ,column_title_side = "top"
                                 ,heatmap_legend_side = "bottom"
                                 ,...
){
        # check that dissim_matrix  has rownames and colnames
        if(length(rownames(dissim_matrix )) == 0) stop("ERROR: dissim_matrix  has no rownames. Please define rownames as ligand-receptor pairs")
        if(length(colnames(dissim_matrix )) == 0) stop("ERROR: dissim_matrix  has no rownames. Please define colnames as ligand-receptor pairs")

        # check that rownames and colnames of the dissim_matrix  are identical
        if(!identical(rownames(dissim_matrix )
                      ,colnames(dissim_matrix ))) stop("ERROR: rownames and colnames of the dissim_matrix  are not identical. Please make sure they are identical")

        # check that the lrp_clusters is named
        if(length(names(lrp_clusters)) == 0) stop("ERROR: lrp_clusters is not a named vector. Please nake sure the names are ligand-receptor pairs")

        # check that rownames and colnames of the dissim_matrix and names of lrp_clusters are identical
        if(!identical(rownames(dissim_matrix)
                      ,names(lrp_clusters))) stop("ERROR: rownamesod dissim_matrix and names of lrp_clusters are not identical. Please make identical")
        if(!identical(colnames(dissim_matrix)
                      ,names(lrp_clusters))) stop("ERROR: colnamesod dissim_matrix and names of lrp_clusters are not identical. Please make identical")

        # assign colours to clusters
        if(cluster_colors == "") cluster_colors <- structure(rainbow(length(levels(as.factor(lrp_clusters)))))
        names(cluster_colors) <- levels(as.factor(lrp_clusters))

        # if there are unclustered ligand-receptor pairs, colour them with a defined colour
        if(as.factor(0) %in% names(cluster_colors)) cluster_colors[names(cluster_colors) == as.factor(0)] <- unclustered_LRP_color

        # mark the as "unclustered"
        if(max(lrp_clusters) != 0) {
                annotation_legend_param <-  list(labels = c("unclustered"
                                                            ,1:max(lrp_clusters))
                )} else { # all ligand-receptor pairs remained unclustered
                        annotation_legend_param <-  list(labels = "unclustered"
                        )
                }

        # define column annotation
        haCol <- ComplexHeatmap::HeatmapAnnotation(
                df = data.frame(cluster = as.vector(lrp_clusters))
                ,col = list(cluster = cluster_colors)
                ,annotation_legend_param = annotation_legend_param
        )

        # heatmap ordered by clusters
        h1 <- ComplexHeatmap::Heatmap(
                dissim_matrix
                ,cluster_rows = F
                ,row_order = rownames(dissim_matrix)[order(as.numeric(lrp_clusters))]
                ,cluster_columns = F
                ,name = "dissimilarity coefficient"
                ,column_order = colnames(dissim_matrix)[order(as.numeric(lrp_clusters))]
                ,top_annotation = haCol
                ,show_column_names = show_column_names
                ,show_row_names = show_row_names
                ,heatmap_legend_param = list(legend_direction = legend_direction
                                             ,legend_width = legend_width
                                             ,title_position = title_position
                )
                ,...
        )

        heatmap_clusters <- ComplexHeatmap::draw(
                h1
                ,column_title = column_title
                ,column_title_side = column_title_side
                ,heatmap_legend_side = heatmap_legend_side
        )
}
