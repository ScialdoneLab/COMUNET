#' @rdname plot_dissimilarity_heatmaps
#' @export
#'
#' @title
#' Plots heatmaps of dissimilarity between two conditions
#'
#' @description
#' Plots heatmaps of dissimilarity between two conditions.
#'
#' In each heatmap, interacting protein pairs from condition 1 are in the rows, and those from condition 2 are in the columns.
#' On the first heatmap, both rows and columns are clustered by the pairwise dissimilarity.
#' On the second heatmap, rows and columns are sorted by presence (shared or unshared) of an interacting protein pairs in both conditions
#' and by decreasing pairwise dissimilarity (in shared interacting pairs).
#'
#' Labels of interacting pairs are colored by their presence in both conditions (shared or unshared).
#' The heatmap is colored by the values of the pairwise dissimilarity.
#'
#' @author
#' Maria Solovey \email{maria.solovey@helmholtz-muenchen.de}
#'
#' @param dissim_cond1_cond2 Numeric matrix: pairwise dissimilarity between all ligand-receptor pairs in the two conditions (condition 1 in rows, condition 2 in columns).
#'
#' @param sorted_LRP_df Data frame with columns:
#'
#' \itemize{
#'
#' \item "pair" (character string): names of ligand-receptor pairs in the same form as they are in ligand_receptor_pair_df$pair.
#' \item "presence" (character string): name of the condition in which the ligand-receptor pair is present or "shared" if present in both conditions.
#' \item "dissimilarity" (numeric): dissimilarity value between the topology of ligand-receptor pair graph in two conditions.
#'
#' The smaller the dissimilarity value, the more similar is the graph topology between the two conditions.
#' If a ligand-receptor pair is present only in one of the conditions, the dissimilarity is equal to 1.
#' }
#'
#' @param cond1_name Character string: sample name for condition 1.
#'
#' @param cond2_name Character string: sample name for condition 2.
#'
#' @param colors_lrp Character string vector: colours for ligand-receptor labels. Default: green for shared, black for unshared.
#'
#' @param show_legend Logical: parameter of the \code{\link[HeatmapAnnotation:ComplexHeatmap]{HeatmapAnnotation}} function. Default is TRUE.
#'
#' @param row_names_fontsize Numeric: parameter of the \code{\link[Heatmap:ComplexHeatmap]{Heatmap}} function. Default is 5.
#'
#' @param colomn_names_fontsize Numeric: parameter of the \code{\link[Heatmap:ComplexHeatmap]{Heatmap}} function. Default is 5.
#'
#' @param show_column_names Logical: parameter of the \code{\link[Heatmap:ComplexHeatmap]{Heatmap}} function. Default is TRUE.
#'
#' @param show_row_names Logical: parameter of the \code{\link[Heatmap:ComplexHeatmap]{Heatmap}} function. Default  is TRUE.
#'
#' @param width Object of class "unit": parameter of the \code{\link[HeatmapAnnotation:ComplexHeatmap]{HeatmapAnnotation}} function. Default is unit(0.1, "cm").
#'
#' @param legend_direction Character string: parameter of the \code{\link[Heatmap:ComplexHeatmap]{Heatmap}} function. Default is "horizontal".
#'
#' @param legend_width Object of class "unit": parameter of the \code{\link[Heatmap:ComplexHeatmap]{Heatmap}} function. Default is unit(5, "cm").
#'
#' @param title_position Character string: parameter of the \code{\link[Heatmap:ComplexHeatmap]{Heatmap}} function. Default is "lefttop".
#'
#' @param row_dend_side Character string: parameter of the \code{\link[Heatmap:ComplexHeatmap]{Heatmap}} function. Default is "left"
#'
#' @param column_dend_side Character string: parameter of the \code{\link[Heatmap:ComplexHeatmap]{Heatmap}} function. Default is "top".
#'
#' @param column_title_side Character string: parameter of the \code{\link[draw:ComplexHeatmap]{draw}} function. Default is "top".
#'
#' @param heatmap_legend_side Character string: parameter of the \code{\link[draw:ComplexHeatmap]{draw}} function. Default is "bottom".
#'
#' @param ... Any other parameters of the \code{\link[Heatmap:ComplexHeatmap]{Heatmap}} function.
#'
#' @return
#' Two heatmaps:
#'
#' \itemize{
#'   \item{clustered heatmap}{
#'
#'   A heatmap with ligand-receptor pairs clustered by their pairwise dissimilarity.
#'   }
#'
#'   \item{sorted heatmap}{
#'
#'   A heatmap with ligand-receptor pairs sorted by decreasing dissimilarity of shared ligand-receptor pairs.
#'   }
#'   }
#'
#' @examples
#' # load comparative analysis
#' data("comparative_analysis")
#'
#' plot_dissimilarity_heatmaps(dissim_cond1_cond2 = comparative_analysis$dissim_cond1_cond2,
#'                             sorted_LRP_df = comparative_analysis$sorted_LRP_df,
#'                             cond1_name = "AML328-d0",
#'                             cond2_name = "AML328-d29")
#'
plot_dissimilarity_heatmaps <- function(dissim_cond1_cond2
                                        ,sorted_LRP_df
                                        ,cond1_name
                                        ,cond2_name
                                        ,colors_lrp = c("green"
                                                        ,"black")
                                        ,show_legend = TRUE
                                        ,row_names_fontsize = 5
                                        ,colomn_names_fontsize = 5
                                        ,show_column_names = TRUE
                                        ,show_row_names = TRUE
                                        ,width = unit(0.1, "cm")
                                        ,legend_direction = "horizontal"
                                        ,legend_width = unit(5, "cm")
                                        ,title_position = "lefttop"
                                        ,row_dend_side = "left"
                                        ,column_dend_side = "top"
                                        ,column_title_side = "top"
                                        ,heatmap_legend_side = "bottom"
                                        ,...
){
        # colour vector for ligand-receptor pairs in condition 1
        cond1_LRP_colors <- sapply(rownames(dissim_cond1_cond2)
                                   ,function(i){
                                           if(as.character(sorted_LRP_df[i,]$presence) == "shared") colors_lrp[1]
                                           else colors_lrp[2]
                                   })
        # colour vector for ligand-receptor pairs in condition 2
        cond2_LRP_colors <- sapply(colnames(dissim_cond1_cond2)
                                   ,function(i){
                                           if(as.character(sorted_LRP_df[i,]$presence) == "shared") colors_lrp[1]
                                           else colors_lrp[2]
                                   })
        # row annotations
        haRow = ComplexHeatmap::HeatmapAnnotation(
                df = data.frame(ligand_receptor_pair = c("shared"
                                                         ,"not shared")
                )
                ,col = list(ligand_receptor_pair = c("shared" = colors_lrp[1]
                                                     ,"not shared" = colors_lrp[2]))
                ,which = "row"
                ,width = width
                ,show_legend = show_legend
        )
        h_clustered <- ComplexHeatmap::Heatmap(
                dissim_cond1_cond2
                ,cluster_rows = T
                ,row_names_gp = gpar(fontsize = row_names_fontsize
                                     ,col = cond1_LRP_colors
                )
                ,cluster_columns = T
                ,column_names_gp = gpar(fontsize = colomn_names_fontsize
                                        ,col = cond2_LRP_colors
                )
                ,name = "dissimilarity coefficient"
                ,show_column_names = show_column_names
                ,show_row_names = show_row_names
                ,heatmap_legend_param = list(legend_direction = legend_direction
                                             ,legend_width = legend_width
                                             ,title_position = title_position
                )
                ,row_title = cond1_name
                ,row_dend_side = row_dend_side
                ,column_title = cond2_name
                ,column_dend_side = column_dend_side
                ,...
        )
        heatmap_comparative_clustered <- ComplexHeatmap::draw(
                h_clustered + haRow
                ,column_title = "Comparative analysis: clustered"
                ,column_title_side = column_title_side
                ,heatmap_legend_side = heatmap_legend_side
        )

        h_sorted <- ComplexHeatmap::Heatmap(
                dissim_cond1_cond2
                ,cluster_rows = F
                ,row_order = {
                        sub_sortedList <- sorted_LRP_df[rownames(sorted_LRP_df) %in% rownames(dissim_cond1_cond2),]
                        rownames(sub_sortedList[order(sub_sortedList$presence),])
                }
                ,row_names_gp = gpar(fontsize = 5
                                     ,col = cond1_LRP_colors
                )
                ,cluster_columns = F
                ,column_order = {
                        sub_sortedList <- sorted_LRP_df[rownames(sorted_LRP_df) %in% colnames(dissim_cond1_cond2),]
                        rownames(sub_sortedList[order(sub_sortedList$presence),])
                }
                ,column_names_gp = gpar(fontsize = 5
                                        ,col = cond2_LRP_colors
                )
                ,name = "dissimilarity coefficient"
                ,show_column_names = show_column_names
                ,show_row_names = show_row_names
                ,heatmap_legend_param = list(legend_direction = legend_direction
                                             ,legend_width = legend_width
                                             , title_position = title_position
                )
                ,row_title = cond1_name
                ,row_dend_side = row_dend_side
                ,column_title = cond2_name
                ,column_dend_side = column_dend_side
                ,...
        )
        heatmap_comparative_sorted <- ComplexHeatmap::draw(
                h_sorted + haRow
                ,column_title = "Comparative analysis: sorted"
                ,column_title_side = column_title_side
                ,heatmap_legend_side = heatmap_legend_side
        )
}
