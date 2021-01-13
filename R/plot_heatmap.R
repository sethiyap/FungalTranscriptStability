#' @title plot_heatmap
#' @description plot heatmap
#' @param dat_mat matrix of input data to be plotted, all the columns should be numeric.
#' @param color_key numeric vector, determine range for color key. Continuous colors for non-negative value. And discriminate color for negative values.
#' Default: c(-1, 0, 1)
#' @param top_annotation logical, to plot top annotation i.e. boxplot or not. Default: TRUE
#' @param additional_column matrix, to plot additional heatmap. Default: NULL
#' @param additional_column_key numeric vector, determine range for color key for additional heatmap. Default: c(-8, 0, 8)
#' @return heatmap object

#' @examples
#' \dontrun{
#' if(interactive()){
#'  dat_mat1 <- readr::read_delim(pipe("pbpaste"), delim="\t", col_names = TRUE)
#'  dd <- dat_mat1 %>% dplyr::select(-c("rank_ratio"))  %>% tibble::column_to_rownames(var="gene_id")
#'  add_col <-  dat_mat1 %>% dplyr::select(c("gene_id","rank_ratio")) %>% tibble::column_to_rownames(var="gene_id")
#'  }
#' }
#' @rdname plot_heatmap
#' @importFrom ComplexHeatmap HeatmapAnnotation anno_boxplot Heatmap
#' @importFrom grid unit gpar
#' @importFrom dplyr summarise_if
#' @importFrom circlize colorRamp2
plot_heatmap <- function(dat_mat,color_key=c(-1,0,1), top_annotation=TRUE,key_name="expression_level",show_row_names=FALSE,additional_column=NULL, additional_column_key=c(-8,0,8)){

  if(top_annotation==TRUE){
  ha = ComplexHeatmap::HeatmapAnnotation(value =
                                           ComplexHeatmap::anno_boxplot(dat_mat, height = grid::unit(3, "cm"),
                                                               box_width = 0.5,
                                                               outline = TRUE,
                                                               gp = grid::gpar(fill = c("grey"),fonsize=12),
                                                               axis_param =list(facing="inside",side="left")))
  }
  else{
    ha=NULL
  }
  max_dat <- max(dat_mat)
  min_dat <- min(dat_mat)
  avg     <- dat_mat %>%
              dplyr::summarise_if(is.numeric, mean, na.rm = TRUE) %>%
              min()

  if(min_dat < 0){
              col <- circlize::colorRamp2(color_key,col=c("#008837","#f7f7f7","#d01c8b"))
  }else{
              col <- circlize::colorRamp2(color_key,col=c("#ffeda0","#feb24c","#bd0026"))
  }

  ht1 <- ComplexHeatmap::Heatmap(dat_mat,
                                 name=(key_name),
                                 border = "black",
                                 col=col,
                                 top_annotation = ha,
                                 show_row_dend = FALSE,
                                 show_row_names = show_row_names,
                                 cluster_columns = FALSE,
                                 cluster_rows = TRUE,
                                 column_names_rot=0,
                                 use_raster = TRUE,
                                 column_names_gp = grid::gpar(fontsize = 12),
  )

  if(is.null(additional_column)==TRUE){
              ht_list = ht1
  }else{

          min_add <-  min(additional_column)
          max_add <- max(additional_column)


          ht_list = ht1+ ComplexHeatmap::Heatmap(additional_column, border = "black",
                                                      row_order = order(additional_column, decreasing = FALSE),
                                                      col = circlize::colorRamp2(additional_column_key,col=c("red","#f7f7f7","blue")),
                                                      name=paste(colnames(additional_column)),show_row_names = FALSE,
                                                      width = grid::unit(4, "mm"))
  }

  return(ht_list)

}
