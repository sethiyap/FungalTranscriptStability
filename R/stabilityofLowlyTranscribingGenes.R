#' @title stabilityofLowlyTranscribingGenes
#' @description Stability of lowly or non-transcribing genes is determined by
#'   comparing RNAP II occupancy levels and steady state mRNA levels
#' @param rna_pol2_tbl a tibble of normalized RNAP II and mRNA values. With
#'   first column (character) of gene-id, second column (numeric) RNAP II values
#'   and third column (numeric) mRNA values.
#' @param percent_cutoff numeric, cutoff to determine actively transcribing
#'   genes, i.e. to consider bottom 85 percent genes. Default: 85
#' @param exp_cutoff  numeric, minimum expression value in either RNAP II or mRNA
#'   to filter. Default: 2
#' @param plot_heatmap  logical, to plot heatmap of the z-score of RNAP II and
#'   mRNA values.If FALSE tibble will be returned. Default: FALSE
#' @param color_key  numeric vector, range for heatmap color key.  Default:
#'   c(-0.75, 0, 0.75)
#' @return quantile analysis plot, heatmap and tibble containing stable and
#'   unstable genes.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'  rna_pol2_tbl <- readr::read_delim(pipe("pbpaste"), delim="\t", col_names = TRUE)
#'  }
#' }
#'
#' @rdname stabilityofLowlyTranscribingGenes
#' @export
#' @importFrom dplyr rename all_of filter mutate ntile if_else arrange select
#' @importFrom forcats as_factor
#' @importFrom ggplot2 ggplot aes geom_point geom_jitter theme_bw
#'   scale_color_manual
#' @importFrom tibble column_to_rownames
#' @importFrom ComplexHeatmap Heatmap draw
#' @importFrom grid unit
#'
stabilityofLowlyTranscribingGenes <- function(rna_pol2_tbl,percent_cutoff=85, exp_cutoff=2, plot_heatmap=FALSE,color_key= c(-0.75, 0, 0.75)){

                column1 <- rna_pol2_tbl %>% colnames() %>% .[1]

                rna_pol2_tbl <- rna_pol2_tbl %>%
                                dplyr::rename(gene_id= dplyr::all_of(column1))

  percentile_cutoff <- (percent_cutoff)/100

  dat_percentile <- rna_pol2_tbl %>%
                      dplyr::filter((.[[2]] < quantile(.[[2]], percentile_cutoff)) &
                                    (.[[3]] < quantile(.[[3]], percentile_cutoff))) %>%
                      dplyr::filter(.[[2]] > exp_cutoff | .[[3]] > exp_cutoff ) %>%
                      dplyr::mutate(log2pol2= log2(.[[2]]+0.01), log2rna= log2(.[[3]]+0.01),
                                    pol2_quantile=dplyr::ntile(log2pol2, 4),
                                    rna_quantile=dplyr::ntile(log2rna, 4))


  dat_stability <- dat_percentile %>%
                      dplyr::mutate(stability=dplyr::if_else(pol2_quantile == 1 & rna_quantile == 4 |
                                                         pol2_quantile == 2 & rna_quantile == 4 |
                                                         pol2_quantile == 1 & rna_quantile == 3, "stable",
                                                       dplyr::if_else(
                                                              pol2_quantile == 4 & rna_quantile == 1 |
                                                              pol2_quantile == 4 & rna_quantile == 2 |
                                                              pol2_quantile == 3 & rna_quantile == 1, "unstable",
                                                              "no-change")))

  gp_plot <- dat_stability %>%
                dplyr::arrange(desc(stability)) %>%
                dplyr::mutate(stability=forcats::as_factor(stability)) %>%
                ggplot2::ggplot(ggplot2::aes(pol2_quantile, rna_quantile, color=stability))+
                ggplot2::geom_point()+
                ggplot2::geom_jitter(size=1.5) +
                ggplot2::theme_bw()+
                ggplot2::scale_color_manual(values=c("#00AFBB","#FC4E07","bisque2"))

  print(gp_plot)

  if(plot_heatmap==TRUE){

   ht_dat <-  dat_stability %>%
                dplyr::select(c(gene_id, log2pol2, log2rna, pol2_quantile)) %>%
                dplyr::mutate(pol2=scale(log2pol2), rna=scale(log2rna)) %>%
                dplyr::select(-c(log2pol2, log2rna)) %>%
                tibble::column_to_rownames(var="gene_id")

   ht1 <- ht_dat %>% dplyr::select(c(pol2,rna))
   ht2 <- ht_dat %>% dplyr::select(pol2_quantile)

   ht_plot <- plot_heatmap(ht1, top_annotation = FALSE, color_key = color_key)

   ht_list <- ComplexHeatmap::Heatmap(ht2, border = "black",
                                      col = structure(1:(length(table(ht2)))),
                                      show_row_names = TRUE,
                                      cluster_rows = FALSE,
                                      name="quantile", width = grid::unit(4, "mm"))
   ht <- ComplexHeatmap::draw(ht_list+ht_plot, split=ht2)
   return(ht)

  }else{
    return(dat_stability)
  }

}
