#' @title orf_len_distribution
#' @description displays distribution of orf lengths for give stable and
#'   unstable genes.
#' @param gff path or connection of gff file
#' @param genelist_stabilityInfo tibble, of genes with stability information.
#'   First_column: gene_id, second_column: stability ie. whether the gene is
#'   stable or unstable, third_column: condition ie. stressed or unstressed
#' @return pvalue of the orf lengths, boxplot of orf length distribution
#' @examples
#' \dontrun{
#' if(interactive()){
#'   gff_file <- "/Users/Pooja/Documents/PhD/thesis_figures/Chapter3/C_glabrata_CBS138_version_s02-m07-r38_features.gff"
#'  }
#' }
#' @rdname orf_len_distribution
#' @export
#' @importFrom GenomicFeatures makeTxDbFromGFF transcriptLengths
#' @importFrom tibble as_tibble
#' @importFrom dplyr select left_join mutate group_by summarize
#' @importFrom forcats as_factor
#' @importFrom ggplot2 ggplot aes facet_wrap geom_boxplot scale_fill_manual
#'   scale_y_log10 theme_bw labs scale_x_discrete theme element_blank
#'   element_text element_rect
#' @importFrom tidyr nest spread
#' @importFrom purrr map2 map_dbl
#' @importFrom broom tidy
orf_len_distribution <- function(gff, genelist_stabilityInfo){

  colnames(genelist_stabilityInfo) <- c("gene_id", "stability", "condition")

  # load gene features
  gff <- GenomicFeatures::makeTxDbFromGFF(gff, metadata = T, format="gff")
  cds_len <- GenomicFeatures::transcriptLengths(gff, with.cds_len = T) %>%
                tibble::as_tibble() %>% dplyr::select(c(gene_id,cds_len))

  # match gene length to the input genes
    gene_width <- dplyr::left_join(genelist_stabilityInfo, cds_len) %>%
                  dplyr::mutate(stability=forcats::as_factor(stability))


    print(gene_width)
  # calculate median length
  median = gene_width %>%
            dplyr::group_by(stability, condition) %>%
            dplyr::summarize(med = median(cds_len)) %>%
            tidyr::spread(stability,med)

  print(median)

  # plot boxplot
  gg <- gene_width %>%
    ggplot2::ggplot(ggplot2::aes(stability, cds_len, fill=stability))+
    ggplot2::facet_wrap(~condition)+
    ggplot2::geom_boxplot(notch = TRUE, notchwidth = 0.4, alpha=0.8, width=0.5) +
    ggplot2::scale_fill_manual(values=c("#FC4E07", "#00AFBB" ))+
    ggplot2::scale_y_log10()+
    ggplot2::theme_bw()+
    ggplot2::labs(x="", y="ORF lenght (in bp)")+
    ggplot2::scale_x_discrete(position="top")+
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_text(color="black",size=12),
          panel.border = ggplot2::element_rect(color="black", size=1),
          axis.title.y = ggplot2::element_text(color="black",size=12, face="bold"),
          strip.text   = ggplot2::element_text(color="black",size=12,face="bold"),
          panel.grid   = ggplot2::element_blank(),
          strip.background  = ggplot2::element_rect(fill = "white", color="black", size=1),
          legend.position  = "bottom"
    )
  print(gg)

  var1 <- unique(gene_width$stability %>% levels() %>% .[1])
  var2 <- unique(gene_width$stability %>% levels() %>% .[2])

  var1 <- rlang::sym(var1)
  var2 <- rlang::sym(var2)

  t_test <- gene_width %>%
    dplyr::group_by(stability, condition) %>%
    tidyr::nest() %>%
    tidyr::spread(key = stability, value = data) %>%
    dplyr::mutate(
      t_test = purrr::map2(!!var1,!!var2, ~ t.test(.x$cds_len, .y$cds_len) %>% broom::tidy())) %>%
    dplyr::mutate(pval = purrr::map_dbl(t_test , ~ .$`p.value`)) %>%
    dplyr::select(condition, pval)

  return(t_test)
}

