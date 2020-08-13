#' @title stabilityofActivelyTranscribingGenes
#' @description Stability of actively transcribed genes (i.e. belonging to top
#'   10-15 percent of total genes by RNA polymerase II) is determined by comparing RNAP
#'   II occupancy levels and steady state mRNA levels.
#' @param rna_pol2_tbl a tibble of normalized RNAP II and mRNA values. With
#'   first column (character) of gene-id, second column (numeric) RNAP II values
#'   and third column (numeric) mRNA values.
#' @param percent_cutoff numeric, cutoff to determine actively transcribing
#'   genes, i.e. to consider top 10 or top 15percent genes. Default: 15
#' @param write_output logical, to write output in text file or not. Default:
#'   FALSE
#' @param output_name character, if \code {write_output=TRUE} then determine the
#'   name of output file. Default: 'sample'
#' @return A tibble of stable and unstable genes determined from mRNA and RNAP
#'   II values.
#' @details mRNA levels should be normalized to gene lengths as mRNA values of
#'   genes within the given condition are considered to determine stability.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  rna_pol2_tbl <- readr::read_delim(pipe("pbpaste"), delim="\t", col_names = TRUE)
#'  }
#' }
#' @rdname stabilityofActivelyTranscribingGenes
#' @export
#' @importFrom dplyr rename all_of filter group_by mutate percent_rank ungroup
#'   if_else select arrange
#' @importFrom tidyr gather nest spread unnest
#' @importFrom ggplot2 ggplot aes geom_histogram theme_bw ggtitle labs theme
#'   element_text
#' @importFrom grid unit
#' @importFrom writexl write_xlsx
stabilityofActivelyTranscribingGenes <- function(rna_pol2_tbl, percent_cutoff= 15,write_output=FALSE,output_name="sample",binwidth=0.4){

  column1 <- rna_pol2_tbl %>% colnames() %>% .[1]

  rna_pol2_tbl <- rna_pol2_tbl %>%
                  dplyr::rename(gene_id= dplyr::all_of(column1))

  percentile_cutoff <- (100-percent_cutoff)/100

  dat_percentile <- rna_pol2_tbl %>%
                dplyr::filter((.[[2]] >= quantile(.[[2]], percentile_cutoff)) |
                              (.[[3]]>= quantile(.[[3]], percentile_cutoff)))

  dat_melt <- dat_percentile %>%
    tidyr::gather(data_type, value,-gene_id) %>%
    dplyr::group_by(data_type) %>%
    dplyr::mutate(log2value=log2(value+0.01), rank=dplyr::percent_rank(log2value)) %>%
    dplyr::ungroup() %>%
    tidyr::nest(value_col=c(value,log2value,rank)) %>%
    tidyr::spread(key = data_type, value = value_col) %>%
    tidyr::unnest(names_sep = '_') %>%
    dplyr::mutate(rank_ratio=log2(.[[7]]/.[[4]]),
                  stability=dplyr::if_else(rank_ratio>=1, "stable",
                                  dplyr::if_else(rank_ratio<=-1, "unstable", "no_change")))

  gg_1 <-  dat_melt %>%
            dplyr::select(c(rank_ratio, stability))%>%
            dplyr::arrange(desc(rank_ratio)) %>%
            dplyr::mutate(stability=factor(stability, levels=unique(stability))) %>%
            ggplot2::ggplot(ggplot2::aes(rank_ratio, fill=stability))+
            ggplot2::geom_histogram(alpha=0.8, color="black", binwidth = binwidth)+
            ggplot2::theme_bw()+
            ggplot2::scale_fill_manual(values=c("#FC4E07","bisque","#00AFBB"))+
            ggplot2::ggtitle(paste0("Percentile rank distribution","(top", percent_cutoff,"%)"))+
            ggplot2::labs(x="log2(percentile_rank_ratio)", y="Number of genes")+
            ggplot2::theme(
                    axis.text  = ggplot2::element_text( color="black",size=12),
                    axis.title   = ggplot2::element_text( color="black",size=12),
                    legend.title = ggplot2::element_text(face="bold", color="black",size=12),
                    legend.key.size = grid::unit(1.5,"line"),
                    legend.text  = ggplot2::element_text(face="bold", color="black",size=12))


  print(gg_1)
  print(table(dat_melt$stability))

        if(write_output==TRUE){
          writexl::write_xlsx(dat_melt,path = paste(output_name, "_stabilityOfActivelyTranscribingGenes.xlsx", sep=""),col_names = T)
        }
        else{
          return(dat_melt)
        }
}
