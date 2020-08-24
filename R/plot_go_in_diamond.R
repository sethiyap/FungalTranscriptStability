#' @title plot_go_in_diamond
#' @description allows to visualize GO term enrichment across different
#'   categories and samples
#' @param data_tbl tibble, a tibble of GO-terms and additional columns as 1.
#'   GO-term, 2. percent of genes enriched over background, 3. pvalue, 4.type
#'   ie. stable or unstable, 5. class ie. ros, infected etc.
#' @return GO term enrichment plot
#' @examples
#' \dontrun{
#' if(interactive()){
#'  data <- readr::read_delim(pipe("pbpaste"), delim="\t", col_names =T)
#'  }
#' }
#' @rdname plot_go_in_diamond
#' @export
#' @importFrom dplyr mutate
#' @importFrom forcats as_factor
#' @importFrom ggplot2 ggplot aes labs facet_grid vars geom_point scale_size
#'   scale_fill_manual theme_bw scale_x_discrete scale_y_discrete guides
#'   guide_legend


plot_go_in_diamond <- function(data_tbl){

          colnames(data_tbl) <- c("term", "percent_over_bkgd", "pvalue", "type", "class")

          print(head(data_tbl))

          data_fct <- data_tbl %>%
                      dplyr::mutate(term=forcats::as_factor(term),
                                    type=forcats::as_factor(type),
                                    class=forcats::as_factor(class))

  gg <-  data_fct %>%
          ggplot2::ggplot(ggplot2::aes(y=term, x=type,ordered=TRUE))+
          ggplot2::labs(x="", y="")+
          ggplot2::facet_grid( ggplot2::vars(class), space="free", scales = "free")+
          ggplot2::geom_point(ggplot2::aes(size=percent_over_bkgd),color="black",shape=23, stroke=0.8)+
          ggplot2::geom_point(ggplot2::aes(size=percent_over_bkgd,alpha=-log10(pvalue), fill=type),shape=23)+
          ggplot2::scale_size(range = c(1,15))+
          ggplot2::scale_fill_manual(values=c("#FC4E07", "#00AFBB" ))+
          ggplot2::theme_bw()+
          ggplot2::scale_x_discrete(position="top")+
          ggplot2::scale_y_discrete(position="right")+
          ggplot2::guides(fill = ggplot2::guide_legend(title="",override.aes = list(size=8)),
           alpha=ggplot2::guide_legend(override.aes=list(size=5)))

  return(gg)


}

