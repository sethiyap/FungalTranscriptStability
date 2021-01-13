#' @title change in stability
#' @description Change in stability determined by applying linear regression on
#'   RNAP II and mRNA levels in two conditions.
#' @param rna_pol2_data tibble, with 5 columns, 1. gene_id, 2. pol2_control 3.
#'   pol2_test 4. rna_control 5. rna_test
#' @param percent_cutoff numeric, cutoff to determine actively transcribing
#'   genes, i.e. to consider top 10 or top 15 percent genes. Default: 15
#'  @param plot_label Title for plots. For eg. "test Vs Control" or "ros Vs unstressed"
#' @return a tibble of all genes denoted as being stabilized or destabilized in
#'   test condition
#' @details Linear regression and studentized residual method is used to classify the genes based on their
#'   fold-change in mRNA and RNAP II levels. Only actively transcribing genes
#'   (by RNAP II) in either test or control are considered.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  rna_pol2_data <- readr::read_delim(pipe("pbpaste"), delim="\t", col_names = TRUE)
#'  }
#' }
#' @rdname change_in_stability
#' @export
#' @importFrom rlang sym
#' @importFrom dplyr filter mutate percent_rank select if_else mutate_if
#'   summarise
#' @importFrom ggplot2 ggplot aes geom_point geom_hline theme_bw
#'   scale_color_manual ggtitle ylab xlab theme element_text xlim ylim geom_text

change_in_stability <- function(rna_pol2_data, percent_cutoff=15, plot_label="sample"){

  gene_id <- rna_pol2_data %>% colnames() %>% .[1] %>% rlang::sym()

  pol2_var1 <-  rna_pol2_data %>% colnames() %>% .[2] %>% rlang::sym()
  pol2_var2 <- rna_pol2_data %>% colnames() %>% .[3] %>% rlang::sym()

  rna_var1 <- rna_pol2_data %>% colnames() %>% .[4] %>% rlang::sym()
  rna_var2 <- rna_pol2_data %>% colnames() %>% .[5] %>% rlang::sym()

  percentile_cutoff <- (100-percent_cutoff)/100


  tbl_pol2_percent <- rna_pol2_data %>%
                          dplyr::filter(.[[2]] >=  quantile(.[[2]], percentile_cutoff) |
                                      .[[3]] >=  quantile(.[[3]], percentile_cutoff)) %>%
                          dplyr::mutate(LFC_pol2=log2(!!pol2_var2/!!pol2_var1),
                                        LFC_rna=log2(!!rna_var2/!!rna_var1)) %>%
                          dplyr::filter(!is.na(LFC_pol2) &
                                          !is.na(LFC_rna)) %>%
                          replace(is.na(.), 0) %>%
                          dplyr::mutate(percentile_rna = dplyr::percent_rank(LFC_rna),
                                        percentile_pol2=dplyr::percent_rank(LFC_pol2))

  mod <- lm(percentile_rna~percentile_pol2, data=tbl_pol2_percent)

  pol2_rna_studentised <- tbl_pol2_percent %>%
                            dplyr::mutate(studentised=rstudent(mod),
                                          pvalue=pnorm(abs(studentised),
                                          sd = sd(studentised),lower.tail = F),
                            category=dplyr::if_else(studentised>=1, "stabilized",
                                     dplyr::if_else(studentised<= -1, "destabilized", "no-change")))

  print(pol2_rna_studentised)

  gg_st <- ggplot2::ggplot(pol2_rna_studentised,ggplot2::aes(x=1:nrow(pol2_rna_studentised),studentised,colour = category))+
    ggplot2::geom_point(size=2.5,alpha=0.7)+ggplot2::geom_hline(yintercept = c(-1,1),lwd=1.1)+
    ggplot2::theme_bw()+
    ggplot2::scale_color_manual(values=c("#00AFBB", "bisque3", "#FC4E07"))+
    ggplot2::ggtitle(plot_label)+
    ggplot2::ylab("Studentised residual")+ggplot2::xlab("No. of Genes")+
    ggplot2::theme( axis.text = ggplot2::element_text(color="black",size=12),
                    axis.title.y=ggplot2::element_text(face="bold", color="black",size=14),
                    axis.title.x=ggplot2::element_text(face="bold", color="black",size=14),
                    legend.text=ggplot2::element_text(face="bold", color="black",size=12))
  print(gg_st)

  cors <- pol2_rna_studentised %>%
              dplyr::select(!!gene_id,LFC_pol2, LFC_rna) %>%
              dplyr::  mutate_if(is.numeric, ~ifelse(abs(.) == Inf,NA,.))%>%
              replace(is.na(.),0) %>%
              dplyr::summarise(cor=round(cor(LFC_pol2, LFC_rna), 2),
                                max=max(c(LFC_pol2, LFC_rna)),
                                min=min(c(LFC_pol2, LFC_rna)))

  max_lim <- max(cors$max)
  min_lim <- min(cors$min)
  x_cor <- max_lim - 1
  y_cor <- min_lim + 1

  gg_lfc <- ggplot2::ggplot(pol2_rna_studentised,ggplot2::aes(x=LFC_pol2, LFC_rna))+
    ggplot2::geom_point(size=1.5,alpha=0.7,ggplot2::aes(color=category))+
    ggplot2::ggtitle(plot_label)+
    ggplot2::scale_color_manual( values=c("#00AFBB", "bisque3", "#FC4E07"))+
    ggplot2::xlim(min_lim, max_lim)+
    ggplot2::ylim(min_lim, max_lim)+
    ggplot2::theme_bw()+
    ggplot2::theme(axis.text = ggplot2::element_text(color="black",size=12),
                   axis.title.y=ggplot2::element_text(face="bold", color="black",size=14),
                   axis.title.x=ggplot2::element_text(face="bold", color="black",size=14),
                   legend.text=ggplot2::element_text(face="bold", color="black",size=12))+
    ggplot2::geom_text(data=cors, ggplot2::aes(label=paste("r=", cor, sep="")), x=x_cor, y=y_cor)

  print(gg_lfc)

  print(table(pol2_rna_studentised$category))

  return(pol2_rna_studentised)

}
