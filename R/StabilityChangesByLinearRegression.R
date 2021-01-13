
tbl_rna <- readr::read_delim(pipe("pbpaste"), delim="\t", col_names = TRUE)

tbl_pol2 <- readr::read_delim(pipe("pbpaste"), delim="\t", col_names = TRUE)

StabilityChangesByLinearRegression <- function(tbl_rna, tbl_pol2, plot_label="sample"){

  column1 <- tbl_rna %>% colnames() %>% .[1]
  tbl_rna <- tbl_rna %>%
              dplyr::rename(Gene= dplyr::all_of(column1))


    rna_LFC <- tbl_rna %>%
                dplyr::mutate(LFC_RNA=log2(.[[3]]/.[[2]]),
                              LFC_RNA=ifelse(is.infinite(LFC_RNA), 0, LFC_RNA)) %>%
                replace(is.na(.), 0)

    column2 <- tbl_pol2 %>% colnames() %>% .[1]
    tbl_pol2 <- tbl_pol2 %>%
                dplyr::rename(Gene= dplyr::all_of(column2))

    pol2_transcribing <- tbl_pol2 %>%
                            dplyr::filter(.[[2]] >=  quantile(.[[2]], 0.85) |
                                          .[[3]] >=  quantile(.[[3]], 0.85))  %>%
                            dplyr::mutate(LFC_Pol2=log2(.[[3]]/.[[2]]))

      pol2_rna <- dplyr::left_join(pol2_transcribing,rna_LFC) %>%
                  replace(is.na(.), 0) %>%
                  dplyr::mutate(percentile_rna = dplyr::percent_rank(LFC_RNA),
                                percentile_pol2=dplyr::percent_rank(LFC_Pol2))

      mod <- lm(percentile_rna~percentile_pol2, data=pol2_rna)

      pol2_rna_studentised <- pol2_rna %>%
                    dplyr::select(c(Gene, LFC_RNA, LFC_Pol2, percentile_rna, percentile_pol2)) %>%
                    dplyr::mutate(studentised=rstudent(mod),
                                  pvalue=pnorm(abs(studentised),sd = sd(studentised),lower.tail = F),
                                  category=dplyr::if_else(studentised>=1, "stabilized",
                                                          dplyr::if_else(studentised<= -1, "destabilized", "no-change")))
      print(pol2_rna)
      #-- Plot Studentised Residual
      gg_st <- ggplot2::ggplot(pol2_rna_studentised,ggplot2::aes(x=1:nrow(pol2_rna_studentised),studentised,colour = category))+
        ggplot2::geom_point(size=2.5,alpha=0.7)+ggplot2::geom_hline(yintercept = c(-1,1),lwd=1.1)+
        ggplot2::theme_bw()+
        ggplot2::scale_color_manual(values=c("cyan3", "blue", "salmon"))+
        ggplot2::ggtitle(plot_label)+
        ggplot2::ylab("Studentised residual")+ggplot2::xlab("No. of Genes")+
        ggplot2::theme( axis.text = ggplot2::element_text(color="black",size=12),
                        axis.title.y=ggplot2::element_text(face="bold", color="black",size=14),
                        axis.title.x=ggplot2::element_text(face="bold", color="black",size=14),
                        legend.text=ggplot2::element_text(face="bold", color="black",size=12))
      print(gg_st)

      cors <- pol2_rna_studentised %>%
              dplyr::summarise(cor=round(cor(LFC_Pol2, LFC_RNA), 2), max=max(c(LFC_Pol2, LFC_RNA)), min=min(c(LFC_Pol2, LFC_RNA)))

      max_lim <- max(cors$max)
      min_lim <- min(cors$min)
      x_cor <- max_lim - 1
      y_cor <- min_lim + 1

      gg_lfc <- ggplot2::ggplot(pol2_rna_studentised,ggplot2::aes(x=LFC_Pol2,LFC_RNA))+
        ggplot2::geom_point(size=2.5,alpha=0.7,ggplot2::aes(color=category))+
        ggplot2::ggtitle(plot_label)+
        ggplot2::scale_color_manual( values=c("cyan3", "blue", "salmon"))+
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
