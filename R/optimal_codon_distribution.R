#' @title optimal_codon_distribution
#' @description Codon distribution based on pre-defined optimal codons for yeast
#'   species (Pechmann, 2013)
#' @param optimal_codon_table Codon table, containing pre-defined optimal and
#'   Non-optimal codons
#' @param genelist_stabilityInfo tibble, of genes with stability information.
#'   First_column: gene_id, second_column: stability ie. whether the gene is
#'   stable or unstable, third_column: condition ie. stressed or unstressed
#' @param ref_sequence_file path or connection of coding nucleotides sequences
#'   fasta file. The fasta should contain all the sequences for genes in
#'   genelist_stabilityInfo tibble. Sequences can be subsetted from all genes
#'   sequences using fastaR::fa_some_records() function or genome-wide sequences
#'   can be provided.
#' @return a tibble of optimal codon percentage for each gene, boxplot of codon
#'   distribution in different class or stability of genes.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  codon_table <- system.file("extdata", "codonOptimalityTable_Yeasts.txt", package = "FungalTranscriptStability")
#'  org_codon_table <- readr::read_tsv(codon_table) %>% dplyr::select(Codon, Cgla)
#'  ref_sequence_file <- system.file("extdata", "CgROS_stability.fa", package = "FungalTranscriptStability")
#'  }
#' }
#' @rdname optimal_codon_distribution
#' @export
#' @importFrom broom tidy
#' @importFrom dplyr mutate left_join group_by filter select
#' @importFrom purrr map map2 map_dbl
#' @importFrom seqinr read.fasta uco getSequence
#' @importFrom tibble as_tibble
#' @importFrom tidyr gather unnest nest spread
#' @importFrom ggplot2 ggplot aes facet_wrap geom_boxplot scale_fill_manual
#'   scale_y_log10 theme_bw labs scale_x_discrete theme element_blank
#'   element_text element_rect
optimal_codon_distribution <- function(optimal_codon_table,genelist_stabilityInfo, ref_sequence_file){

  colnames(genelist_stabilityInfo) <- c("gene_id", "stability", "condition") # rename columns of stability info table

  colnames(optimal_codon_table) <- c("codon", "organism") # rename columns of codon table

  seq_dat <- broom::tidy(ref_sequence_file)

  codon_calci <- seq_dat %>%
    dplyr::mutate(cond_freq=purrr::map(x, function(ii){
      seqinr::read.fasta(ii,as.string = TRUE,forceDNAtolower = FALSE) %>% # read input fasta file
        purrr::map(~seqinr::uco(seqinr::getSequence(.x, as.string = FALSE),frame = 0,index="eff")) %>% # separate the sequences as codons and calculate freq of each codon
        do.call("cbind",.)  %>%
        tibble::as_tibble(rownames = "codon") %>%
        dplyr::mutate(codon=toupper(codon)) %>% # convert lower case to upper to match codon table
        dplyr::left_join(optimal_codon_table, by = c("codon"))  %>% # combine codon table with codon freq
        tidyr::gather("gene_id","freq", -codon, -organism) %>%
        dplyr::group_by(gene_id) %>%
        dplyr::mutate(total = sum(freq)) %>% # calculate total N and O frequencies
        dplyr::filter(organism %in% "O") %>%
        dplyr::mutate(optimal=sum(freq), percent=100*(optimal/total)) })) %>% # calculate percentage of optimal codons
    dplyr::select(-c(x)) %>%
    tidyr::unnest(cols = c(cond_freq)) %>%
    dplyr::select(c(gene_id, total, optimal, percent)) %>%
    unique()

  combine_stability <- dplyr::left_join(genelist_stabilityInfo, codon_calci) %>%
                          dplyr::mutate(stability=forcats::as_factor(stability))

  if(any(is.na(combine_stability))){

    print("no. of genes in stabilityInfo", nrow(genelist_stabilityInfo))
    print("no. of sequence", nrow(codon_calci))
    warning("sequence file and stability does not match!")

  }

  gp <- combine_stability %>%
    ggplot2::ggplot(ggplot2::aes(stability, percent, fill=stability)) +
    ggplot2::facet_wrap(~condition)+
    ggplot2::geom_boxplot(notch = TRUE, notchwidth = 0.4, alpha=0.8, width=0.5) +
    ggplot2::scale_fill_manual(values=c("#FC4E07", "#00AFBB" ))+
    ggplot2::scale_y_log10()+
    ggplot2::theme_bw()+
    ggplot2::labs(x="", y="Optimal codons (in %)")+
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

  print(gp)

  # compute p-value

  var1 <- unique(combine_stability$stability %>% levels() %>% .[1])
  var2 <- unique(combine_stability$stability %>% levels() %>% .[2])

  var1 <- rlang::sym(var1)
  var2 <- rlang::sym(var2)

  t_test <- combine_stability %>%
    dplyr::select(c(stability, condition, percent)) %>%
    dplyr::group_by(stability, condition) %>%
    tidyr::nest() %>%
    tidyr::spread(key = stability, value = data) %>%
    dplyr::mutate(
      t_test = purrr::map2(!!var1, !!var2, ~ t.test(.x$percent, .y$percent) %>% broom::tidy())) %>%
    dplyr::mutate(pval = purrr::map_dbl(t_test , ~ .$`p.value`)) %>%
    dplyr::select(condition, pval)

  print(t_test)

  return(combine_stability)

}
