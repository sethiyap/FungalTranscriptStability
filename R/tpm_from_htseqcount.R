#' @title tpm_from_htseqcount
#' @description TPM can be calculated from count matrix. Count matrix can be
#'   output of htseq-count or prepDE in Hisat2-Stringtie pipeline.
#' @param dir string, path of the text files output from htseq-count
#' @param pattern string, determining pattern to detect the count file and to
#'   remove extensions from the sample names
#' @param gff_file string, path or connection to gff file
#' @param HTSeqOutput logical, whether the output is Htseq-count or not. If it's
#'   htseq-count output last sixlines are filtered from further process.
#'   Default: TRUE
#' @param drop_genes string, to drop genes with specific pattern like
#'   mitochondrial, Default: NULL
#' @param minimum_read_cutoff numeric, filter genes which has less than cut-off
#'   reads in all of the samples , Default: 2
#' @param write_output logical, whether to return output on console or write in
#'   a file. If true, output file is generated. Default: FALSE
#' @param output_name string, name or path of output file name.
#' @return  a correlation plot and tpm matrix of input files.
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  gff_file <- "/Users/Pooja/Documents/PhD/thesis_figures/Chapter3/C_glabrata_CBS138_version_s02-m07-r38_features.gff"
#'  pattern="_12.*_HTSeqCount.txt"
#'  dir="/Users/Pooja/Documents/PhD/thesis_figures/Chapter3/rna"
#'  tpm_from_htseqcount(dir = dir, pattern = pattern,gff_file=gff_file, HTSeqOutput = TRUE, write_output = TRUE, output_name = "tpm_ros_rna_repeat1",drop_genes = "CaglfM*" )
#'  }
#' }
#' @rdname tpm_from_htseqcount
#' @export
#' @importFrom tibble tibble column_to_rownames
#' @importFrom dplyr mutate slice filter rename all_of filter_at vars any_vars
#'   inner_join select group_by
#' @importFrom purrr map
#' @importFrom data.table fread
#' @importFrom tidyr unnest spread gather
#' @importFrom stringr str_detect
#' @importFrom GenomicFeatures makeTxDbFromGFF transcriptLengths
#' @importFrom GGally ggpairs wrap
#' @importFrom ggplot2 theme_bw theme element_text element_rect
#' @importFrom readr write_delim

tpm_from_htseqcount <- function(dir, pattern, gff_file, HTSeqOutput=TRUE, drop_genes=NULL, minimum_read_cutoff=2, write_output=FALSE, output_name){

  #--- Load files
  dd <- list.files(path = dir,
                   pattern=pattern,
                   full.names = T)

  if(length(dd)==0L){
    warning("file name or pattern not correct")
  }else{
    print(dd)
  }

  #--- get count matrix
  df_count_matrix <- tibble::tibble(file_name = dd) %>%
                              dplyr::mutate(file_cont = purrr::map(file_name,data.table::fread,data.table = F))  %>%
                              tidyr::unnest(cols = c(file_cont), names_repair="universal") %>%
                              dplyr::mutate(file_name = basename(gsub(pattern=pattern,replacement="",file_name))) %>%
                              tidyr::spread(key = file_name , value = V2)

  print(colnames(df_count_matrix))

       if(HTSeqOutput==TRUE){
                  df_count_matrix <- df_count_matrix %>%
                  dplyr::slice(6:nrow(.))#remove first fivelines of statistics
      } else{
                  df_count_matrix=df_count_matrix
      }

      if(is.null(drop_genes)){
                  df_count_matrix <- df_count_matrix
      }
      else{
                  df_count_matrix <- df_count_matrix  %>% # remove mitochondrial mRNA
                                          dplyr::filter(!stringr::str_detect(V1, drop_genes))
      }
          column1 <- df_count_matrix %>% colnames() %>% .[1]
          df_count_matrix <- df_count_matrix %>%
                             dplyr::rename(gene_id= dplyr::all_of(column1))

        df_count_matrix_filtered <- df_count_matrix %>%
                                    dplyr::filter_at(dplyr::vars(-gene_id),
                                    dplyr::any_vars(. >= minimum_read_cutoff))


        gff <- GenomicFeatures::makeTxDbFromGFF(gff_file, metadata = T, format="gff")
        trans_len <- GenomicFeatures::transcriptLengths(gff, with.cds_len = T)

        tpm_matrix <- df_count_matrix_filtered %>%
                      dplyr::inner_join(trans_len) %>%
                      dplyr::select(c(colnames(df_count_matrix_filtered),cds_len)) %>%
                      tidyr::gather(names,rcount,-cds_len,-gene_id) %>%
                      dplyr::group_by(names) %>%
                      dplyr::mutate(rpk=(rcount/cds_len), pm=sum(rpk)/(10^6), tpm=round((rpk/pm),2)) %>%
                      dplyr::select(-c(rpk, pm, cds_len, rcount)) %>%
                      tidyr::spread(key=names, value=tpm)

        tpm_matrix_1 <- tpm_matrix %>% tibble::column_to_rownames("gene_id")

        gg <- GGally::ggpairs(log2(tpm_matrix_1+0.01),
                              upper = list(continuous = GGally::wrap("cor", size = 5, color="black"))) +
          ggplot2::theme_bw() +
          ggplot2::theme(legend.text = ggplot2::element_text(size = 18),
                          legend.title = ggplot2::element_text(size = 20),
                          strip.text = ggplot2::element_text(size = 8,color="black"),
                          axis.title = ggplot2::element_text(size = 10,color="black"),
                          strip.background = ggplot2::element_rect(color="white"),
                          axis.text = ggplot2::element_text(size = 10,color="black"))
      print(gg)

        if(write_output==TRUE){
          readr::write_delim(tpm_matrix, paste(output_name,"_tpm_matrix.tab",sep=""),delim="\t",col_names = TRUE)
        }else{
          return(tpm_matrix)
        }
}


