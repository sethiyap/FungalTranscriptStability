#' @title deseq_from_htseqcount
#' @description Computes normalized gene expression value and differential gene
#'   expression using DESeq2 from read counts matrix obtained from HTSeqCount or
#'   any other method.
#' @param dir string, absolute path of directory containing read count file for
#'   each sample (a HTSeqCount output)
#' @param pattern string, pattern of files to be selected as input from the dir.
#' For example `pattern="*._HTSeqCount.txt`
#'   The pattern is also used to replace the pattern name from file, which is
#'   used as column names in further processing. The pattern replaced column
#'   names should match with first column of metadata file.
#' @param metadata_file string, absolute path of metadata file in tab-delimited
#'   format. metadata file should contain two columns, First: names of the
#'   samples exactly matching with the pattern replaced file names. Second:
#'   sample type i.e. treated or untreated
#'sample1_set1	treated
# sample1_set2	treated
# sample2_set1	untreated
# sample2_set2	untreated
#' @param HTSeqOutput logical, whether the read count files are HTSeq-count
#'   output or not. If TRUE, last six rows from containing the statistics are
#'   ignored. If FALSE, the file is considered as it is. Default: TRUE
#' @param drop_genes string, pattern of the genes to be dropped (Like
#'   mitochondrial rRNA, rRNA, tRNAs) from calculation.For example:'CaglfM.*'
#'   mitochrondrial rRNA in C. glabrata. Default:NULL
#' @param minimum_read_cutoff numeric, remove genes containing reads less than
#'   the threshold in all samples from further processing, Default: 2
#' @param write_output logical, whether to return output on console or write in
#'   a file. If true, output file is generated. Default: FALSE
#' @param outfile string, name or path of output file name.
#' @return combined count matrix of all samples, log2 correlation plot and
#'   expression matrix
#' @details DETAILS
#' @rdname deseq_from_htseqcount
#' @export
#' @importFrom tibble tibble
#' @importFrom dplyr mutate slice filter filter_at vars any_vars rename all_of select_all tibble select
#' @importFrom purrr map
#' @importFrom data.table fread
#' @importFrom tidyr unnest spread
#' @importFrom stringr str_detect
#' @importFrom readr write_delim read_delim
#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateSizeFactors DESeq results counts
#' @importFrom GGally ggpairs wrap ggcorr
#' @importFrom ggplot2 theme_bw theme element_text element_blank
#' @importFrom writexl write_xlsx
#' @examples
#' \dontrun{
#' if(interactive()){
#'  dir <- system.file(file.path('extdata'), package='FungalTranscriptStability')
#'  pattern="_S.*_HTSeqCount.txt"
#'  metadata <- system.file(file.path('extdata', 'metadata.txt'), package='FungalTranscriptStability')
#'  deseq_from_htseqcount(dir = dir, pattern = pattern, metadata_file = metadata_file, drop_genes = "CaglfM*", outfile = "ROS")
#'  }
#' }
#'
deseq_from_htseqcount <- function(dir, pattern, metadata_file,HTSeqOutput=TRUE,drop_genes=NULL, minimum_read_cutoff=2,write_output=FALSE, outfile) {

          #--- Load DESeq2 files
          dd <- list.files(path = dir,
                           pattern=pattern,
                           full.names = T)
          print(dd)

          #--- get count matrix
          df_count_matrix <- tibble::tibble(file_name = dd) %>%
                    dplyr::mutate(file_cont = purrr::map(file_name,data.table::fread,data.table = F))  %>%
                    tidyr::unnest(cols = c(file_cont), names_repair="universal") %>%
                    dplyr::mutate(file_name =gsub(pattern=pattern,replacement="", basename(file_name)))  %>%
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
                    df_count_matrix=df_count_matrix  %>%
                                        dplyr::filter(!stringr::str_detect(V1, drop_genes))
          }



          #--- filter genes if read count is less than 10 in all columns
          df_count_matrix_filtered <- df_count_matrix %>%
                                          dplyr::filter_at(dplyr::vars(-V1),
                                                           dplyr::any_vars(. >= minimum_read_cutoff))

          df_count_matrix_filtered <-  data.frame(df_count_matrix_filtered)
          row.names(df_count_matrix_filtered) <- df_count_matrix_filtered$V1
          df_count_matrix_filtered <- subset(df_count_matrix_filtered, select=-c(V1))

          #--- add conditions
          condition <- readr::read_delim(metadata_file, delim="\t", col_names = FALSE)

          columnname <- condition %>% colnames()

          condition <-  condition %>%
                              dplyr::rename(sample= dplyr::all_of(columnname[1]),
                                        type=dplyr::all_of(columnname[2]))

          df_count_matrix_filtered_1 <- df_count_matrix_filtered %>%
                    #dplyr::select_( .dots =condition$X1[1:2] ) %>% # for no replicates
                     dplyr::select_all( .dots =condition$sample )

            rownames(df_count_matrix_filtered_1) <-    row.names(df_count_matrix_filtered)
          colData <- as.data.frame(colnames(df_count_matrix_filtered_1))
          colData$condition <- condition$type[match(colData[,1], condition$sample)]

          colnames(colData) <- c("colData", "condition")
          print(ncol(df_count_matrix_filtered_1) == nrow(colData))

          colData$condition=factor(colData$condition, levels=c("untreated", "treated"))
          print(colData)

          #--- make matrix for deseq2
          dds <- DESeq2::DESeqDataSetFromMatrix(countData = df_count_matrix_filtered_1,
                                        colData = colData,design = ~ condition)

          dds <- DESeq2::estimateSizeFactors(dds) # run this or DESeq() first

          dds$condition <- factor(dds$condition)

          #--- Compute DESeq2
          dds  <- DESeq2::DESeq(dds) # for replicates

          #--- For samples with no replicates
          # dds <- estimateDispersionsGeneEst(dds)
          # dispersions(dds) <- mcols(dds)$dispGeneEst
          # dds <- nbinomWaldTest(dds)

          #--- get DEG
          res <- DESeq2::results(dds)
          cat(summary(res))

          resdata <- merge(as.data.frame(res),
                           as.data.frame(DESeq2::counts(dds, normalized=TRUE)),
                           by="row.names", sort=FALSE)

          names(resdata)[1] <- "Gene"
          print(head(resdata))

          deseq_matrix <- resdata %>%
                                        dplyr::tibble() %>%
                                        dplyr::select(condition$sample)


          #--- get correlation matrix
          gg <-
                    GGally::ggpairs(log2(deseq_matrix+1), upper = list(continuous = GGally::wrap("cor", size = 5, color="black"))) +
                    ggplot2::theme_bw() +
                    ggplot2::theme(legend.text = ggplot2::element_text(size = 25),
                          legend.title = ggplot2::element_text(size = 20),
                          axis.title.x = ggplot2::element_text(size = 15),
                          axis.title.y = ggplot2::element_text(size = 15),
                          axis.text.x = ggplot2::element_text(size = 10),
                          axis.text.y = ggplot2::element_text(size = 10),
                          panel.grid = ggplot2::element_blank())


          #--- corrplot
          rcorr = cor(as.matrix(log2(deseq_matrix+1)))
          gc <- GGally::ggcorr(log2(deseq_matrix+1),label=TRUE,label_size = 8, hjust = 0.75,
                               size = 5,
                               label_color = "white",
                               label_round = 2,
                               layout.exp = 1)



          print(gc)


          if(write_output==TRUE){

            #--- write count matrix in a file
            readr::write_delim(df_count_matrix, paste(outfile, "_count_matrix.tab", sep=""),col_names = TRUE,delim="\t")

            # correlation plot
            pdf(paste(outfile, "_correlation2_plot.pdf", sep=""),width = 12, height = 12)
            print(gc)
            dev.off()

            png(paste(outfile, "_correlation1_plot.png", sep=""),width = nrow(condition)*100, height = nrow(condition)*100, units = "px", pointsize = 12)
            print(gg )
            dev.off()

            # DEG output
            up_deg <- subset(resdata, resdata$log2FoldChange> 0.6 & pvalue <= 0.05)
            down_deg <- subset(resdata, resdata$log2FoldChange< -0.6 & pvalue <= 0.05)
            writexl::write_xlsx(list(deseq_out=resdata,up=up_deg, down=down_deg),path = paste(outfile, "_deseq_output.xlsx", sep=""),col_names = T)

            }else{
            return(resdata)
          }

}

