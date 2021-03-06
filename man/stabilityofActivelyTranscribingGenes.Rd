% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stabilityofActivelyTranscribingGenes.R
\name{stabilityofActivelyTranscribingGenes}
\alias{stabilityofActivelyTranscribingGenes}
\title{stabilityofActivelyTranscribingGenes}
\usage{
stabilityofActivelyTranscribingGenes(
  rna_pol2_tbl,
  percent_cutoff = 15,
  write_output = FALSE,
  output_name = "sample",
  binwidth = 0.4
)
}
\arguments{
\item{rna_pol2_tbl}{a tibble of normalized RNAP II and mRNA values. With
first column (character) of gene-id, second column (numeric) RNAP II values
and third column (numeric) mRNA values.}

\item{percent_cutoff}{numeric, cutoff to determine actively transcribing
genes, i.e. to consider top 10 or top 15percent genes. Default: 15}

\item{write_output}{logical, to write output in text file or not. Default:
FALSE}

\item{output_name}{character, if \code{write_output=TRUE} then determine the
name of output file. Default: 'sample'}
}
\value{
A tibble of stable and unstable genes determined from mRNA and RNAP
II values.
}
\description{
Stability of actively transcribed genes (i.e. belonging to top
10-15 percent of total genes by RNA polymerase II) is determined by comparing RNAP
II occupancy levels and steady state mRNA levels.
}
\details{
mRNA levels should be normalized to gene lengths as mRNA values of
genes within the given condition are considered to determine stability.
}
\examples{
\dontrun{
if(interactive()){
 rna_pol2_tbl <- readr::read_delim(pipe("pbpaste"), delim="\t", col_names = TRUE)
 }
}
}
