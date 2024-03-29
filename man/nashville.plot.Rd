% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nashvillePlot.R
\name{nashville.plot}
\alias{nashville.plot}
\title{Generate a Nashville plot}
\usage{
nashville.plot(
  data1,
  data2 = NULL,
  data1_direction = 1,
  map_df = "37",
  chr = NULL,
  zoom_ensg = NULL,
  zoom_gene = NULL,
  zoom_left = 0,
  zoom_right = Inf,
  draw_genes = FALSE,
  labels_cat = NULL,
  group_color = NULL,
  sig_line1 = NULL,
  sig_line2 = NULL,
  sig_line1_color = "black",
  sig_line2_color = "black",
  gene_tag_1 = -Inf,
  gene_tag_2 = -Inf,
  config = data.frame(V1 = NA, V2 = NA, V3 = NA),
  ...
)
}
\arguments{
\item{data1}{data to be plotted generated by read_metaXcan_folder or read_gwas_file.}

\item{data2}{data to be plotted generated by read_metaXcan_folder or read_gwas_file. Optional}

\item{data1_direction}{direction for data1 to be drawn, 1 is up and -1 is down}

\item{draw_genes}{draw lines along the length of each gene instead of a dot at the midpoint}

\item{sig_line1}{numeric, draw a horizontal line at -log(sig_line1)}

\item{sig_line2}{numeric, draw a horizontal line at -log(sig_line2)}

\item{sig_line1_color}{color for sig_line1}

\item{sig_line2_color}{color for sig_line2}

\item{gene_tag_1}{numeric, annotate the SNPs with P-values more extreme than log(gene_tag_p)}

\item{gene_tag_2}{numeric, annotate the SNPs with P-values more extreme than log(gene_tag_p)}

\item{map_df:}{37 or 38 depending on the Human Genome assembly reference number, degfault = 37}

\item{chr:}{chromosome number. If provided, only data from that chromosome will be graphed}

\item{zoom_ensg:}{if `chr` is set this will graph around the gene described by endembl ID}

\item{zoom_gene:}{if `chr` is set this willgraph around the gene described by name}

\item{zoom_left:}{if `chr` is set this will graph points to the right of this base pair number}

\item{zoom_right:}{if `chr` is set this will graph points to the left of this base pair number}
}
\description{
This function returns a ggplot2 graph object of the Nashville plot
}
