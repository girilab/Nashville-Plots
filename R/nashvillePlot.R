require('ggplot2') #these are not supposed to be used for package
#require('ggrepel')

#' Validate Input
#'
#' This function internally renames the column names for Chromosome, Base Pair
#' Location, P-value, and SNP name. Cases where Chromosome, Base Pair Location
#' or P-value are not numeric values are dropped.
#' @param x data frame of SNPs
#' @param chr_name column name of chromosome number
#' @param pos_name column name of SNP position along a chromosome given as number of base pairs
#' @param snp_name column name containing SNP names
validate_input <- function(x, chr_name="CHR", pos_name="BP", p_val_name="P", snp_name = "SNP", ...)
{
  if ((chr_name %in% names(x)) & (pos_name %in% names(x)) & (p_val_name %in% names(x))) {
    if (is.null(x[[snp_name]])) {
      colnames(x)[colnames(x) %in% c(chr_name)] <- c("CHR")
      colnames(x)[colnames(x) %in% c(pos_name)] <- c("BP")
      colnames(x)[colnames(x) %in% c(p_val_name)] <- c("P")
    } else {
      colnames(x)[colnames(x) %in% c(chr_name)] <- c("CHR")
      colnames(x)[colnames(x) %in% c(pos_name)] <- c("BP")
      colnames(x)[colnames(x) %in% c(p_val_name)] <- c("P")
      colnames(x)[colnames(x) %in% c(snp_name)] <- c("SNP")
    }
  }
  else {
    stop("Please double check columns to ensure column names for chromosome, BP and P-value match")
  }
  gwas <- subset(x, is.numeric(CHR) & is.numeric(BP) & is.numeric(P))
  gwas <- gwas[order(gwas$CHR, gwas$BP), ]
  gwas$log_p = log10(gwas$P)
  gwas
}

#' Generate  X axis
#'
#' This function creates an absolute order of all SNPs across all chromosomes
#' and optionally trims the number of points that will be graphed by removing
#' points that are unlikely to be visible.
#'
#' @param x data frame of SNPs
#' @param
#new_xaxis_func_v2 <- function(x, chr_name="CHR", pos_name="BP", p_val_name="P", snp_name = "SNP", ...){
#  gwas <- validate_input(x, chr_name, pos_name, p_val_name, snp_name, ...)
#  gwas <- gwas[order(gwas$CHR, gwas$BP), ]
#  gwas$log_p = log10(gwas$P)
#  gwas
#}


#' Read MetaXcan Folder
#'
#' This function reads in metaXcan or prediXcan results files by tissue
#' @param directory a string representation of a file path containing predixcan output
#' @param map a file path containing a mapping gene names to ensembl ID. The header should be "Gene ENSG CHR START_POS END_POS"
#' @param label a list of labels to use instead of file names for graphing. These should probably be in alphabetical order.
#' @param pattern a regex describing files to be read from "directory"
#' @export
read_metaXcan_folder <-  function(directory, map, label=c(), pattern='*.csv$') {
  files <- list.files(directory, pattern=pattern)
  appended_gene_d <- read.table(map, header=TRUE)
  tissues <- data.frame()
  if(length(label)!=0){
    stopifnot(length(label)==length(files))
  }
  for (i in 1:length(files)){
    tissue_file <- read.csv(file=as.character(paste0(directory, files[i])), header = T)
    tissue_file$tiss_number <- i
    tissue_file$dtype <- ifelse(length(label)==length(files), label[i], files[i])
    tissues <- rbind(tissues, tissue_file)
  }
  appended_gene_d <- merge(tissues, appended_gene_d, by.x = "gene", by.y = "ENSG")
  names(appended_gene_d)[names(appended_gene_d) == "gene"] <- "ENSG"
  t1 <- which(duplicated(appended_gene_d$ENSG))
  if(length(t1) == 0) {appended_gene_d <- appended_gene_d} else {appended_gene_d <- appended_gene_d[-c(t1), ]}
  appended_gene_d$MID_POS <- (appended_gene_d$END_POS + appended_gene_d$START_POS)/2
  appended_gene_d$log_p <- -log10(appended_gene_d$pvalue)
  appended_gene_d
}

#' Generate the nashville plot
#'
#' This function returns a ggplot2 graph of the nashville plot
#'
#' @param gwas data frame of gwas summary stats
#' @param metaxcan data frame from read_metaxcan_folder
#' @param left numeric, leftmost point to plot
#' @param right numeric, rightmost point to plot
#' @param y_min numeric, ggplot ymin aes
#' @param y_max ggplot ymax aes
#' @param y_ticks numeric, number of y breaks
#' @param test_ylab string, y axis label
#' @param x_axis_name string, x axis label
#' @param gene_tag_p numeric, annotate the SNPs with P-values more extreme than log(gene_tag_p)
#' @param color_tissue list of colors to label tissues with
#' @param sig_line1 numeric, draw a horizontal line at log(sig_line1)
#' @param sig_line2 numeric, draw a horizontal line at log(sig_line2)
#' @param sig_line1_color color for sig_line1
#' @param sig_line2_color color for sig_line2
#' @param labels_cat list, tissue labels
#' @param xbreaks numeric, number of x breaks
#' @param xlabels list of labels for genes
#' @param draw_genes draw lines along the length of each gene
#' @param draw_bottom bool whether or not to plot data from gwas
#' @param draw_top bool whether or not to plot data from metaxcan
#' @importFrom ggplot2 ggplot aes theme_bw guides geom_label_repel geom_hline guide_legend scale_x_continuous scale_y_continuous scale_colour_manual theme element_text element_line element_blank xlab ylab expand_limits geom_hline geom_point
plot_whole_genome <- function(gwas=gwas, metaxcan=metaxcan, y_min = NULL, y_max = NULL, y_ticks = NULL, test_ylab = test_ylab, x_axis_name, gene_tag_p = NULL,
                              color_tissue = NULL, sig_line1 = NULL, sig_line2 = NULL, sig_line1_color = NULL, sig_line2_color = NULL, labels_cat,
                              xbreaks, xlabels, left, right, draw_genes, draw_bottom, draw_top)
{
  if(is.null(y_min) == TRUE) {y_min <- round(min(gwas$log_p, na.rm=TRUE) - 1)} else {y_min <- y_min}
  if(is.null(y_max) == TRUE) {y_max <- round(max(metaxcan$log_p, na.rm=TRUE) + 1)} else {y_max <- y_max}
  #y_ticks <- if(is.null(y_ticks) == TRUE && y_max - y_min > 16) 16 else round(y_max - y_min)
  if(is.null(y_ticks) == TRUE) {
    if((y_max - y_min) > 16) {
      y_ticks <- 16
    } else {y_ticks <- round(y_max - y_min)}
  } else {y_ticks <- y_ticks}

  break_length <- round(round(y_max - y_min)/y_ticks)
  for_tag <- metaxcan[metaxcan$log_p > gene_tag_p, ]                    #make sure we know that log_p & Gene exist!
  for_tag <- for_tag[order(for_tag$Gene, -for_tag$log_p), ]
  for_tag <- for_tag[!duplicated(for_tag$Gene),]
  if(is.null(sig_line1) == TRUE) {sig_line1 <- 0} else{sig_line1 <- sig_line1}
  if(is.null(sig_line2) == TRUE) {sig_line2 <- 0} else{sig_line2 <- sig_line2}
  if(is.null(sig_line1_color) == TRUE) {sig_line1_color <- "black"} else {sig_line1_color <- sig_line1_color}
  if(is.null(sig_line2_color) == TRUE) {sig_line2_color <- "black"} else {sig_line2_color <- sig_line2_color}
  p1 <- ggplot() +
    theme_bw() +
    #scale_x_continuous(breaks=xbreaks, labels=xlabels, limits=c(left, right)) + #, limits=c(left, right)) +
    scale_x_continuous(breaks=xbreaks, labels=xlabels) +
    theme(axis.text.x=element_text(size=12, color='black'),
          axis.text.y=element_text(size=12, color='black'),
          axis.title.x=element_text(size = 12, face="bold", color="black"),
          axis.title.y=element_text(size = 12, face="bold", color="black"),
          axis.ticks.x=element_line()) +
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank()) +
    xlab(x_axis_name) +
    ylab(test_ylab) +
    scale_y_continuous(breaks=seq(y_min, y_max, break_length)) +
    expand_limits(y=c(y_min, y_max)) +
    geom_hline(aes(yintercept = 0), size = 1)

  if(draw_bottom) {
    if(!draw_top) {
      p1 <- p1 + geom_point(data=gwas, aes(x=absolutePos, y=-log_p, size=1.0, color=factor(col_cat2)))
    } else {
      p1 <- p1 + geom_point(data=gwas, aes(x=absolutePos, y=log_p, size=1.0, color=factor(col_cat2)))
    }
    p1 <- p1 + scale_colour_manual(name = "Tissue Type",
                                   values = c("darkgray", "black", color_tissue),
                                   labels = c("Single SNP", "Single SNP", labels_cat)) +
    guides(shape = "none", size = "none", colour = guide_legend(reverse = TRUE, override.aes = list(size=6))) +
    ggrepel::geom_label_repel(data = for_tag, aes(x = absolutePos, y = log_p, label = Gene)) +
    geom_hline(yintercept = sig_line1, color = c(sig_line1_color), size = 1) +
    geom_hline(yintercept = sig_line2, color = c(sig_line2_color), size = 1)
  }
  if(draw_top) {
    if(draw_genes) {
      temp_data <- metaxcan
      metaxcan$abs <- metaxcan$absolutePos - (metaxcan$MID_POS - metaxcan$START_POS )   #start_pos absolute
      temp_data$abs <-  metaxcan$absolutePos + (metaxcan$END_POS - metaxcan$MID_POS)    #end_pos absolute
      metaxcan <- rbind(metaxcan, temp_data)
      p1 <- p1 + geom_line(data = metaxcan, aes(x = abs, y = log_p, size = 1.0, color = dtype, group = interaction(Gene, dtype)))
    } else {
      p1 <- p1 + geom_point(data = metaxcan, aes(x = absolutePos, y = log_p, size = 1.0, color = dtype, shape = factor(col_cat2)))
    }
    p1 <- p1 + guides(shape = "none", size = "none")
  }
  if(!draw_bottom | !draw_top) {
    #p1 <- p1 + ylim(0, NA)
  }
  #p1 <- p1 + guides(size = FALSE, shape = FALSE)
  p1
}

#' Get gene bounds
#'
#' Determines boundaries for plotting a gene by gene name
#'
#' @param gene_name name of a gene
#' @param map_df data frame mapping gene names to ensembl ID.
get_gene_bounds <- function(gene_name, map_df) {
  lower <- unique(map_df[which(map_df$Gene==gene_name), "START_POS"])
  upper <- unique(map_df[which(map_df$Gene==gene_name), "END_POS"])
  if (length(lower) > 1) { stop(paste0("Multiple genes with name: ", gene_name, " Try giving an ENSG."))}
  lower <- lower - (upper - lower)* 0.05
  upper <- upper + (upper - lower)* 0.05
  bounds <- c(lower, upper)
}

#' Get gene bounds
#'
#' Determines boundaries for plotting a gene  by ensembl ID
#'
#' @param ensg ensembl ID
#' @param map_df data frame mapping gene names to ensembl ID.
get_gene_bounds_ensg <- function(ensg, map_df) {
  lower <- unique(map_df[which(map_df$ENSG==ensg), "START_POS"])
  upper <- unique(map_df[which(map_df$ENSG==ensg), "END_POS"])
  if (length(lower) > 1) { stop(paste0("Multiple genes with name: ", ensg, " Try giving a gene name"))}
  lower <- lower - (upper - lower) * 0.05
  upper <- upper + (upper - lower) * 0.05
  bounds <- c(lower, upper)
}


#' Create a Nashville Plot
#'
#' This function makes two plots of GWAS summary data which share the same
#' x-axis, one pointed up and the other pointed down.
#'
#' @param gwas data frame of gwas summary stats
#' @param metaxcan data frame from read_metaxcan_folder
#' @param left numeric, leftmost point to plot
#' @param right numeric, rightmost point to plot
#' @param y_min numeric, ggplot ymin aes
#' @param y_max ggplot ymax aes
#' @param y_ticks numeric, number of y breaks
#' @param test_ylab string, y axis label
#' @param x_axis_name string, x axis label
#' @param gene_tag_p numeric, annotate the SNPs with P-values more extreme than log(gene_tag_p)
#' @param color_tissue list of colors to label tissues with
#' @param sig_line1 numeric, draw a horizontal line at log(sig_line1)
#' @param sig_line2 numeric, draw a horizontal line at log(sig_line2)
#' @param sig_line1_color color for sig_line1
#' @param sig_line2_color color for sig_line2
#' @param labels_cat list, tissue labels
#' @param xbreaks numeric, number of x breaks
#' @param xlabels list of labels for genes
#' @param draw_genes draw lines along the length of each gene
#' @param draw_bottom bool whether or not to plot data from gwas
#' @param draw_top bool whether or not to plot data from metaxcan
#' @param chr Limit the graph to the given chromosome
#' @param zoom_left plot genes within a gene after this base pair
#' @param zoom_right plot genes within a gene before this base pair
#' @param sample draw what percent of points with P-value greater than 0.1
#' @export
nashville.plot <- function(gwas=gwas, metaxcan=NULL, map_df, chr = NULL, y_min = NULL, y_max = NULL, y_ticks = NULL, zoom_left = 0,
                        zoom_right = 0, zoom_gene = NULL, zoom_ensg = NULL, gene_tag_p = NULL, color_tissue = NULL, sig_line1 = NULL, sig_line2 = NULL,
                        sig_line1_color = NULL, sig_line2_color = NULL, draw_bottom = TRUE, draw_top = TRUE, samp = NULL)
{
  labels_cat <- c(sort(unique(as.character(metaxcan$dtype))))
  test_ylab <- expression(log["10"]*italic((p))~"&"~-log["10"]*italic((p)))
  gwas <- validate_input(gwas, chr_name = "CHR", pos_name = "BP", p_val_name = "P", snp_name = "SNP")
  if(!is.null(samp)) {
    greater <- metaxcan[ sample(which(metaxcan$pvalue>0.1), round(samp*length(which(metaxcan$pvalue>0.1)))), ]
    lesser <- metaxcan[which(metaxcan$pvalue<=0.1), ]
    metaxcan <- rbind(lesser, greater)
  }
  gwas$col_cat2 <- as.numeric(gwas$CHR) %% 2
  metaxcan$col_cat2 <- as.numeric(metaxcan$CHR) %% 2

  min_pos <- tapply(gwas$BP, gwas$CHR, min)                                             #make sure we know that BP & CHR exist!
  max_pos <- tapply(gwas$BP, gwas$CHR, max)
  chr_shift <- head(c(0,cumsum(as.numeric(max_pos))),-1)
  if (all(chr_shift==0)) {
    chr_shift <- rep(0, times=max(rbind(gwas$CHR, max(as.numeric(metaxcan$CHR)))))
  }
  gwas$absolutePos <- gwas$BP + chr_shift[gwas$CHR]
  metaxcan <- metaxcan[which(metaxcan$CHR %in% row.names(min_pos)), ]
  metaxcan$absolutePos <- metaxcan$MID_POS + chr_shift[as.numeric(metaxcan$CHR)]    #make sure we know that MID_POS exists!

  left <- NA
  right <- NA
  draw_genes <- F
  if(!is.null(chr)) {
    message('Chromosome specific plot\n')
    x_axis_name <- paste0("Chromosome ",  chr, " (MB)")
    if (!is.null(zoom_gene)) {
      message('Zoom to gene by name\n')
      bounds <- get_gene_bounds(zoom_gene, map_df)
      left <- chr_shift[chr] + bounds[1]
      right <- chr_shift[chr] + bounds[2]
      draw_genes <- T
    }
    else if (!is.null(zoom_ensg)) {
      message('zoom to gene by ENSG\n')
      bounds <- get_gene_bounds_ensg(zoom_ensg, map_df)
      left <- chr_shift[chr] + bounds[1]
      right <- chr_shift[chr] + bounds[2]
      draw_genes <- T
    }
    else if ( zoom_left !=0 | zoom_right !=0) {
      message('Manual zoom\n')
      left <- chr_shift[chr] + zoom_left
      right <- chr_shift[chr] + zoom_right
      draw_genes <- T
    } else {  #zoom about a chromosome
      left <- chr_shift[chr]
      if (chr_shift[chr+1] != 0) {
        right <- chr_shift[chr+1]
      }
    }

    #pad by 2 kilobase
    left <- left - 2000
    right <- right + 2000

    gwas <- gwas[which(gwas$absolutePos >= left & gwas$absolutePos <= right), ]                   #send plotted data only
    metaxcan <- metaxcan[which(metaxcan$absolutePos >= left & metaxcan$absolutePos <= right), ]             #send plotted data only
    #convert to relative position in a chromosome
    gwas$absolutePos <- gwas$absolutePos - left
    metaxcan$absolutePos <- metaxcan$absolutePos - left

    xbreaks <- waiver()
    xlabels <- waiver()
    } else {
    message('Whole genome plot\n')
    x_axis_name <- "Chromosome"
    xbreaks <- (chr_shift + (max_pos - min_pos)/2) * 1e-6                                 #megabase
    xlabels <- unique(gwas$CHR)
  }
  if(is.null(color_tissue) == TRUE) {
    color_tissue <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(labels_cat))
  }
  gwas$absolutePos <- gwas$absolutePos * 1e-6                                                 #megabase
  metaxcan$absolutePos <- metaxcan$absolutePos * 1e-6                                     #megabase
  metaxcan$START_POS <- metaxcan$START_POS * 1e-6
  metaxcan$END_POS <- metaxcan$END_POS * 1e-6
  metaxcan$MID_POS <- metaxcan$MID_POS * 1e-6
  left <- left * 1e-6                                                                     #megabase
  right <- right * 1e-6                                                                   #megabase

  p1 <- plot_whole_genome(gwas, metaxcan, y_min, y_max, y_ticks, test_ylab, x_axis_name, gene_tag_p, color_tissue, sig_line1, sig_line2, sig_line1_color, sig_line2_color, labels_cat, xbreaks, xlabels, left, right, draw_genes, draw_bottom, draw_top)
}



