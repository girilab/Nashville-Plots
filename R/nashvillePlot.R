library('ggplot2')

#' Read config file
#'
#' Reads in configuration info from tab separated values file
#' @param file a file path to a 3 column tab separated values file without headers where the first column contains file names, the second column contains group names, and the 3rd column contains colors for those groups to be plotted in.
#' @export
read_config <- function(file) {
  config <- read.table(file, header=FALSE, stringsAsFactors = FALSE)
  if (length(unique(config$V2)) != length(config$V2)) {
    stop("Each element of the group column must be unique.")
  }
  config
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
  lower <- lower - (upper - lower) * 0.05
  upper <- upper + (upper - lower) * 0.05
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

#' Read GWAS File
#'
#' @param file path to gwas results as obtained from plink --assoc
#' @export
read_gwas_file <- function(file) {
  gwas <- read.table(file, header=T)
  gwas.obj <- make.valid.object(CHR = gwas$CHR, P = gwas$P, BP = as.numeric(gwas$BP), group=file)
}

#' Read MetaXcan Folder
#'
#' This function reads in metaXcan or prediXcan results files by tissue
#' @param directory a string representation of a file path containing predixcan output
#' @param map_df a file path containing a mapping gene names to ensembl ID. The header should be "Gene ENSG CHR START_POS END_POS"
#' @param label a list of labels to use instead of file names for graphing. These should probably be in alphabetical order.
#' @param pattern a regex describing files to be read from "directory"
#' @export
read_metaXcan_folder <-  function(directory, map_df = '37', pattern='*.csv$') {
  if (map_df == "37") {
    appended_gene_d <- gene.build.37
  } else if (map_df == "38") {
    appended_gene_d <- gene.build.38
  } else {
    stop("map_df must be '37' or '38'")
  }
  files <- list.files(directory, pattern=pattern)
  tissues <- data.frame()
  for (i in 1:length(files)){
    tissue_file <- read.csv(file=as.character(paste0(directory, files[i])), header = TRUE)
    tissue_file$dtype <- files[i]
    tissues <- rbind(tissues, tissue_file)
  }
  appended_gene_d <- merge(tissues, appended_gene_d, by.x = "gene", by.y = "ENSG")
  names(appended_gene_d)[names(appended_gene_d) == "gene"] <- "ENSG"
  # TODO: test functionality
  t1 <- which(duplicated(appended_gene_d$ENSG))
  if(length(t1) == 0) {appended_gene_d <- appended_gene_d} else {appended_gene_d <- appended_gene_d[-c(t1), ]}
  #appended_gene_d$log_p <- -log10(appended_gene_d$pvalue)
  appended_gene_d$CHR <- as.numeric(appended_gene_d$CHR)
  meta.obj <- make.valid.object(CHR = appended_gene_d$CHR,
                                P = appended_gene_d$pvalue,
                                gene.start = appended_gene_d$START_POS,
                                gene.end = appended_gene_d$END_POS,
                                group = appended_gene_d$dtype,
                                gene.name = appended_gene_d$Gene)
  meta.obj
}

#' Builds an object ready to plot
#'
#' This function take data.frame columns and returns an object ready to plot with nashville.plot()
#' @param CHR numeric column of chromosome numbers
#' @param P numeric column of P values to be plotted
#' @param BP numeric column of base pair locations within a chromosome, either this or gene.start and gene.end are required
#' @param gene.start numeric column of base pair locations for the start of a gene
#' @param gene.end numeric column of base pair locations for the end of a gene
#' @param group label for plotting, often this will be a tissue name
#' @param gene.name name of each gene, this will be shown for P values greater than the tag threshold
make.valid.object <- function(CHR, P, group=NA, BP=NA, gene.start=NA, gene.end=NA, gene.name=NA) {
  #Check types
  ifelse(all(is.numeric(CHR)), NA, stop("CHR must be numeric"))
  ifelse(all(is.numeric(P)), NA, stop("P must be numeric"))
  ifelse(all(is.numeric(BP)) | all(is.na(BP)), NA, stop("BP must be numeric"))
  ifelse(all(is.numeric(gene.start)) | all(is.na(gene.start)), NA, stop("gene.start must be numeric"))
  ifelse(all(is.numeric(gene.end)) | all(is.na(gene.end)), NA, stop("gene.end must be numeric"))
  #make sure that either BP or gene.start+gene.end exists
  if(all(is.na(gene.start)) && all(is.na(gene.end))) {
    stopifnot(is.numeric(BP))
  }
  if(all(is.na(BP))) {
    stopifnot(is.numeric(gene.start) && is.numeric(gene.end))
    BP <- (gene.start + gene.end)/2  # TODO this should be delivered from map_df
  }
  #set group if it wasn't provided
  if (missing(group)) {
    group <- ifelse(CHR %% 2 == 0, group, paste0(group, "\u200b"))
    #group <- paste0("Single SNP ", as.character(CHR %% 2))
  }
  object <- data.frame(gene.name, group, CHR, BP, P, gene.start, gene.end)
}

#' plot manhattan
#'
#' draws a manhattan plot onto the nashville plot
#' @param data data frame from make.valid.object()
#' @param direction direction for data to be drawn
#' @param draw_genes bool plot with geam_line
#' @importFrom ggplot2 ggplot aes geom_point
plot.mh <- function(data, direction, draw_genes) {
  if(draw_genes & any(!is.na(data$gene.name) & any(!is.na(data$group)))) {
    data_copy <- data
    data_copy$gene.start <- data_copy$gene.end
    data <- rbind(data, data_copy)
    mh <- geom_line(data = data, aes(x=gene.start,
                                     y=direction * log(P, 10),
                                     color=color,#color=as.factor(group),
                                     group=interaction(gene.name, group)))
  }
  else {
    mh <- geom_point(data=data, aes(x=absolute,
                                    y=direction * log(P, 10),
                                    color=color,#color=as.factor(group),
                                    shape=as.factor(CHR %% 2)))
  }
}

#' Generate a Nashville plot
#'
#' This function returns a ggplot2 graph object of the Nashville plot
#'
#' @param data1 data to be plotted generated by read_metaXcan_folder or read_gwas_file.
#' @param data2 data to be plotted generated by read_metaXcan_folder or read_gwas_file. Optional
#' @param map_df: 37 or 38 depending on the Human Genome assembly reference number, degfault = 37
#' @param chr: chromosome number. If provided, only data from that chromosome will be graphed
#' @param zoom_ensg: if `chr` is set this will graph around the gene described by endembl ID
#' @param zoom_gene: if `chr` is set this willgraph around the gene described by name
#' @param zoom_left: if `chr` is set this will graph points to the right of this base pair number
#' @param zoom_right: if `chr` is set this will graph points to the left of this base pair number
#' @param gene_tag_1 numeric, annotate the SNPs with P-values more extreme than log(gene_tag_p)
#' @param gene_tag_2 numeric, annotate the SNPs with P-values more extreme than log(gene_tag_p)
#' @param sig_line1 numeric, draw a horizontal line at -log(sig_line1)
#' @param sig_line2 numeric, draw a horizontal line at -log(sig_line2)
#' @param sig_line1_color color for sig_line1
#' @param sig_line2_color color for sig_line2
#' @param draw_genes draw lines along the length of each gene instead of a dot at the midpoint
#' @param data1_direction direction for data1 to be drawn, 1 is up and -1 is down
#' @importFrom ggplot2 ggplot aes theme_bw guides geom_label_repel geom_hline guide_legend scale_x_continuous scale_y_continuous scale_colour_manual theme element_text element_line element_blank xlab ylab expand_limits geom_hline geom_point
#' @export
nashville.plot <- function(data1, data2 = NULL, data1_direction = 1, map_df = "37",
                           chr=NULL, zoom_ensg=NULL, zoom_gene=NULL, zoom_left = 0, zoom_right = Inf,
                           draw_genes = FALSE, labels_cat=NULL, group_color=NULL, sig_line1 = NULL,
                           sig_line2 = NULL, sig_line1_color = 'black', sig_line2_color = 'black',
                           gene_tag_1 = -Inf, gene_tag_2 = -Inf, config = data.frame(V1=NA,V2=NA,V3=NA), ...) {
  if (data1_direction == 1) {
    data2_direction = -1
  } else {
    data1_direction = -1
    data2_direction = 1
  }

  if (map_df == "37") {
    map_df <- gene.build.37
  } else if (map_df == "38") {
    map_df <- gene.build.38
  } else {
    stop("map_df must be '37' or  '38'")
  }

  # step 1) Combine the data
  data1$source <- 1
  if(is.null(data2)) {
    full.obj <- data1
  } else {
    data2$source <- 2
    full.obj <- rbind(data1, data2)
  }

  # step 2) Determine absolute positions
  max_pos <- tapply(full.obj$BP, full.obj$CHR, max)
  chr_shift <- head(c(0,cumsum(as.numeric(max_pos))),-1)
  if (all(chr_shift==0)) {
    chr_shift <- rep(0, times=max(full.obj$CHR))
  }
  ifelse(is.null(chr),
         full.obj$absolute <- full.obj$BP + chr_shift[full.obj$CHR],
         full.obj$absolute <- full.obj$BP)

  #step 3) Get tick mark positions
  full.obj$absolute <- full.obj$absolute * 1e-6                                                 #megabase
  full.obj$gene.start <- full.obj$gene.start * 1e-6                                                 #megabase
  full.obj$gene.end <- full.obj$gene.end * 1e-6                                                 #megabase
  full.obj$BP <- full.obj$BP * 1e-6
  axis_set <- aggregate(full.obj$absolute, by = list(full.obj$CHR), mean)
  names(axis_set) <- c("CHR", "center")

  # step) decide what will get tagged
  for_tag1 <- full.obj[which(full.obj$P < gene_tag_1 & full.obj$source == 1), ]
  for_tag1 <- for_tag1[order(for_tag1$gene.name, -log(for_tag1$P, 10)), ]
  for_tag1 <- for_tag1[!duplicated(for_tag1$gene.name),]
  for_tag2 <- full.obj[which(full.obj$P < gene_tag_2 & full.obj$source == 2), ]
  for_tag2 <- for_tag2[order(for_tag2$gene.name, -log(for_tag2$P, 10)), ]
  for_tag2 <- for_tag2[!duplicated(for_tag2$gene.name),]
  for_tag <- rbind(for_tag1, for_tag2)

  #step 4) Decide what to plot
  x_axis_name <- "Chromosome"
  if(!is.null(chr)) {
    full.obj <- full.obj[which(full.obj$CHR == chr), ]
    x_axis_name <- paste0("Chromosome ",  chr, " (MB)")
    if(!is.null(zoom_gene)) {
      bounds <- get_gene_bounds(zoom_gene, map_df)
      bounds[1] <- bounds[1] - 3000
      bounds[2] <- bounds[2] + 3000
    } else if(!is.null(zoom_ensg)) {
      bounds <- get_gene_bounds_ensg(zoom_ensg, map_df)
      bounds[1] <- bounds[1] - 3000
      bounds[2] <- bounds[2] + 3000
    } else if(zoom_right != 0) {
      bounds <- c(zoom_left, zoom_right)
    } else {
      bounds <- c(zoom_left, zoom_right)
    }
    bounds <- bounds * 1e-6       #megabase
    full.obj <- full.obj[which(full.obj$BP >= bounds[1] & full.obj$BP <= bounds[2]), ]
  }
  #TODO build labels cat from Label if available
  #step 5) Set colors and legend labels
  #if(is.null(labels_cat)){
    labels_cat <- c(unique(as.character(full.obj$group)))
    labels_cat[grep("^Single SNP [0-1]$", labels_cat)] <- "Single SNP"
  #}
  #if(is.null(group_color)) {
    ##########################data.frame(filename, tissue, color)
    group_color <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(labels_cat))
    #names(group_color) <- labels_cat
    #print(group_color)

    color_df <- data.frame(color = group_color, group = labels_cat)

    #group_color[grep("^Single SNP$", labels_cat)] <- 'black'
    #group_color[grep("Single \u200bSNP", labels_cat)] <- 'darkgrey'
  #}

    #### we have defined group color and labels_cat
    #### we must replace these values with values provided by config if they exist
    colnames(config) <- c('source', 'group', 'color')
    full.obj <- merge(full.obj, color_df, by.x = 'group', by.y = 'group', all.x = T)
    full.obj <- merge(full.obj, config, by.x = "group", by.y = "source", all.x = TRUE)
    full.obj$group <- ifelse(is.na(full.obj$group.y), full.obj$group, full.obj$group.y)
    full.obj$color <- ifelse(is.na(full.obj$color.y), full.obj$color.x, full.obj$color.y)

    for(i in 1:length(full.obj$color)) {
      s <- unlist(strsplit(full.obj[i, 'color'], "|", fixed = TRUE))
      if(length(s) == 2) {
        if(full.obj[i, 'CHR'] %% 2 == 0) {
          full.obj[i, 'color'] <- s[1]
        } else {
          full.obj[i, 'color'] <- s[2]
        }
      }
    }
    #full.obj$color2 <- ifelse(length(unlist(strsplit(full.obj$color, "|", fixed=T)) == 2),
    #                          print(full.obj$color),
    #                          full.obj$color)


    #  for(label in labels_cat) {
    #    print(paste0("label: ", label))
    #    line <- config[which(config$source == label), ]
    #    print(paste0("line: ", line))
    #    labels_cat[which(labels_cat == label)] <- line$group
    #    group_color[label] <- line$color
    #  }
    # names(group_color) <- labels_cat
    # print(group_color)

    #df <- data.frame(group = labels_cat, color = group_color)
    #colnames(config) <- c('source', 'group', 'color')
    #config <- merge(config, df, by.x = 'filename', by.y = )
    #group_color
    #labels_cat

  #step 6) Plot
  plt <- ggplot() + theme_bw()
  #plt <- plt + scale_y_continuous() #breaks
  plt <- plt + xlab(x_axis_name)
  plt <- plt + ylab(expression(log["10"]*italic((p))~"&"~-log["10"]*italic((p))))
  #plt <- plt + expand_limits()
  plt <- plt + ggrepel::geom_label_repel(data = for_tag,
                                         aes(x = absolute, y = log(P, 10), label = gene.name))
  plt <- plt + theme(axis.text.x=element_text(size=12, color='black'),
                     axis.text.y=element_text(size=12, color='black'),
                     axis.title.x=element_text(size = 12, face="bold", color="black"),
                     axis.title.y=element_text(size = 12, face="bold", color="black"),
                     axis.ticks.x=element_line())
  plt <- plt + theme(panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank())
  #plt <- plt + scale_colour_manual(name = "Tissue",
  #                                 values = c(group_color))
  #                                 labels = c(labels_cat))
  plt <- plt + scale_color_identity("Tissue", guide = "legend",
                                    labels = full.obj$group,
                                    breaks = full.obj$color)
  plt <- plt + guides(shape = "none",
                      size = "none",
                      colour = guide_legend(override.aes = list(size=6))) #reverse = TRUE
  plt <- plt + plot.mh(data = subset(full.obj, source == 1),
                       direction = data1_direction,
                       draw_genes)
  if(!is.null(data2)) {
    plt <- plt + plot.mh(data = subset(full.obj, source == 2),
                         direction = data2_direction,
                         draw_genes)
  }
  if(!is.null(sig_line1)) {
    plt <- plt + geom_hline(aes(yintercept = data1_direction * sig_line1), color = sig_line1_color)
  }
  if(!is.null(sig_line2)) {
    plt <- plt + geom_hline(aes(yintercept = data2_direction * sig_line2), color = sig_line2_color)
  }
  if(is.null(chr)) {
    plt <- plt + scale_x_continuous(label = axis_set$CHR, breaks = axis_set$center)
  } else {
    plt <- plt + scale_x_continuous(label = waiver(), breaks = waiver())
  }
  plt <- plt + geom_hline(aes(yintercept = 0), size = 1)
  plt
}
