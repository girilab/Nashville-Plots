suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggbreak))

#' Read config file
#'
#' Reads in configuration info from tab separated values file
#' @param file a file path to a 3 column tab separated values file without headers where the first column contains file names, the second column contains group names, and the 3rd column contains colors for those groups to be plotted in.
#' @export
read_config <- function(file) {
  lines <- readLines(file)
  lines <- lines[nzchar(trimws(lines))]  # drop blank lines

  parsed <- lapply(lines, function(ln) {
    parts <- strsplit(trimws(ln), "\\s+")[[1]]
    v1 <- if (length(parts) >= 1) parts[1] else NA_character_
    v2 <- if (length(parts) >= 2) parts[2] else NA_character_
    v3 <- if (length(parts) >= 3) parts[3] else NA_character_
    data.frame(V1 = v1, V2 = v2, V3 = v3, stringsAsFactors = FALSE)
  })

  config <- do.call(rbind, parsed)

  if (length(unique(config$V2)) != length(config$V2)) {
    stop("Each element of the group column must be unique.")
  }
  config
}

# converts CHR as string to number and handles X, Y, and some other less common chromosome designations
chr_as_numeric <- function(chr) {
  chr <- as.character(chr)
  chr <- sub("^chr", "", chr, ignore.case = TRUE)
  chr[chr %in% c("X", "x")] <- "23"
  chr[chr %in% c("Y", "y")] <- "24"
  chr[chr %in% c("M", "MT", "m", "mt")] <- "25"
  as.numeric(chr)
}

#' Get gene bounds
#'
#' Determines boundaries for plotting a gene by gene name
#'
#' @param gene_name name of a gene
#' @param map_df data frame mapping gene names to ensembl ID.
#' @param chr chromosome to find the gene in, optional
get_gene_bounds <- function(gene_name, map_df, chr = NULL) {
  hits <- map_df[!is.na(map_df$Gene) & map_df$Gene == gene_name, , drop = FALSE]
  if (!is.null(chr)) {
    hits <- hits[!is.na(chr_as_numeric(hits$CHR)) & chr_as_numeric(hits$CHR) == chr_as_numeric(chr), , drop = FALSE]
    if (nrow(hits) == 0) {
      stop(paste0(gene_name, " not found in Chromosome ", chr))
    }
  }
  if (nrow(hits) == 0) {
    stop(paste0("Gene not found: ", gene_name))
  }
  lower <- unique(hits$START_POS)
  upper <- unique(hits$END_POS)
  if (length(lower) > 1 || length(upper) > 1) {
    stop(paste0("Multiple genes with name: ", gene_name, " Try giving an ENSG."))
  }
  pad <- (upper - lower) * 0.05
  c(lower - pad, upper + pad)
}

#' Get gene bounds
#'
#' Determines boundaries for plotting a gene  by ensembl ID
#'
#' @param ensg ensembl ID
#' @param map_df data frame mapping gene names to Ensembl ID.
get_gene_bounds_ensg <- function(ensg, map_df, chr = NULL) {
  hits <- map_df[!is.na(map_df$ENSG) & map_df$ENSG == ensg, , drop = FALSE]
  if (!is.null(chr)) {
    hits <- hits[!is.na(chr_as_numeric(hits$CHR)) & chr_as_numeric(hits$CHR) == chr_as_numeric(chr), , drop = FALSE]
    if (nrow(hits) == 0) {
      stop(paste0(ensg, " not found in Chromosome ", chr))
    }
  }
  if (nrow(hits) == 0) {
    stop(paste0("ENSG not found: ", ensg))
  }
  lower <- unique(hits$START_POS)
  upper <- unique(hits$END_POS)
  if (length(lower) > 1 || length(upper) > 1) {
    stop(paste0("Multiple genes with name: ", ensg, " Try giving a gene name"))
  }
  pad <- (upper - lower) * 0.05
  c(lower - pad, upper + pad)
}

#' Read GWAS File
#'
#' @param file path to gwas results as obtained from plink --assoc
#' @param samp a number between 0 and 1 describing what percent of entries with p-value > 0.1 to randomly keep for plotting
#' @param chr column of chromosome of SNP location, numeric
#' @param p column of p-values of SNP, numeric
#' @param bp column of base pair of SNP location, numeric
#' @param snp.name column of SNP names
#' @param config data frame of plot configuration
#' @param color color for points from this file
#' @param shape shape for points from this file
#' @importFrom data.table fread
#' @export
read_gwas_file <- function(file, config=data.frame(V1=NA,V2=NA,V3=NA), samp = NULL, chr = "CHR", p = "P", bp = "BP", snp.name = NA, color=NA, shape=NA) {
  gwas <- data.table::fread(file)
  gwas[[p]] <- as.numeric(gwas[[p]])
  if(!is.null(samp)) {
    gwas <- gwas[which(gwas[[p]] < 0.5), ]
    greater <- gwas[sample(which(gwas[[p]] > 0.1), round(samp*length(which(gwas[[p]] > 0.1)))), ]
    lesser <- gwas[which(gwas[[p]] <= 0.1), ]
    lesser <- lesser[which(lesser[[p]] < 0.5), ]
    gwas <- rbind(lesser, greater)
  }
  snp.values <- rep(NA_character_, nrow(gwas))
  if (!is.na(snp.name) && snp.name %in% names(gwas)) {
    snp.values <- gwas[[snp.name]]
  }
  gwas.obj <- make.valid.object(CHR = as.numeric(gwas[[chr]]),
                                P = as.numeric(gwas[[p]]),
                                BP = as.numeric(gwas[[bp]]),
                                group=file,
                                color=color,
                                shape=shape,
                                datatype="gwas",
                                config=config,
                                snp.name = snp.values)
  gwas.obj
}
#' Read MetaXcan Folder
#'
#' This function reads in MetaXcan or PrediXcan results files by tissue
#' @param directory a string representation of a file path containing PrediXcan output
#' @param map_df a file path containing a mapping gene names to Ensembl ID. The header should be "Gene ENSG CHR START_POS END_POS"
#' @param label a list of labels to use instead of file names for graphing. These should probably be in alphabetical order.
#' @param pattern a regex describing files to be read from "directory"
#' @param samp a number between 0 and 1 describing what percent of entries with p-value > 0.1 to randomly keep for plotting
#' @export
read_metaXcan_folder <-  function(directory, config=data.frame(V1=NA,V2=NA,V3=NA), map_df = '37', pattern='*.csv$', samp = NULL, color=NA, shape=NA) {
  if (map_df == "37") {
    appended_gene_d <- gene.build.37
  } else if (map_df == "38") {
    appended_gene_d <- gene.build.38
  } else {
    stop("map_df must be '37' or '38'")
  }
  files <- list.files(directory, pattern=pattern)
  if (length(files) == 0) {
    stop(paste0("No files matching pattern '", pattern, "' found in ", directory))
  }
  tissues <- data.frame()
  for (i in seq_along(files)){
    tissue_file <- read.csv(file=file.path(directory, files[i]), header = TRUE)
    tissue_file$dtype <- files[i]
    tissues <- rbind(tissues, tissue_file)
  }
  tissues$ENSG <- sub("\\..*", "", tissues[["gene"]])
  appended_gene_d <- merge(tissues, appended_gene_d, by.x = "ENSG", by.y = "ENSG")

  if(!is.null(samp)) {
    greater <- appended_gene_d[sample(which(appended_gene_d[["pvalue"]]>0.1), round(samp*length(which(appended_gene_d[["pvalue"]]>0.1)))), ]
    lesser <- appended_gene_d[which(appended_gene_d[["pvalue"]]<=0.1), ]
    appended_gene_d <- rbind(lesser, greater)
  }

  appended_gene_d$CHR <- as.numeric(appended_gene_d[["CHR"]])
  meta.obj <- make.valid.object(CHR = appended_gene_d[["CHR"]],
                                P = appended_gene_d[["pvalue"]],
                                gene.start = appended_gene_d[["START_POS"]],
                                gene.end = appended_gene_d[["END_POS"]],
                                group = appended_gene_d[["dtype"]],
                                gene.name = appended_gene_d[["gene_name"]],
                                color=color,
                                shape=shape,
                                datatype="meta",
                                config=config)
  meta.obj
}

#############
#############
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
make.valid.object <- function(CHR, P, group, config, BP=NA, gene.start=NA, gene.end=NA, gene.name=NA, color=NA, shape=NA, snp.name = NA, datatype=NA) {
  ifelse(all(is.numeric(CHR)), NA, stop("CHR must be numeric"))
  ifelse(all(is.numeric(P)), NA, stop("P must be numeric"))
  ifelse(all(is.numeric(BP)) | all(is.na(BP)), NA, stop("BP must be numeric"))
  ifelse(all(is.numeric(gene.start)) | all(is.na(gene.start)), NA, stop("gene.start must be numeric"))
  ifelse(all(is.numeric(gene.end)) | all(is.na(gene.end)), NA, stop("gene.end must be numeric"))
  if (length(shape) != 1 && length(shape) != length(P)) {
    stop("shape must be the same length as P or a single value")
  }
  #make sure that either BP or gene.start+gene.end exists
  if(all(is.na(gene.start)) && all(is.na(gene.end))) {
    stopifnot(is.numeric(BP))
  }
  if(all(is.na(BP))) {
    stopifnot(is.numeric(gene.start) && is.numeric(gene.end))
    BP <- (gene.start + gene.end)/2  # TODO this should be delivered from map_df
  }

  # Step 1: auto-assign a color to every row based on group/CHR
  if (identical(datatype, "gwas")) {
    group_color <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(unique(CHR)))
    color <- group_color[as.numeric(as.factor(CHR %% 2))]
  } else if (identical(datatype, "meta")) {
    #group_color <- colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(length(unique(group)))
    # unused for now
    other_palette <- c('#e6194B', '#3cb44b', '#ffe119', '#4363d8',
                       '#f58231', '#911eb4', '#42d4f4', '#f032e6',
                       '#bfef45', '#fabed4', '#469990', '#dcbeff',
                       '#9A6324', '#fffac8', '#800000', '#aaffc3',
                       '#808000', '#ffd8b1', '#000075', '#a9a9a9')
    # meta palette is paletteer_dynamic("cartography::multi.pal", 20)
    meta_palette <- c("#CB7C77FF", "#68D359FF", "#6B42C8FF", "#C9D73DFF",
                      "#C555CBFF", "#AED688FF", "#502E71FF", "#C49A3FFF",
                      "#6A7DC9FF", "#D7652DFF", "#7CD5C8FF", "#C5383CFF",
                      "#507D41FF", "#CF4C8BFF", "#5D8D9CFF", "#722E41FF",
                      "#C8B693FF", "#33333CFF", "#C6A5CCFF", "#674C2AFF")
    group_color <- colorRampPalette(meta_palette)(length(unique(group)))
    color <- group_color[as.numeric(as.factor(group))]
  }
  # Step 2: override with user-supplied colors from config where V3 is non-NA/non-empty.
  # Applied after auto-colors so only explicitly specified rows are overridden.

  #set shape if not provided
  if(any(is.na(shape))) {
    shape <- as.factor(CHR %% 2)
  }
  gene.name <- ifelse(is.na(gene.name), gene.name, as.character(gene.name))
  temp <- data.frame(gene.name, group, CHR, BP, P, gene.start, gene.end, color, shape, datatype, snp.name)
  object <- merge(temp, config, by.x='group', by.y = 'V1', all.x=T)
  object$names <- object$V2
  object[which(is.na(object$names)), "names"] <- object[which(is.na(object$names)), "group"]

  # Override auto-assigned colors with any user-supplied colors from config (V3)
  has_user_color <- !is.na(object$V3) & nzchar(trimws(object$V3))
  if (any(has_user_color)) {
    # Work group-by-group so each group's palette is applied independently
    for (grp in unique(object$group[has_user_color])) {
      grp_rows <- object$group == grp & has_user_color
      palette  <- trimws(strsplit(object$V3[grp_rows][1], ",")[[1]])
      chr_vals <- object$CHR[grp_rows]
      # Map CHR to palette index: CHR 1->1, 2->2, ..., wrapping by palette length
      object$color[grp_rows] <- palette[(chr_vals - 1) %% length(palette) + 1]
    }
  }
  object$V2 <- NULL
  object$V3 <- NULL
  object[which(object$P != 0), ]
}

#' plot manhattan
#'
#' draws a manhattan plot onto the nashville plot
#' @param data data frame from make.valid.object()
#' @param direction direction for data to be drawn
#' @param draw_genes bool plot with geam_line
#' @importFrom ggplot2 ggplot aes geom_point
plot.mh <- function(data, direction, draw_genes) {

  if(draw_genes && any(!is.na(data$gene.name)) && any(!is.na(data$group))) {
    data_copy <- data
    data_copy$gene.start <- data_copy$gene.end
    data <- rbind(data, data_copy)
    mh <- geom_line(data = data,
                    aes(x=gene.start,
                        y=direction * -log(P, 10),
                        color=color,
                        group=interaction(gene.name, group)))
  } else {
    # data[sample(1:nrow(data)), ]
    mh <- geom_point(data=data[sample(1:nrow(data)), ],#  |> sample_n(size = nrow(df), replace = FALSE), # the sampling is done to randomize the order of rows to mix colors better
                     aes(x=absolute,
                         y=direction * -log(P, 10),
                         color=color,
                         shape=as.factor(CHR %% 2)))
  }
}

#' Determine absolute position
#'
#' converts from chromosome and base pair position to absolute base pair position
#' @param full.obj data.frame from make.valid.object()
determine_abs_position <- function(full.obj) {
  max_pos <- tapply(full.obj[["BP"]], full.obj[["CHR"]], max)
  chr_shift <- head(c(0,cumsum(as.numeric(max_pos))),-1)
  names(chr_shift) <- names(max_pos)
  if (all(chr_shift==0)) {
    chr_shift <- setNames(rep(0, length(max_pos)), names(max_pos))
  }
  chr_shift
}

#' converts to megabase values from base values
#'
#' converts from chromosome and base pair position to absolute base pair position
#' @param full.obj data.frame from make.valid.object()
to_megabase <- function(df) {
  df$absolute <- df$absolute * 1e-6
  df$gene.start <- df$gene.start * 1e-6
  df$gene.end <- df$gene.end * 1e-6
  df$BP <- df$BP * 1e-6
  df
}

# find_y_break <- function(log10p_vals,
#                          direction        = 1,            # 1 for "up" data, -1 for "down" data
#                          sig_threshold    = -log10(5e-8),  # ~7.3
#                          min_gap          = 1,             # minimum gap size to bother breaking
#                          max_outlier_frac = 0.10) {         # never isolate more than this fraction of hits
#
#   vals <- direction * log10p_vals
#   sig_vals <- sort(vals[vals > sig_threshold])
#   n <- length(sig_vals)
#
#   if (n < 2) {
#     message("Too few significant hits to determine a break.")
#     return(NULL)
#   }
#
#   gaps <- diff(sig_vals)
#   max_outliers  <- max(1, floor(n * max_outlier_frac))
#   candidate_idx <- seq(from = max(1, n - max_outliers), to = n - 1)
#
#   local_gaps <- gaps[candidate_idx]
#   best       <- which.max(local_gaps)
#   gap_idx    <- candidate_idx[best]
#   max_gap    <- gaps[gap_idx]
#
#   break_lo <- sig_vals[gap_idx]        # bottom of the break, in magnitude space
#   break_hi <- sig_vals[gap_idx + 1]    # top of the break, in magnitude space
#
#   if (max_gap < min_gap) {
#     message(sprintf(
#       "Largest gap near the tail is %.1f (below min_gap = %g); no break applied.",
#       max_gap, min_gap
#     ))
#     return(NULL)
#   }
#   message(sprintf(
#     "Break placed between %.1f and %.1f (gap = %.1f -log10 units), isolating %d/%d hits",
#     break_lo, break_hi, max_gap, n - gap_idx, n
#   ))
#
#   # Convert back out of magnitude space and re-sort, since flipping the sign
#   # for "down" data reverses the numeric order of break_lo/break_hi.
#   as.list(sort(direction * c(break_lo, break_hi)))
# }

#' Parse a tag specification into plain names and gene/tissue tuples
#'
#' `tag_genes`/`tag_names1`/`tag_names2` can be supplied two ways:
#'  \itemize{
#'    \item a plain character vector, e.g. \code{c("BRCA1", "TP53")}. Each
#'      element is treated as an unqualified gene (or SNP/group) name
#'      if the same gene appears in multiple tissues, all of its
#'      rows are eligible to be tagged, but they're still subject to the usual
#'      per-gene/per-peak collapsing (see \code{\link{build_tag_subset}}).
#'    \item a \code{list} whose elements are either a single name (character
#'      scalar, same as above) or a length-2 character vector
#'      \code{c(gene_name, tissue_name)} that uniquely identifies one
#'      tissue's point for that gene, e.g.
#'      \code{list("TP53", c("BRCA1", "Liver"), c("BRCA1", "Muscle"))}.
#'      Tuple-specified entries are exempted from the per-gene collapsing, so
#'      this is how to force a gene to be tagged in more than one tissue.
#'  }
#' @param tag_genes see above
parse_tag_spec <- function(tag_genes) {
  plain <- character(0)
  tuple_gene <- character(0)
  tuple_tissue <- character(0)

  if (is.null(tag_genes) || length(tag_genes) == 0) {
    return(list(plain = plain, tuple_gene = tuple_gene, tuple_tissue = tuple_tissue))
  }

  if (is.list(tag_genes)) {
    for (item in tag_genes) {
      item <- as.character(item)
      if (length(item) == 1) {
        plain <- c(plain, item)
      } else if (length(item) == 2) {
        tuple_gene <- c(tuple_gene, item[1])
        tuple_tissue <- c(tuple_tissue, item[2])
      } else {
        stop("Each tag_genes/tag_names element must be a single name, or a gene/tissue pair of length 2.")
      }
    }
  } else {
    # plain atomic vector: every element is a standalone gene/name
    plain <- as.character(tag_genes)
  }

  list(plain = plain, tuple_gene = tuple_gene, tuple_tissue = tuple_tissue)
}

#' Build a tag subset for one data source
#'
#' @param obj subset of full.obj belonging to a single datasource
#' @param tag_genes gene/tissue tag specification; see \code{\link{parse_tag_spec}}
#'   for the accepted formats. Matched against \code{obj$gene.name} (plain
#'   entries) or, for a tuple, against \code{obj$gene.name} together with
#'   \code{obj$group} (the tissue/dataset a row came from); plain entries also
#'   match against \code{obj$names} (SNP/group names for GWAS data).
#' @param gene_tag numeric P-value threshold; rows with P < gene_tag are tagged
#' @param peak_window x range to be considered a "peak". Only 1 point is labeled
#'   per peak window, regardless of tissue. 
#'   To force a gene to be labeled in more than one tissue 
#'   (even at/near the same peak), pass explicit \code{c(gene_name, tissue_name)} 
#'   tuples in \code{tag_genes}; tuple-specified points always bypass peak 
#'   suppression, each uniquely identified by its own gene+tissue pair.
build_tag_subset <- function(obj, tag_genes, gene_tag, peak_window = 500000) {

  if (nrow(obj) == 0) return(obj[0, , drop = FALSE])

  has_gene <- !is.na(obj$gene.name)
  has_name <- !is.na(obj$names)

  spec <- parse_tag_spec(tag_genes)
  has_tags <- length(spec$plain) > 0 || length(spec$tuple_gene) > 0
  tuple_key_spec <- paste(spec$tuple_gene, spec$tuple_tissue, sep = "\r")

  # rows matching a plain (tissue-agnostic) name
  list_hit_plain <- rep(FALSE, nrow(obj))
  if (length(spec$plain) > 0) {
    list_hit_plain <- (has_gene & obj$gene.name %in% spec$plain) |
      (has_name & obj$names %in% spec$plain) |
      (!is.na(obj$snp.name) & obj$snp.name %in% spec$plain)
  }

  # rows matching an explicit (gene, tissue) tuple
  list_hit_tuple <- rep(FALSE, nrow(obj))
  if (length(spec$tuple_gene) > 0) {
    obj_key <- paste(obj$gene.name, obj$group, sep = "\r")
    list_hit_tuple <- has_gene & obj_key %in% tuple_key_spec
  }

  list_hit <- list_hit_plain | list_hit_tuple

  # rows exceeding the significance threshold
  thresh_hit <- rep(FALSE, nrow(obj))
  if (is.numeric(gene_tag) && length(gene_tag) == 1 && !is.na(gene_tag)) {
    thresh_hit <- !is.na(obj$P) & obj$P < gene_tag & (has_gene | has_name)
  }

  keep <- list_hit | thresh_hit
  for_tag <- obj[keep, , drop = FALSE]
  tuple_flag <- list_hit_tuple[keep]

  if (nrow(for_tag) == 0) return(for_tag)

  # De-duplication: rows uniquely identified via a (gene, tissue) tuple keep
  # one row PER tissue (their dedup key includes `group`). 
  tuple_rows     <- for_tag[tuple_flag, , drop = FALSE]
  non_tuple_rows <- for_tag[!tuple_flag, , drop = FALSE]

  if (nrow(non_tuple_rows) > 0) {
    non_tuple_rows <- non_tuple_rows[order(non_tuple_rows$gene.name, -log(non_tuple_rows$P, 10), decreasing = TRUE), ]
    non_tuple_rows <- non_tuple_rows[!duplicated(non_tuple_rows[, c("gene.name", "snp.name")]), ]
  }
  if (nrow(tuple_rows) > 0) {
    tuple_rows <- tuple_rows[order(tuple_rows$gene.name, tuple_rows$group, -log(tuple_rows$P, 10), decreasing = TRUE), ]
    tuple_rows <- tuple_rows[!duplicated(tuple_rows[, c("gene.name", "group", "snp.name")]), ]
  }
  for_tag <- rbind(non_tuple_rows, tuple_rows)

  # Sort by significance so the top hit per region comes first
  for_tag <- for_tag[order(for_tag$P), ]

  kept <- rep(FALSE, nrow(for_tag))
  peak_regions <- data.frame(chr = character(), start = numeric(), end = numeric())

  for (i in seq_len(nrow(for_tag))) {
    row <- for_tag[i, ]

    # Named hits are always kept regardless of proximity to a peak
    is_named <- has_tags &&
      (((!is.na(row$gene.name) && row$gene.name %in% spec$plain) |
          (!is.na(row$names)     && row$names     %in% spec$plain)) ||
         (!is.na(row$gene.name) &&
            paste(row$gene.name, row$group, sep = "\r") %in% tuple_key_spec) || 
            (!is.na(row$snp.name) && row$snp.name %in% spec$plain))

    if (is_named) {
      kept[i] <- TRUE
      next
    }

    # Check if this variant falls within an already-claimed peak window
    if (nrow(peak_regions) > 0) {
      in_peak <- peak_regions$chr == row$CHR &
        peak_regions$start <= row$BP &
        peak_regions$end   >= row$BP
      if (any(in_peak)) next  # suppressed — inside an existing peak
    }

    # This is the top hit for a new peak; claim a window around it
    kept[i] <- TRUE
    peak_regions <- rbind(peak_regions, data.frame(
      chr   = row$CHR,
      start = row$BP - peak_window,
      end   = row$BP + peak_window
    ))
  }

  for_tag <- for_tag[kept, , drop = FALSE]
  if (nrow(for_tag) > 100) {
    message(paste0("Warning: ", nrow(for_tag), "data points tagged.") )
  }
  for_tag
}

#' Generate a Nashville plot
#'
#' This function returns a ggplot2 graph object of the Nashville plot
#'
#' @param data1 data to be plotted generated by read_metaXcan_folder or read_gwas_file.
#' @param data2 data to be plotted generated by read_metaXcan_folder or read_gwas_file. Optional
#' @param map_df: 37 or 38 depending on the Human Genome assembly reference number, degfault = 37
#' @param chr: chromosome number. If provided, only data from that chromosome will be graphed.
#' @param zoom_ensg: if `chr` is set this will graph around the gene described by endembl ID
#' @param zoom_gene: if `chr` is set this will graph around the gene described by name
#' @param zoom_left: if `chr` is set this will graph points to the right of this base pair number
#' @param zoom_right: if `chr` is set this will graph points to the left of this base pair number
#' @param tag_names1 gene/tissue tag specification for data1: either a plain character
#'     vector of gene names (or group/tissue names for GWAS SNPs) to label on the plot
#'     regardless of significance, or a \code{list} mixing plain names with explicit
#'     \code{c(gene_name, tissue_name)} tuples, e.g.
#'     \code{list("TP53", c("BRCA1", "Liver"), c("BRCA1", "Muscle"))}. A plain name
#'     tags that gene's single most significant point. a tuple
#'     uniquely tags that gene's point in that specific tissue, which is how to label
#'     the same gene in more than one tissue.
#' @param tag_names2 same as \code{tag_names1}, applied to data2.
#' @param tag_threshold1 numeric, in data1 annotate items with P-values more extreme than tag_threshold1
#'    (default -Inf, i.e. no automatic threshold tagging). Use \code{tag_names1}
#'    tuples to label additional tissues for the same gene.
#' @param tag_threshold2 numeric, in data2 annotate items with P-values more extreme than tag_threshold2
#'     (default -Inf, i.e. no automatic threshold tagging). Use \code{tag_names2}
#'     tuples to label additional tissues for the same gene.
#' @param sig_line1 numeric, draw a horizontal line at -log(sig_line1)
#' @param sig_line2 numeric, draw a horizontal line at -log(sig_line2)
#' @param sig_line1_color color for sig_line1
#' @param sig_line2_color color for sig_line2
#' @param draw_genes draw lines along the length of each gene instead of a dot at the midpoint
#' @param config config object as produced by read_config()
#' @param y_min minimum y-value to plot
#' @param y_max maximum y-value to plot
#' @param y_ticks number of breaks on the graph y-axis
#' @param data1_direction direction for data1 to be drawn, "up" or "down"
#' @param axis_break_scale multiplier for middle chunk of the plot if there is an axis_break
#' @param axis_breaks where to cut the y-axis for shrinking extremes of the plot,
#'  accepts FALSE, 'auto', or a list of y values of places to adjust the relative scale
#' @importFrom ggplot2 ggplot aes theme_bw guides geom_hline guide_legend scale_x_continuous scale_y_continuous scale_colour_manual theme element_text element_line element_blank xlab ylab expand_limits geom_hline geom_point
#' @importFrom ggrepel geom_label_repel
#' @export

nashville.plot <- function(data1, data2=NULL, map_df="37", chr=NULL, zoom_ensg=NULL, zoom_gene=NULL,
                           zoom_left=0, zoom_right=Inf,
                           tag_names1=NULL, tag_names2=NULL,
                           tag_threshold1=-Inf, tag_threshold2=-Inf,
                           sig_line1=NULL, sig_line2=NULL,
                           sig_line1_color='black', sig_line2_color='black', draw_genes=FALSE,
                           config=data.frame(V1=NA,V2=NA,V3=NA), y_min=NA, y_max=NA, y_ticks=NA,
                           data1_direction="up", axis_breaks = FALSE, axis_break_scale = 5, ...) {
  validate_axis_breaks <- function(axis_breaks, full.obj) {
    if (identical(axis_breaks, FALSE) || identical(axis_breaks, "auto") || length(axis_breaks) == 0) return(invisible(NULL))
    if (!is.numeric(axis_breaks)) {
      stop("axis_breaks must be FALSE, 'auto', or a numeric vector of length 1 or 2", call. = FALSE)
    }

    logP_range <- range(full.obj$logP, na.rm = TRUE)
    data_min   <- logP_range[1]
    data_max   <- logP_range[2]

    for (b in axis_breaks) {
      if (b <= data_min || b >= data_max) {
        cat(sprintf(
          "axis_breaks value %g is outside the data range [%g, %g].\n",
          b, data_min, data_max
        ), sprintf(
          "  Break points must be within the plotted logP range.\n"
        ), sprintf(
          "  Your data spans %.1f to %.1f -- try axis_breaks = c(%g, %g) or set axis_breaks = FALSE.",
          data_min, data_max,
          round(data_min * 0.75),   # suggest 75% of range as sensible defaults
          round(data_max * 0.75)
        ), call. = FALSE)
      }
    }
    invisible(NULL)
  }
  map_df <- switch(map_df,
                   "37" = {map_df = gene.build.37},
                   "38" = {map_df = gene.build.38},
                   stop("map_df must be '37' or '38'"))
  infer_chr_from_map <- function(id, id_col) {
    hits <- map_df[!is.na(map_df[[id_col]]) & map_df[[id_col]] == id, , drop = FALSE]
    if (nrow(hits) == 0) {
      (paste0(id, " not found in map_df; cannot infer chromosome. ",
              "Pass `chr` explicitly."))
    }
    chrs <- unique(chr_as_numeric(hits$CHR))
    if (length(chrs) > 1) {
      stop(paste0(id, " maps to multiple chromosomes (", paste(chrs, collapse=", "),
                  "); pass `chr` explicitly to disambiguate."))
    }
    chrs[1]
  }
  if (is.null(chr) && !is.null(zoom_gene)) {
    chr <- infer_chr_from_map(zoom_gene, "Gene")
    message(paste0("chr not specified; inferred chr = ", chr, " from zoom_gene = '", zoom_gene, "'"))
  }else if (is.null(chr) && !is.null(zoom_ensg)) {
    chr <- infer_chr_from_map(zoom_ensg, "ENSG")
    message(paste0("chr not specified; inferred chr = ", chr, " from zoom_ensg = '", zoom_ensg, "'"))
  }
  data1$datasource <- 1
  if (!is.null(data2)) { data2$datasource <- 2 }
  full.obj <- rbind(data1, data2)
  full.obj <- merge(full.obj, config, by.x="group", by.y="V2", all.x = TRUE)
  chr_shift <- determine_abs_position(full.obj)
  ifelse(is.null(chr),
         full.obj$absolute <- full.obj$BP + chr_shift[full.obj$CHR],
         full.obj$absolute <- full.obj$BP)

  if(data1_direction == 'up') {
    data1_direction = 1
  } else
    if(data1_direction == 'down') {
      data1_direction = -1
    } else {
      stop("direction should be either 'up' or 'down'")
    }
  data2_direction = data1_direction * -1

  full.obj[which(full.obj$datasource == 1), "logP"] <- -log10(full.obj[which(full.obj$datasource == 1), "P"]) * data1_direction
  full.obj[which(full.obj$datasource == 2), "logP"] <- -log10(full.obj[which(full.obj$datasource == 2), "P"]) * data2_direction

  validate_axis_breaks(axis_breaks, full.obj)

  # x bounds and name
  x_axis_name <- "Chromosome"
  if(is.null(chr) && !is.null(zoom_gene)) {
    hits <- map_df[!is.na(map_df$Gene) & map_df$Gene == zoom_gene, , drop = FALSE]
    chr_candidates <- unique(chr_as_numeric(hits$CHR))
    if (length(chr_candidates) != 1) {
      stop(paste0(zoom_gene, " maps to zero or multiple chromosomes; please pass chr explicitly."))
    }
    chr <- chr_candidates
  } else if(is.null(chr) && !is.null(zoom_ensg)) {
    hits <- map_df[!is.na(map_df$ENSG) & map_df$ENSG == zoom_ensg, , drop = FALSE]
    chr_candidates <- unique(chr_as_numeric(hits$CHR))
    if (length(chr_candidates) != 1) {
      stop(paste0(zoom_ensg, " maps to zero or multiple chromosomes; please pass chr explicitly."))
    }
    chr <- chr_candidates
  }
  if(!is.null(chr)) {
    full.obj <- full.obj[which(full.obj$CHR == chr), ]
    x_axis_name <- paste0("Chromosome ",  chr, " (MB)")
    if(!is.null(zoom_gene)) {
      bounds <- get_gene_bounds(zoom_gene, map_df, chr)
      bounds[1] <- bounds[1] - 250000
      bounds[2] <- bounds[2] + 250000
    } else if(!is.null(zoom_ensg)) {
      bounds <- get_gene_bounds_ensg(zoom_ensg, map_df, chr)
      bounds[1] <- bounds[1] - 250000
      bounds[2] <- bounds[2] + 250000
    } else if(zoom_right != 0) {
      bounds <- c(zoom_left, zoom_right)
    } else {
      bounds <- c(zoom_left, zoom_right)
    }
    full.obj <- full.obj[which(full.obj$BP >= bounds[1] & full.obj$BP <= bounds[2]), ]
  }
  full.obj <- to_megabase(full.obj)
  axis_set <- aggregate(full.obj$absolute, by = list(full.obj[["CHR"]]), mean)
  names(axis_set) <- c("CHR", "center")

  ### simple axis breaks try
  if(is.numeric(axis_breaks)) {
    cat("using manual breaks")
  } else if(identical(axis_breaks, 'auto')) {
    axis_breaks <- c()
    if(min(full.obj$logP) < -30) {
      axis_breaks <- append(axis_breaks, -30)
    }
    if(max(full.obj$logP) > 30) {
      axis_breaks <- append(axis_breaks, 30)
    }
    cat('using auto breaks: ', axis_breaks, "\n")
  } else {
    axis_breaks = NULL
    cat('no break')
  }

  tag1 <- build_tag_subset(full.obj[which(full.obj$datasource == 1), ], tag_names1, tag_threshold1)
  tag2 <- build_tag_subset(full.obj[which(full.obj$datasource == 2), ], tag_names2, tag_threshold2)

  for_tag <- rbind(tag1, tag2)

  # y tick marks
  ymax_pad <- 1
  ymin_pad <- 1
  if(is.null(data2)) {
    ifelse(data1_direction == 1, ymin_pad <- 0, ymax_pad <- 0)
  }

  if(is.na(y_min) == TRUE) {y_min <- round(min(full.obj$logP, na.rm=TRUE) - 1*ymin_pad)} # pad only in data direction
  if(is.na(y_max) == TRUE) {y_max <- round(max(full.obj$logP, na.rm=TRUE) + 1*ymax_pad)} # pad only in data direction
  y_range <- y_max - y_min

  signed_log <- function(x) sign(x) * log1p(abs(x))
  signed_log_inv <- function(x) sign(x) * (expm1(abs(x)))

  if (y_range > 60) {
    # Generate candidate tick positions in signed-log space, then back-transform
    log_min <- signed_log(y_min)
    log_max <- signed_log(y_max)
    log_breaks_raw <- seq(log_min, log_max, length.out = 16)
    log_breaks <- unique(round(signed_log_inv(log_breaks_raw)))
    break_length <- NULL  # not used in log scale mode
    y_log_scale <- TRUE
  } else {
    # Original linear tick mark logic
    if(is.na(y_ticks) == TRUE) { if((y_range) > 16) {y_ticks <- 16} else {y_ticks <- round(y_range)}}
    break_length <- round(y_range / y_ticks)
    y_log_scale <- FALSE
  }
  # plotting
  plt <- ggplot() + theme_bw()
  plt <- plt + plot.mh(data = full.obj[which(full.obj$datasource == 1), ],
                       direction = data1_direction,
                       draw_genes=draw_genes)
  plt <- plt + plot.mh(data = full.obj[which(full.obj$datasource == 2), ],
                       direction = data2_direction,
                       draw_genes=draw_genes)
  only_gwas <- all(full.obj$datatype == "gwas", na.rm = TRUE)
  plt <- plt + scale_color_identity("Tissue", guide = if (only_gwas) "none" else "legend",
                                    labels = full.obj$name,
                                    breaks = full.obj$color)
  plt <- plt + xlab(x_axis_name)
  has_neg_logP <- any(full.obj$logP < 0, na.rm = TRUE)
  has_pos_logP <- any(full.obj$logP > 0, na.rm = TRUE)
  y_axis_label <- if (has_neg_logP && has_pos_logP) {
    expression(log["10"]*italic((p))~"&"~-log["10"]*italic((p)))
  } else if (has_neg_logP) {
    expression(log["10"]*italic((p)))
  } else {
    expression(-log["10"]*italic((p)))
  }
  plt <- plt + ylab(y_axis_label)
  plt <- plt + theme(axis.text.x=element_text(color='black'),
                     axis.text.y=element_text(color='black'),
                     axis.title.x=element_text(face="bold", color="black"),
                     axis.title.y=element_text(face="bold", color="black"),
                     axis.ticks.x=element_line())
  plt <- plt + theme(panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank())
  plt <- plt + ggplot2::guides(shape = "none",
                               size = "none",
                               colour = if (only_gwas) "none" else guide_legend(override.aes = list(size=6)))
  if(!is.null(sig_line1)) {
    plt <- plt + geom_hline(aes(yintercept = data1_direction * -log10(sig_line1)), color = sig_line1_color, linetype='dashed')
  }
  if(!is.null(sig_line2)) {
    plt <- plt + geom_hline(aes(yintercept = data2_direction * -log10(sig_line2)), color = sig_line2_color, linetype='dashed')
  }
  plt <- plt + geom_hline(aes(yintercept = 0), linewidth = 1)
  if(is.null(chr)) {
    plt <- plt + scale_x_continuous(label = axis_set$CHR, breaks = axis_set$center)
  } else {
    plt <- plt + scale_x_continuous(label = waiver(), breaks = waiver())
  }
  if (y_log_scale) {
    plt <- plt + scale_y_continuous(breaks = log_breaks, limits = c(y_min, y_max))
  } else {
    plt <- plt + scale_y_continuous(breaks=seq(y_min, y_max, break_length),
                                    limits=c(y_min, y_max))
  }
  plt <- plt + ggrepel::geom_label_repel(data = for_tag,
                                         fontface = "bold",
                                         aes(x = absolute,
                                             y = logP,
                                             label = ifelse(datatype == "gwas",
                                                            snp.name,
                                                            paste0(gene.name, "-", names))))
  plt <- plt + theme(
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    panel.spacing = ggplot2::unit(0, "pt")
  )

  #use no break
  # --- axis break handling ---------------------------------------------------
  # axis_breaks: FALSE (no break) | numeric vector of length 1 or 2
  # axis_break_scale: user-supplied scale factor (default should be 5 in function signature)
  if (identical(axis_breaks, FALSE) || length(axis_breaks) == 0) {
    # no breaks — nothing to do
  } else if (length(axis_breaks) == 1) {
    # single break
    # below 0 → bottom segment compressed, use scale as-is
    # above 0 → top segment compressed, invert the scale
    if (axis_breaks < 0) {
      plt <- plt + ggbreak::scale_y_break(
        breaks = c(axis_breaks, axis_breaks + 0.01),
        scales = axis_break_scale,
        space  = 0.0, symbol = "slash", expand = c(0,0)
      )
    } else {
      plt <- plt + ggbreak::scale_y_break(
        breaks = c(axis_breaks, axis_breaks + 0.01),
        scales = 1 / axis_break_scale,
        space  = 0.0, symbol = "slash", expand = c(0,0)
      )
    }

  } else if (length(axis_breaks) == 2) {
    bottom_break <- min(axis_breaks)
    top_break    <- max(axis_breaks)

    message(sprintf("bottom break interval: c(%g, %g)", bottom_break - 0.01, bottom_break))
    message(sprintf("top break interval:    c(%g, %g)", top_break, top_break + 0.01))
    message(sprintf("axis_break_scale: %g", axis_break_scale))
    message(sprintf("logP range: [%g, %g]", min(full.obj$logP, na.rm=TRUE), max(full.obj$logP, na.rm=TRUE)))

    plt <- plt + ggbreak::scale_y_break(breaks = c(bottom_break - 0.01, bottom_break),
                                        scales = axis_break_scale,
                                        space = 0.0, symbol = "slash",
                                        expand = c(0,0))
    plt <- plt + ggbreak::scale_y_break(breaks = c(top_break, top_break + 0.01),
                                        scales = 1/axis_break_scale,  # should be 1?
                                        space = 0.0, symbol = "slash",
                                        expand = c(0,0))
  } else {
    stop("axis_breaks must be FALSE, or a numeric vector of length 1 or 2")
  }
  plt
}
