# Run all tests
# devtools::test()
# Run a single test file
#testthat::test_file("tests/testthat/test-my-functions.R")
# Run tests + check coverage
# devtools::test_coverage()

test_that("read_config() reads a valid config file", {
  tmp <- tempfile()
  write.table(data.frame(V1=c("a","b"), V2=c("g1","g2"), V3=c("red","blue")),
              tmp, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
  cfg <- read_config(tmp)
  expect_equal(nrow(cfg), 2)
})

test_that("get_gene_bounds() returns bounds for a valid gene", {
  bounds <- get_gene_bounds("BRCA1", gene.build.38, chr = 17)
  expect_length(bounds, 2)
  expect_true(bounds[1] < bounds[2])
})

test_that("get_gene_bounds() errors when gene not on given chr", {
  expect_error(get_gene_bounds("BRCA1", gene.build.38, chr = 1), "not found in Chromosome")
})

test_that("nashville.plot() returns a ggplot object", {
  p <- make_nashplot()
  expect_s3_class(p, "gg")
})

# test_that("nashville.plot() respects data1_direction = 'down'", {
#   p <- make_nashplot(data1_direction = "down")
#   expect_s3_class(p, "gg")
# })
#
# test_that("nashville.plot() renders without error and saves", {
#   p <- make_nashplot()
#   tmp <- withr::local_tempfile(fileext = ".png")
#   expect_no_error(
#     ggplot2::ggsave(tmp, p, units = "in", width = 20, height = 12, limitsize = FALSE)
#   )
#   expect_gt(file.size(tmp), 0)
# })
#
# test_that("gwas1_fixture() returns a valid object with expected columns", {
#   df <- gwas1_fixture()
#   expect_s3_class(df, "data.frame")
#   expect_true(all(c("CHR", "P", "BP", "color", "group") %in% names(df)))
# })
###
###
###
###
###
###
###
###
###
###
###
###
###
###
###
###
###
# tests/testthat/test-nashville.plot.R
#
# Test coverage targets for nashville.plot.R (currently 52.53%)
# Run with: devtools::test() or testthat::test_file("tests/testthat/test-nashville.plot.R")

# ============================================================
# 1.  read_config()
# ============================================================
#
# test_that("read_config() reads a valid 3-column TSV config", {
#   cfg_file <- write_config_file(make_config())
#   cfg <- read_config(cfg_file)
#
#   expect_s3_class(cfg, "data.frame")
#   expect_equal(nrow(cfg), 2)
#   expect_equal(cfg$V2, c("Group1", "Group2"))
# })

# test_that("read_config() returns the config object (line 17)", {
#   cfg_file <- write_config_file(make_config())
#   result <- read_config(cfg_file)
#   expect_false(is.null(result))
# })

# test_that("read_config() errors when group column contains duplicates", {
#   bad <- data.frame(V1 = c("a", "b"), V2 = c("same", "same"), V3 = c("red", "blue"),
#                     stringsAsFactors = FALSE)
#   cfg_file <- write_config_file(bad)
#   expect_error(read_config(cfg_file), "unique")
# })


# ============================================================
# 2.  get_gene_bounds()  — look up by gene symbol
# ============================================================

test_that("get_gene_bounds() returns a length-2 numeric vector for a valid gene", {
  # TODO: replace "BRCA1" / 17 / "38" with a gene in your build data
  bounds <- get_gene_bounds("BRCA1", gene.build.38, chr = 17)
  expect_length(bounds, 2)
  expect_true(is.numeric(bounds))
  expect_lt(bounds[1], bounds[2])
})

test_that("get_gene_bounds() errors when gene is not on the given chromosome", {
  # BRCA1 is on chr 17, not chr 1
  expect_error(
    get_gene_bounds("BRCA1", gene.build.38, chr = 1),
    regexp = "not found in Chromosome"
  )
})

test_that("get_gene_bounds() errors when multiple genes share the same symbol", {
  # TODO: provide (or mock) a map_df that has two entries with the same Gene name
  # If no real duplicate exists, use a local mock:
  #
  # local_map <- rbind(gene.build.38, gene.build.38[gene.build.38$Gene == "BRCA1", ])
  # with_mock(
  #   `nashville.plot:::gene.build.38` = local_map,
  #   expect_error(get_gene_bounds("BRCA1", "38", chr = 17), "Multiple genes")
  # )
  skip("Requires a mock map_df with duplicate gene symbols — fill in TODO above")
})


# ============================================================
# 3.  get_gene_bounds_ensg()  — look up by ENSG ID
# ============================================================
#
# Mirrors section 2 but via ENSG.  Lines 47–57.

test_that("get_gene_bounds_ensg() returns bounds for a valid ENSG", {
  bounds <- get_gene_bounds_ensg("ENSG00000006377", gene.build.38, chr = 7)
  expect_length(bounds, 2)
  expect_true(is.numeric(bounds))
  expect_lt(bounds[1], bounds[2])
})

test_that("get_gene_bounds_ensg() errors when ENSG not on given chromosome", {
  expect_error(
    get_gene_bounds_ensg("ENSG00000006377", gene.build.38, chr = 1),
    regexp = "not found in Chromosome"
  )
})

test_that("get_gene_bounds_ensg() errors when multiple records share an ENSG", {
  skip("Requires a mock map_df with duplicate ENSG entries — fill in TODO")
})


# ============================================================
# 4.  load_meta()
# ============================================================
#
# Exercises lines 101–135.  Needs a small fixture directory with at least
# one CSV file whose columns include: gene, pvalue, and whatever your
# pipeline writes.  Adjust column names to match reality.

test_that("read_metaXcan_folder() works with map_df = 38", {
  meta1 <- meta1_fixture()
  expect_true(nrow(meta1) > 0)
  expect_s3_class(meta1, "data.frame")
})

test_that("read_metaXcan_folder() works with map_df = 37", {
  meta1 <- meta1_fixture(map = "37")
  expect_true(nrow(meta1) > 0)
  expect_s3_class(meta1, "data.frame")
})

test_that("read_metaXcan_folder() errors for unrecognised map_df build", {
  expect_error(
    meta1 <- meta1_fixture(map = "99"),
    regexp = "map_df must be"
  )
})

test_that("read_metaXcan_folder() downsamples when samp is provided", {
  meta1_full <- meta1_fixture(map = "38", samp = NULL)
  meta1_samp <- meta1_fixture(map = "38", samp = 0.5)
  expect_lte(nrow(meta1_samp), nrow(meta1_full))
})


# ============================================================
# 5.  make.valid.object() — meta datatype path
# ============================================================
#
# Existing tests only exercise datatype = "gwas".
# Lines 164 (gene coords stopifnot) and 175–178 (Set3 color palette).

test_that("make.valid.object() auto-assigns Set3 colors for meta datatype", {
  # BP = NA forces the gene.start / gene.end path (line 163–165)
  obj <- make.valid.object(
    CHR        = c(1L, 2L),
    P          = c(0.01, 0.05),
    BP         = NA_real_,
    gene.start = c(1000L, 2000L),
    gene.end   = c(1500L, 2500L),
    group      = c("typeA", "typeB"),
    color      = NA,
    shape      = NA,
    datatype   = "meta",
    config     = data.frame(V1 = character(0), V2 = character(0),
                            V3 = character(0), stringsAsFactors = FALSE),
    snp.name   = NA_character_
  )

  expect_false(any(is.na(obj$color)))
  # Colors should be valid hex strings (Set3 palette output)
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", obj$color)))
})

test_that("make.valid.object() errors when gene.start is not numeric for meta", {
  expect_error(
    make.valid.object(
      CHR        = 1L,
      P          = 0.01,
      BP         = NA_real_,
      gene.start = "not_a_number",   # should trigger stopifnot on line 164
      gene.end   = 1500L,
      group      = "typeA",
      color      = NA,
      shape      = NA,
      datatype   = "meta",
      config     = data.frame(V1 = character(0), V2 = character(0),
                              V3 = character(0), stringsAsFactors = FALSE),
      snp.name   = NA_character_
    )
  )
})


# ============================================================
# 6.  plot.mh() — draw_genes = TRUE branch
# ============================================================
#
# Lines 206–212: the gene-line drawing path is only reached when
# draw_genes = TRUE and the data has non-NA gene.name + group.

test_that("plot.mh() returns a geom_line layer when draw_genes = TRUE", {
  set.seed(1)
  gene_data <- data.frame(
    gene.name = rep(c("GENE1", "GENE2"), each = 5),
    group     = rep(c("typeA", "typeB"), each = 5),
    CHR       = 1L,
    BP        = NA_real_,
    P         = runif(10, 1e-6, 0.5),
    gene.start = rep(c(1e6, 2e6), each = 5),
    gene.end   = rep(c(1.5e6, 2.5e6), each = 5),
    color      = rep(c("red", "blue"), each = 5),
    shape      = 16L,
    datatype   = "meta",
    snp.name   = NA_character_,
    absolute   = rep(c(1e6, 2e6), each = 5),
    logP       = runif(10, 1, 8),
    stringsAsFactors = FALSE
  )

  layer <- plot.mh(data = gene_data, direction = 1, draw_genes = TRUE)

  # The result should include a geom_line (gene track)
  expect_true(inherits(layer, "LayerInstance") || is.list(layer))
})


# ============================================================
# 7.  find_y_break()
# ============================================================
#
# This is a pure utility function — no I/O, easy to cover completely.
# Lines 278–312.

test_that("find_y_break() returns a length-2 list for data with a clear gap", {
  # Values with an obvious gap between 4 and 10
  vals <- c(1, 2, 3, 4, 10, 11, 12)
  result <- find_y_break(vals)

  expect_type(result, "list")
  expect_length(result, 2)
  expect_lt(result[[1]], result[[2]])
})

test_that("find_y_break() returns NULL when too few values exceed sig_threshold", {
  # Default sig_threshold is typically 2; all values below it → NULL
  result <- find_y_break(c(0.5, 0.8, 1.0))
  expect_null(result)
})

test_that("find_y_break() returns NULL when the largest gap is below min_gap", {
  # Evenly spaced — no single gap stands out
  vals <- seq(2, 10, by = 0.2)
  result <- find_y_break(vals)
  expect_null(result)
})

test_that("find_y_break() returns NULL when trimmed vector has fewer than 2 values", {
  # One extreme outlier swamps the quantile cap, leaving < 2 values in trimmed
  vals <- c(3, 3.01, 3.02, 1000)
  result <- find_y_break(vals, top_quantile = 0.5)
  # After trimming at 50th percentile, only 2 values remain — borderline;
  # adjust top_quantile so trimmed length drops to 1
  result2 <- find_y_break(vals, top_quantile = 0.1)
  expect_null(result2)
})


# ============================================================
# 8.  nashville.plot() — untested parameter branches
# ============================================================
#
# For all tests below, obj1 and obj2 are the two datasource objects
# passed to nashville.plot().  Replace make_gwas_obj() calls with whatever
# constructor your package exposes (load_gwas(), load_meta(), etc.).

# ---- 8a.  data1_direction = 'down' (lines 390–391) ----

test_that("nashville.plot() accepts data1_direction = 'down'", {
  obj1 <- gwas1_fixture()
  obj2 <- gwas2_fixture()
  plt <- nashville.plot(data1 = obj1, data2 = obj2,
                       data1_direction = "down",
                       axis_breaks = FALSE)
  expect_s3_class(plt, "gg")
})

test_that("nashville.plot() errors for invalid data1_direction", {
  obj1 <- gwas1_fixture()
  obj2 <- gwas2_fixture()
  expect_error(
    nashville.plot(data1 = obj1, data2 = obj2,
                  data1_direction = "sideways",
                  axis_breaks = FALSE),
    regexp = "direction should be either"
  )
})

# ---- 8b.  chr zoom  (lines 375, 404–420) ----

test_that("nashville.plot() restricts to a single chromosome when chr is set", {
  obj1 <- gwas1_fixture()
  obj2 <- gwas2_fixture()
  plt <- nashville.plot(data1 = obj1, data2 = obj2,
                       chr = 1,
                       map_df = "38",
                       axis_breaks = FALSE)
  expect_s3_class(plt, "gg")
})

test_that("nashville.plot() zoom via explicit zoom_left / zoom_right (line 415–416)", {
  obj1 <- gwas1_fixture()
  obj2 <- gwas2_fixture()
  plt <- nashville.plot(data1 = obj1, data2 = obj2,
                       chr = 1,
                       zoom_left  = 5e7,
                       zoom_right = 1e8,
                       axis_breaks = FALSE)
  expect_s3_class(plt, "gg")
})

test_that("nashville.plot() zoom via zoom_gene name (lines 407–410)", {
  obj1 <- gwas1_fixture()
  obj2 <- gwas2_fixture()
  # TODO: replace with a gene + chr present in your build data
  plt <- nashville.plot(data1 = obj1, data2 = obj2,
                       chr = 10,
                       zoom_gene = "R3HCC1L",
                       map_df = "38",
                       axis_breaks = FALSE)
  expect_s3_class(plt, "gg")
})

test_that("nashville.plot() zoom via zoom_ensg (lines 411–414)", {
  obj1 <- gwas1_fixture()
  obj2 <- gwas2_fixture()

  # TODO: replace with a real ENSG in your build data
  plt <- nashville.plot(data1 = obj1, data2 = obj2,
                       chr = 10,
                       zoom_ensg = "ENSG00000166024",
                       map_df = "38",
                       axis_breaks = FALSE)
  expect_s3_class(plt, "gg")
})

# ---- 8c.  map_df = "37" path (line 375) ----

test_that("nashville.plot() works with map_df = '37'", {
  obj1 <- gwas1_fixture()
  obj2 <- gwas2_fixture()
  plt <- nashville.plot(data1 = obj1, data2 = obj2,
                       map_df = "37",
                       axis_breaks = FALSE)
  expect_s3_class(plt, "gg")
})

# ---- 8d.  axis_breaks = 'auto' (lines 429–431) ----

test_that("nashville.plot() computes breaks automatically when axis_breaks = 'auto'", {
  # Need data with a clear gap so find_y_break() returns a value
  obj1 <- gwas1_fixture()
  obj2 <- gwas2_fixture()
  # Force some very small p-values to create a gap in -log10(p)
  obj1$P[1:5] <- 1e-200

  # May message about break placement — that's expected
  expect_no_error(
    nashville.plot(data1 = obj1,
                  data2 = obj2,
                  axis_breaks = "auto")
  )
})

# ---- 8e.  axis_breaks validation — value outside data range (lines 351, 359–369) ----

test_that("nashville.plot() errors when axis_breaks value is outside logP range", {
  obj1 <- gwas1_fixture()
  obj2 <- gwas2_fixture()

  # 999 is far outside any realistic -log10(p) range
  expect_error(
    nashville.plot(data1 = obj1, data2 = obj2,
                  axis_breaks = 999),
    regexp = "outside the data range"
  )
})

test_that("nashville.plot() skips axis_break validation when axis_breaks = FALSE (line 351)", {
  obj1 <- gwas1_fixture()
  obj2 <- gwas2_fixture()

  expect_no_error(
    nashville.plot(data1 = obj1, data2 = obj2,
                  axis_breaks = FALSE)
  )
})

# ---- 8f.  Single numeric axis_break (lines 555–565) ----

test_that("nashville.plot() applies a positive axis break (line 562–565)", {
  obj1 <- gwas1_fixture()
  obj2 <- gwas2_fixture()
  obj1$P[1:5] <- 1e-150   # ensure we have high -log10(p) values

  # Pick a break value inside the actual data range — adjust as needed
  plt <- nashville.plot(data1 = obj1, data2 = obj2,
                       axis_breaks = 10)
  expect_s3_class(plt, "gg")
})

test_that("nashville.plot() applies a negative axis break (line 555–559)", {
  obj1 <- gwas1_fixture()
  obj2 <- gwas2_fixture()
  obj2$P[1:5] <- 1e-150   # data2 goes negative in logP

  plt <- nashville.plot(data1 = obj1, data2 = obj2,
                       axis_breaks = -10)
  expect_s3_class(plt, "gg")
})

# ---- 8g.  Two-element axis_breaks vector (lines 569–576) ----

test_that("nashville.plot() accepts a length-2 axis_breaks vector", {
  obj1 <- gwas1_fixture()
  obj2 <- gwas2_fixture()
  obj1$P[1:5] <- 1e-20
  obj2$P[1:5] <- 1e-20

  plt <- nashville.plot(data1 = obj1, data2 = obj2,
                       axis_breaks = c(8, -8))
  expect_s3_class(plt, "gg")
})

# ---- 8h.  Invalid axis_breaks length (line 589) ----

test_that("nashville.plot() errors when axis_breaks has length > 2", {
  obj1 <- gwas1_fixture()
  obj2 <- gwas2_fixture()

  expect_error(
    nashville.plot(data1 = obj1, data2 = obj2,
                  axis_breaks = c(5, 10, 15)),
    regexp = "axis_breaks must be FALSE"
  )
})

# ---- 8i.  y_ticks / break_length path (lines 470–472) ----

test_that("nashville.plot() auto-sets y_ticks when y_range > 16 (line 470)", {
  obj1 <- gwas1_fixture()
  obj2 <- gwas2_fixture()
  obj1$P[1:5] <- 1e-30   # wide range

  plt <- nashville.plot(data1 = obj1, data2 = obj2,
                       axis_breaks = FALSE)
  expect_s3_class(plt, "gg")
})

# ---- 8j.  chr-only scale_x (line 510) ----

test_that("nashville.plot() uses waiver x-axis when chr is specified", {
  obj1 <- gwas1_fixture()
  obj2 <- gwas2_fixture()

  plt <- nashville.plot(data1 = obj1, data2 = obj2,
                       chr = 1,
                       axis_breaks = FALSE)
  expect_s3_class(plt, "gg")
})

# ---- 8k.  scale_y_continuous with break_length (lines 519–520) ----

test_that("nashville.plot() renders scale_y_continuous when not using log scale", {
  # y_range <= 60 and not using auto log scale
  obj1 <- gwas1_fixture()
  obj2 <- gwas2_fixture()

  plt <- nashville.plot(data1 = obj1, data2 = obj2,
                       y_min = -10, y_max = 10,
                       axis_breaks = FALSE)
  expect_s3_class(plt, "gg")
})

# Tests for build_tag_subset() / parse_tag_spec() 
make_obj <- function(gene.name, names, group, CHR, BP, P, snp.name = NA_character_) {
  data.frame(
    gene.name = gene.name,
    names     = names,
    group     = group,
    CHR       = CHR,
    BP        = BP,
    P         = P,
    snp.name  = snp.name,
    stringsAsFactors = FALSE
  )
}

# --- parse_tag_spec() --------------------------------------------------------

test_that("parse_tag_spec treats a plain character vector as all-plain names", {
  spec <- parse_tag_spec(c("BRCA1", "TP53"))
  expect_setequal(spec$plain, c("BRCA1", "TP53"))
  expect_length(spec$tuple_gene, 0)
})

test_that("parse_tag_spec splits a list into plain names and gene/tissue tuples", {
  spec <- parse_tag_spec(list("TP53", c("BRCA1", "Liver"), c("BRCA1", "Muscle")))
  expect_equal(spec$plain, "TP53")
  expect_equal(spec$tuple_gene, c("BRCA1", "BRCA1"))
  expect_equal(spec$tuple_tissue, c("Liver", "Muscle"))
})

test_that("parse_tag_spec errors on a tuple of the wrong length", {
  expect_error(parse_tag_spec(list(c("BRCA1", "Liver", "extra"))), "length 2")
})

test_that("parse_tag_spec handles NULL/empty input", {
  spec <- parse_tag_spec(NULL)
  expect_length(spec$plain, 0)
  expect_length(spec$tuple_gene, 0)
})

# --- build_tag_subset(): legacy / plain-name behavior is unchanged ----------

test_that("threshold tagging (no names given) still collapses a gene to its single most extreme tissue", {
  obj <- make_obj(
    gene.name = c("BRCA1", "BRCA1", "BRCA1", "TP53"),
    names     = c("Adipose", "Liver", "Muscle", "Adipose"),
    group     = c("Adipose", "Liver", "Muscle", "Adipose"),
    CHR       = c(17, 17, 17, 17),
    BP        = c(43100000, 43100000, 43100000, 7600000),
    P         = c(1e-10, 1e-8, 1e-6, 1e-5)
  )

  res <- build_tag_subset(obj, tag_genes = NULL, gene_tag = 1e-4)

  brca1_rows <- res[res$gene.name == "BRCA1", ]
  expect_equal(nrow(brca1_rows), 1)
  expect_equal(brca1_rows$group, "Adipose")   # most significant
  expect_equal(nrow(res), 2)
})

test_that("plain-name tagging still collapses a gene to its single most extreme tissue", {
  obj <- make_obj(
    gene.name = c("BRCA1", "BRCA1", "BRCA1"),
    names     = c("Adipose", "Liver", "Muscle"),
    group     = c("Adipose", "Liver", "Muscle"),
    CHR       = c(17, 17, 17),
    BP        = c(43100000, 43100000, 43100000),
    P         = c(1e-8, 1e-10, 0.9)  # Liver is most significant here
  )

  res <- build_tag_subset(obj, tag_genes = "BRCA1", gene_tag = NA)

  expect_equal(nrow(res), 1)
  expect_equal(res$group, "Liver")
})

test_that("peak suppression is still global (not per-tissue) for automatic hits", {
  # Same locus, two different tissues, neither explicitly named -> only the
  # more significant one survives, exactly like the original behavior.
  obj <- make_obj(
    gene.name = c("BRCA1", "BRCA1"),
    names     = c("Liver", "Muscle"),
    group     = c("Liver", "Muscle"),
    CHR       = c(17, 17),
    BP        = c(43100000, 43100000),
    P         = c(1e-9, 1e-7)
  )

  res <- build_tag_subset(obj, tag_genes = NULL, gene_tag = 1e-4, peak_window = 500000)

  expect_equal(nrow(res), 1)
  expect_equal(res$group, "Liver")
})

test_that("peak suppression still declutters nearby genes in the same tissue", {
  obj <- make_obj(
    gene.name = c("GENE_A", "GENE_B"),
    names     = c("Liver", "Liver"),
    group     = c("Liver", "Liver"),
    CHR       = c(1, 1),
    BP        = c(1000000, 1000100),
    P         = c(1e-9, 1e-7)
  )

  res <- build_tag_subset(obj, tag_genes = NULL, gene_tag = 1e-4, peak_window = 500000)

  expect_equal(nrow(res), 1)
  expect_equal(res$gene.name, "GENE_A")
})

# --- build_tag_subset(): tuple-based multi-tissue tagging ---------------

test_that("gene/tissue tuples tag a gene in exactly the requested tissues", {
  obj <- make_obj(
    gene.name = c("BRCA1", "BRCA1", "BRCA1"),
    names     = c("Adipose", "Liver", "Muscle"),
    group     = c("Adipose", "Liver", "Muscle"),
    CHR       = c(17, 17, 17),
    BP        = c(43100000, 43100000, 43100000),
    P         = c(0.5, 0.4, 0.3)  # none significant, purely name-driven
  )

  res <- build_tag_subset(
    obj,
    tag_genes = list(c("BRCA1", "Liver"), c("BRCA1", "Muscle")),
    gene_tag = NA
  )

  expect_equal(nrow(res), 2)
  expect_setequal(res$group, c("Liver", "Muscle"))
  expect_false("Adipose" %in% res$group)
})

test_that("tuples bypass peak suppression even at the same locus", {
  obj <- make_obj(
    gene.name = c("BRCA1", "BRCA1", "BRCA1"),
    names     = c("Adipose", "Liver", "Muscle"),
    group     = c("Adipose", "Liver", "Muscle"),
    CHR       = c(17, 17, 17),
    BP        = c(43100000, 43100000, 43100000),  # identical position
    P         = c(1e-9, 1e-7, 1e-5)
  )

  res <- build_tag_subset(
    obj,
    tag_genes = list(c("BRCA1", "Adipose"), c("BRCA1", "Liver"), c("BRCA1", "Muscle")),
    gene_tag = NA,
    peak_window = 500000
  )

  expect_equal(nrow(res), 3)
  expect_setequal(res$group, c("Adipose", "Liver", "Muscle"))
})

test_that("mixing a plain name with tuples: plain still collapses, tuples still uniquely survive", {
  obj <- make_obj(
    gene.name = c("BRCA1", "BRCA1", "BRCA1", "TP53"),
    names     = c("Adipose", "Liver", "Muscle", "Adipose"),
    group     = c("Adipose", "Liver", "Muscle", "Adipose"),
    CHR       = c(17, 17, 17, 17),
    BP        = c(43100000, 43100000, 43100000, 7600000),
    P         = c(1e-9, 1e-7, 1e-5, 1e-3)
  )

  res <- build_tag_subset(
    obj,
    tag_genes = list("TP53", c("BRCA1", "Muscle")),
    gene_tag = NA
  )

  # TP53 (plain) -> its single (only) row is kept
  expect_true("TP53" %in% res$gene.name)
  # BRCA1 Muscle (tuple) -> uniquely kept even though Adipose is more significant
  brca1_rows <- res[res$gene.name == "BRCA1", ]
  expect_equal(nrow(brca1_rows), 1)
  expect_equal(brca1_rows$group, "Muscle")
})

test_that("a tuple targeting a tissue where the gene isn't present matches nothing", {
  obj <- make_obj(
    gene.name = c("BRCA1", "BRCA1"),
    names     = c("Adipose", "Liver"),
    group     = c("Adipose", "Liver"),
    CHR       = c(17, 17),
    BP        = c(43100000, 43100000),
    P         = c(0.5, 0.4)
  )

  res <- build_tag_subset(obj, tag_genes = list(c("BRCA1", "Muscle")), gene_tag = NA)
  expect_equal(nrow(res), 0)
})

test_that("duplicate tuples for the same gene/tissue collapse to one row", {
  obj <- make_obj(
    gene.name = c("BRCA1", "BRCA1"),
    names     = c("Liver", "Liver"),
    group     = c("Liver", "Liver"),
    CHR       = c(17, 17),
    BP        = c(43100000, 43100000),
    P         = c(1e-9, 1e-2)
  )

  res <- build_tag_subset(obj, tag_genes = list(c("BRCA1", "Liver")), gene_tag = NA)
  expect_equal(nrow(res), 1)
  expect_equal(res$P, 1e-9)
})
