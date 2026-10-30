# =============================================================================
# Shared test fixtures for nashvilleplots tests
# Auto-loaded by testthat before any test files run
# =============================================================================

library(ggplot2)
library(ggbreak)

# --- Paths to test data files ------------------------------------------------
# chr 1 data only
GWAS_FILE_1 <- testthat::test_path("fixtures", "HbA1c.CURSMK.COMBINED.JOINT.metal1.miami.txt")
GWAS_FILE_2 <- testthat::test_path("fixtures", "HbA1c.CURSMK.COMBINED.INT.metal1.miami.txt")
META_DIR_1 <- testthat::test_path("fixtures", "meta/")

# --- Shared GWAS read arguments ----------------------------------------------
GWAS_READ_ARGS <- list(
  samp      = 0.01,
  chr       = "chr",
  p         = "P-value",
  bp        = "pos",
  snp.name  = "MarkerName"
)

# --- Lazy-loaded fixtures (only read from disk when first accessed) -----------
#     Wrapped in functions so data isn't loaded if a test file doesn't need it

gwas1_fixture <- function() {
  read_gwas_file(
    file     = GWAS_FILE_1,
    samp     = GWAS_READ_ARGS$samp,
    chr      = GWAS_READ_ARGS$chr,
    p        = GWAS_READ_ARGS$p,
    bp       = GWAS_READ_ARGS$bp,
    snp.name = GWAS_READ_ARGS$snp.name
  )
}

gwas2_fixture <- function() {
  read_gwas_file(
    file     = GWAS_FILE_2,
    samp     = GWAS_READ_ARGS$samp,
    chr      = GWAS_READ_ARGS$chr,
    p        = GWAS_READ_ARGS$p,
    bp       = GWAS_READ_ARGS$bp,
    snp.name = GWAS_READ_ARGS$snp.name
  )
}

meta1_fixture <- function(map = "38", samp = NULL) {
  read_metaXcan_folder(
    directory = META_DIR_1,
    map_df = map,
    samp = samp
  )
}

# --- Default nashville.plot() args -------------------------------------------
DEFAULT_NASH_ARGS <- list(
  data1_direction  = "up",
  map_df           = "38",
  sig_line1        = 1e-9,
  sig_line2        = 1e-9,
  sig_line1_color  = "darkred",
  sig_line2_color  = "darkred",
  axis_breaks      = c(-5, 5)
)

# --- Helper to build a plot with defaults -------
make_nashplot <- function(...) {
  overrides <- list(...)
  args <- utils::modifyList(
    c(list(data1 = gwas1_fixture(), data2 = gwas2_fixture()), DEFAULT_NASH_ARGS),
    overrides
  )
  do.call(nashville.plot, args) + ggplot2::theme(legend.position = "none")
}
