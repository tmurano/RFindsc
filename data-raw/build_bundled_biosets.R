# ==============================================================================
# data-raw/build_bundled_biosets.R
# ------------------------------------------------------------------------------
# Regenerate the bundled reference biosets from the references of the RFind-sc
# paper (0.3.0, 2026-09-15; `ranked` field added in 0.3.4, 2026-09-17; `eligible` added in 0.3.5, 2026-09-18;
# only the axes the paper uses since 0.3.6, 2026-09-29):
#   inst/extdata/biosets.rds          <- data/Hammond2019/output/references/references_panel9.rds
#       (Hammond 4 axes: DAM, Maturity, Interferon, G2M; other entries of the
#        file are not bundled; built by analysis/hammond/06_build_references.R;
#        the file name is kept from the earlier 9-axis version)
#   inst/extdata/biosets_dedonno.rds  <- data/DeDonno2026/output/references_dedonno.rds
#       (De Donno 4 axes: Glucocorticoid [CORT], Seizure [KA6h], UPR, NF-kB;
#        built by analysis/dedonno/03_references.R)
# Both source paths are relative to the root of the RFind-sc paper repository
# (https://github.com/tmurano/RFind-sc). The Hammond file is tracked there; the
# De Donno file is not, so run analysis/dedonno/03_references.R first (its
# header lists its inputs).
# Not run at package install time. Run with Rscript from the paper repository
# root (or set RFIND_PROJECT_ROOT to it); the bundles are written to
# inst/extdata of the package that contains this script:
#
#   Rscript <path to RFindsc>/data-raw/build_bundled_biosets.R
# ==============================================================================
PROJECT_ROOT <- Sys.getenv("RFIND_PROJECT_ROOT", ".")
SCRIPT <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(SCRIPT) != 1L) stop("run with Rscript: Rscript <path to RFindsc>/data-raw/build_bundled_biosets.R")
SCRIPT <- gsub("~+~", " ", SCRIPT, fixed = TRUE)   # Rscript passes spaces in the path as "~+~"
PKG_EXTDATA <- file.path(dirname(dirname(normalizePath(SCRIPT))), "inst", "extdata")
SRC_H <- file.path(PROJECT_ROOT, "data/Hammond2019/output/references/references_panel9.rds")
SRC_D <- file.path(PROJECT_ROOT, "data/DeDonno2026/output/references_dedonno.rds")
for (f in c(SRC_H, SRC_D)) if (!file.exists(f)) stop("reference file not found: ", f)

# The references store genes in rank order (no fold changes). RFind-sc scoring
# (to_b2) uses only the order of the gene column, so `diff` is left NA.
# Type C references are MSigDB Hallmark membership sets: they have no order (the
# genes are stored alphabetically), so they are bundled with `ranked = FALSE` and
# scored whole (0.3.4). Type B references are ranked lists: `signed` (UP + DOWN)
# or `ranked_up` (UP only, e.g. Seizure).
to_side_df <- function(genes) {
  genes <- as.character(genes)
  genes <- genes[!is.na(genes) & nzchar(genes)]
  data.frame(gene = genes, diff = rep(NA_real_, length(genes)),
             rank = seq_along(genes), stringsAsFactors = FALSE)
}
bundle <- function(ref, rename, source_label) {
  keys <- names(rename)
  bs <- lapply(keys, function(k) {
    e <- ref[[k]]
    ranked <- !identical(e$type, "C")
    p <- list(up = to_side_df(e$up), down = to_side_df(e$down), label = rename[[k]],
              type = if (!ranked) "gene_group" else if (length(e$down)) "signed" else "ranked_up",
              ranked = ranked,
              source = if (!is.null(e$desc)) e$desc else NA_character_)
    # eligible = genes the reference's study reported (finite, non-zero change; the pipeline's references carry it since
    # 2026-09-18). score_cells() intersects it with the dataset's genes: the pair universe. Gene groups have none.
    if (!is.null(e$eligible)) p$eligible <- sort(unique(as.character(e$eligible)))
    p
  })
  names(bs) <- unname(rename)
  bs$meta <- list(source = source_label, axes = unname(rename),
                  types = vapply(bs[unname(rename)], `[[`, "", "type"),
                  ranked = vapply(bs[unname(rename)], `[[`, NA, "ranked"),
                  delta = "P2 Pearson residual (obs - exp)/sqrt(exp + 1), top/bottom 500",
                  kernel = paste("Running Fisher, smallest p over every top-G x top-K cut-off pair;",
                                 "gene groups (ranked = FALSE) are used whole, over the top-G cut-offs only"),
                  universe = paste("pair universe: genes of the dataset in `eligible` (ranked references),",
                                   "all genes of the dataset (gene groups); the cell's lists restricted to it before the top 500,",
                                   "the reference's top 500 per side restricted to it (0.3.8); gene groups used whole"),
                  organism = "mouse", processed_at = format(Sys.time(), "%Y-%m-%d"))
  bs
}

H <- readRDS(SRC_H)
hammond <- bundle(H, c(DAM = "DAM", Maturity = "Maturity", Interferon = "Interferon", G2M = "G2M"),
                  "data/Hammond2019/output/references/references_panel9.rds (RFind-sc paper, Hammond 4 axes)")
D <- readRDS(SRC_D)
dedonno <- bundle(D, c(CORT = "Glucocorticoid", KA6h = "Seizure", UPR = "UPR", NFkB = "NF-kB"),
                  "data/DeDonno2026/output/references_dedonno.rds (RFind-sc paper, De Donno 4 axes)")

dir.create(PKG_EXTDATA, showWarnings = FALSE, recursive = TRUE)
saveRDS(hammond, file.path(PKG_EXTDATA, "biosets.rds"))
saveRDS(dedonno, file.path(PKG_EXTDATA, "biosets_dedonno.rds"))
for (x in list(hammond, dedonno))
  message(sprintf("%s: %s", x$meta$source, paste(sprintf("%s (%d/%d)", x$meta$axes,
    vapply(x$meta$axes, function(a) nrow(x[[a]]$up), 0L), vapply(x$meta$axes, function(a) nrow(x[[a]]$down), 0L)), collapse = ", ")))
