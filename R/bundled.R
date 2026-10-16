#' Load bundled Hammond reference biosets (4 axes of the RFind-sc paper)
#'
#' Returns the 4 named axes used for the Hammond et al. microglia atlas in the
#' RFind-sc paper (Fig. 2):
#' * Signed references (500 up- and 500 down-regulated genes by fold change):
#'   `DAM` (Keren-Shaul 2017), `Maturity` (Matcovitch-Natan 2016; adult vs
#'   embryonic microglia, high = mature).
#' * Gene groups (MSigDB Hallmark membership, UP only): `Interferon`, `G2M`. A
#'   gene group has no order (`ranked = FALSE`): it is used whole at scoring
#'   (0.3.4; before, it was scanned in alphabetical order).
#'
#' The two signed references carry `eligible`, the genes their source study
#' reported (finite, non-zero fold change); [score_cells()] uses it for the pair
#' universe (0.3.5). All mouse gene symbols, in reference rank order. Changed in 0.3.0: the
#' former 9-axis set (with the internal `Development`, `Aging`, `LPC` and the
#' `Immaturity` axis) was replaced. Regenerate with
#' `data-raw/build_bundled_biosets.R`.
#'
#' @return A named list of 4 panels plus `$meta`. Each panel has `$up`, `$down`
#'   data frames (columns `gene`, `diff`, `rank`) plus `label`, `type`
#'   (`"signed"` or `"gene_group"`), `ranked` (`FALSE` for a gene group; added
#'   in 0.3.4), `eligible` (ranked references only; added in 0.3.5) and
#'   `source`.
#' @export
default_biosets <- function() {
  readRDS(system.file("extdata", "biosets.rds", package = "RFindsc",
                      mustWork = TRUE))
}

#' Load bundled De Donno reference biosets (4 axes of the RFind-sc paper)
#'
#' Returns the 4 axes used for the De Donno et al. hippocampus dataset in the
#' RFind-sc paper (Fig. 3): `Glucocorticoid` (corticosterone vs vehicle, ventral
#' hippocampus, GSE184924; 500 up + 500 down), `Seizure` (kainate
#' seizure 6 h vs saline, hippocampal excitatory neurons, GSE125068; top 200
#' up-regulated genes, ranked; type `"ranked_up"`), `UPR` and `NF-kB` (MSigDB
#' Hallmark unfolded protein response and TNF-alpha signalling via NF-kB, UP
#' only; gene groups with `ranked = FALSE`, used whole). Added in 0.3.0. In
#' 0.3.4 the `ranked` field was added and the type of `Seizure` changed from
#' `"gene_group"` to `"ranked_up"`. In 0.3.5 `Glucocorticoid` and `Seizure`
#' gained `eligible` (the genes their source study reported), used for the pair
#' universe.
#'
#' @return A named list of 4 panels plus `$meta`, in the same form as
#'   [default_biosets()].
#' @export
load_dedonno_biosets <- function() {
  readRDS(system.file("extdata", "biosets_dedonno.rds", package = "RFindsc",
                      mustWork = TRUE))
}
