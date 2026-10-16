#' RFindsc: direction-aware panel scoring for single-cell RNA-seq
#'
#' @description
#' RFind-sc scores each cell against case-vs-control gene panels with the
#' Running Fisher test (smallest hypergeometric p over every top-G x top-K
#' cut-off pair; four up/down pairs combined into one signed score). Per-cell
#' UP/DOWN signatures are ranked by the P2 Pearson residual
#' `(obs - exp)/sqrt(exp + 1)` of observed vs baseline-expected counts. Panels
#' consist of paired UP and DOWN gene lists (or one gene group); scores are
#' signed (positive = case-aligned, negative = control-aligned). The package
#' ships the reference axes of the RFind-sc paper (4 Hammond axes, 4 De Donno
#' axes).
#'
#' @section Typical workflow:
#' ```r
#' library(RFindsc)
#' library(Seurat)
#' mg <- readRDS("my_microglia_seurat.rds")
#' out <- run_pipeline(mg, biosets = default_biosets())
#' head(out$scores)
#' ```
#'
#' @section Key functions:
#' * [run_pipeline()]: end-to-end wrapper (Seurat or counts matrix in).
#' * [score_cells()]: lower-level scoring.
#' * [running_fisher()]: single-term kernel ([bid_running_fisher()] = former name).
#' * [build_panel_from_deg()]: construct a signed panel from a DEG table.
#' * [default_biosets()], [load_dedonno_biosets()]: bundled references.
#' * [msigdbr_to_panels()], [c2_cgp_auto_pair()]: pull MSigDB collections.
#' * [rfindsc_config()]: hyperparameter defaults.
#'
#' @keywords internal
"_PACKAGE"
