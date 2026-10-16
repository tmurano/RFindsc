#' End-to-end RFind-sc scoring pipeline
#'
#' Top-level entry: accepts a Seurat object (or a raw counts matrix) + a
#' bioset list, runs per-cell signature construction and Running Fisher
#' scoring (four up/down terms, signed), returns a cell x panel score matrix.
#'
#' @param input Either a Seurat object (if the `Seurat` package, >= 5.0.0, is
#'   installed) or a gene x cell counts matrix. From a Seurat object the raw
#'   counts are read from the `counts` layer of the `RNA` assay; if that layer
#'   is split (e.g. after `merge()`), join it first with
#'   `SeuratObject::JoinLayers(input, assay = "RNA")`.
#' @param biosets A bioset list. Default [default_biosets()] (the 4 Hammond axes
#'   of the RFind-sc paper). Combine multiple sources with `c(default_biosets(),
#'   load_dedonno_biosets())`.
#' @param config Configuration from [rfindsc_config()].
#' @param return_signatures Logical. If TRUE, also return per-cell signatures.
#' @inheritParams score_cells
#'
#' @return A list. See [score_cells()].
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' mg <- readRDS("processed/seurat_mg_small.rds")
#' out <- run_pipeline(mg)
#' head(out$scores)
#' }
#'
#' @export
run_pipeline <- function(input,
                         biosets = default_biosets(),
                         config = rfindsc_config(),
                         return_signatures = FALSE,
                         libsize = NULL,
                         baseline_groups = NULL) {
  counts <- extract_counts(input)
  score_cells(counts, biosets = biosets, config = config,
              return_signatures = return_signatures, libsize = libsize,
              baseline_groups = baseline_groups)
}

#' Extract a counts matrix from a Seurat object or passthrough
#' @keywords internal
extract_counts <- function(input) {
  if (inherits(input, "Seurat")) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("Seurat package required to score Seurat objects.")
    }
    # The `layer` argument below needs SeuratObject >= 5.0.0; an older version
    # would not select the raw counts.
    if (utils::packageVersion("SeuratObject") < "5.0.0") {
      stop("Scoring a Seurat object needs SeuratObject >= 5.0.0; ",
           "pass the raw counts matrix instead.")
    }
    return(Seurat::GetAssayData(input, assay = "RNA", layer = "counts"))
  }
  if (is.matrix(input) || inherits(input, "Matrix")) return(input)
  stop("`input` must be a Seurat object or a gene x cell counts matrix.")
}

#' Score cells against single-list gene sets (AUCell-compatible interface)
#'
#' Drop-in replacement for [AUCell::AUCell_calcAUC()] that takes the same
#' input shape (a gene x cell counts matrix or Seurat object, plus a named
#' list of gene character vectors) and returns a cell x set RFind-sc score
#' matrix using the unsigned two-term reduction (DOWN list empty). Each gene
#' set is a membership set with no order (`ranked = FALSE`): it is used whole,
#' and only the cut-offs of the cell's ranked lists are scanned (see
#' [running_fisher()]), so the order of the genes in a set does not matter
#' (0.3.4; before, sets were scanned as if ranked in the order given). Use this
#' when porting an AUCell-based workflow; use [run_pipeline()] or
#' [score_panel_pairs()] when you have UP and DOWN lists and want the signed
#' four-term direction-aware score. To score a ranked UP-only list, pass it to
#' [score_panel_pairs()] without a `down` side.
#'
#' @param input Either a Seurat object (see [run_pipeline()]) or a gene x cell
#'   counts matrix.
#' @param gene_sets A named list of character vectors. Each element is a gene
#'   set (a single UP list, AUCell style; order ignored). The list names become
#'   panel names.
#' @param config Configuration from [rfindsc_config()].
#' @param return_signatures Logical. Default FALSE.
#' @inheritParams score_cells
#'
#' @return A list with `scores` (cell x set numeric matrix) and `universe`
#'   (genes used). See [score_cells()] for full structure.
#'
#' @examples
#' \dontrun{
#' library(RFindsc)
#' # Direct counterpart to AUCell::AUCell_calcAUC():
#' gene_sets <- list(
#'   microglia_marker = c("Tmem119", "P2ry12", "Cx3cr1"),
#'   inflammation     = c("Tnf", "Il1b", "Nfkb1")
#' )
#' mg <- readRDS("processed/seurat_mg_small.rds")
#' out <- score_genesets(mg, gene_sets)
#' head(out$scores)
#' }
#'
#' @seealso [run_pipeline()] for direction-aware scoring with UP/DOWN panels,
#'   [score_panel_pairs()] for explicit (UP, DOWN) panel input.
#'
#' @export
score_genesets <- function(input,
                           gene_sets,
                           config = rfindsc_config(),
                           return_signatures = FALSE,
                           libsize = NULL,
                           baseline_groups = NULL) {
  if (!is.list(gene_sets) || length(gene_sets) == 0L) {
    stop("`gene_sets` must be a non-empty named list of gene vectors.")
  }
  if (is.null(names(gene_sets)) || any(names(gene_sets) == "")) {
    stop("`gene_sets` must have unique non-empty names for each gene set.")
  }
  biosets <- lapply(gene_sets, function(g) {
    g <- as.character(g)
    g <- g[!is.na(g) & nzchar(g)]
    list(
      up   = data.frame(gene = g, diff = rep(NA_real_, length(g)),
                         rank = seq_along(g), stringsAsFactors = FALSE),
      down = data.frame(gene = character(0), diff = numeric(0),
                         rank = integer(0)),
      type = "gene_set",
      ranked = FALSE,          # membership set: used whole (0.3.4)
      label = NA_character_,
      source = NA_character_
    )
  })
  names(biosets) <- names(gene_sets)
  counts <- extract_counts(input)
  score_cells(counts, biosets = biosets, config = config,
              return_signatures = return_signatures, libsize = libsize,
              baseline_groups = baseline_groups)
}

#' Score cells against (UP, DOWN) panel pairs (direction-aware)
#'
#' Convenience wrapper around [run_pipeline()] for explicit `(UP, DOWN)` panel
#' input. Use this when you have DEG-derived case-vs-control panel pairs and
#' want the signed four-term direction-aware score.
#'
#' @param input Either a Seurat object (see [run_pipeline()]) or a gene x cell
#'   counts matrix.
#' @param panels A named list. Each element must contain at least `up` (a
#'   character vector or data.frame with column `gene`) and may contain
#'   `down` (same shape; default empty). The lists are taken as ranked, most
#'   extreme first; add `ranked = FALSE` to a panel whose lists are unordered
#'   membership sets, so that they are used whole (see [running_fisher()]).
#'   An optional `eligible` (the genes the reference's study measured) sets
#'   the pair universe of the panel (see [score_cells()]; 0.3.5).
#' @param config Configuration from [rfindsc_config()].
#' @param return_signatures Logical. Default FALSE.
#' @inheritParams score_cells
#'
#' @return Same as [score_cells()].
#'
#' @examples
#' \dontrun{
#' panels <- list(
#'   my_contrast = list(up = case_up_genes, down = case_down_genes)
#' )
#' out <- score_panel_pairs(mg, panels)
#' }
#'
#' @seealso [score_genesets()] for AUCell-style single-list input.
#'
#' @export
score_panel_pairs <- function(input,
                              panels,
                              config = rfindsc_config(),
                              return_signatures = FALSE,
                              libsize = NULL,
                              baseline_groups = NULL) {
  if (!is.list(panels) || length(panels) == 0L) {
    stop("`panels` must be a non-empty named list.")
  }
  if (is.null(names(panels)) || any(names(panels) == "")) {
    stop("`panels` must have unique non-empty names.")
  }
  to_df <- function(x, sign = 1) {
    if (is.null(x) || length(x) == 0L) {
      return(data.frame(gene = character(0), diff = numeric(0),
                         rank = integer(0)))
    }
    if (is.character(x)) {
      x <- as.character(x); x <- x[!is.na(x) & nzchar(x)]
      return(data.frame(gene = x, diff = rep(sign * 1, length(x)),
                         rank = seq_along(x), stringsAsFactors = FALSE))
    }
    if (is.data.frame(x) && "gene" %in% names(x)) return(x)
    stop("Each panel `up`/`down` must be a character vector or a data.frame with `gene`.")
  }
  biosets <- lapply(panels, function(p) {
    up_df <- to_df(p$up, sign = 1)
    dn_df <- to_df(p$down, sign = -1)
    b <- list(
      up = up_df,
      down = dn_df,
      type = if (nrow(dn_df) > 0L) "deg" else "gene_set",
      label = NA_character_,
      source = NA_character_
    )
    if (!is.null(p$ranked)) b$ranked <- p$ranked
    if (!is.null(p$eligible)) b$eligible <- as.character(p$eligible)
    b
  })
  names(biosets) <- names(panels)
  counts <- extract_counts(input)
  score_cells(counts, biosets = biosets, config = config,
              return_signatures = return_signatures, libsize = libsize,
              baseline_groups = baseline_groups)
}
