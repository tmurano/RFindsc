#' Add RFind-sc scores to a Seurat object
#'
#' Writes the cell x panel scores of [run_pipeline()], [score_panel_pairs()],
#' [score_genesets()] or [score_cells()] into the cell metadata, one column per
#' panel, so that they can be plotted and grouped like any other metadata
#' (for example `Seurat::FeaturePlot(object, "RFind_DAM")`). This is the
#' counterpart of `Seurat::AddModuleScore()`; the scores themselves are not
#' changed. Cells are matched by name.
#'
#' @param object A Seurat object, or a data frame with one row per cell and the
#'   cell names as row names (for example the `meta.data` of a Seurat object,
#'   or the `colData` of a SingleCellExperiment converted with
#'   `as.data.frame()`).
#' @param scores The list returned by a scoring function (its `scores` element
#'   is used), or a numeric matrix with cell names as row names and panel names
#'   as column names.
#' @param prefix Character: prefix of the new column names. Default
#'   `"RFind_"` (column `RFind_DAM` for the panel `DAM`).
#' @param overwrite Logical: replace columns of the same name that are already
#'   there. Default `TRUE`; with `FALSE`, such columns stop the function.
#'
#' @return `object` with one new column per panel. Cells of `object` that are
#'   not in `scores` get `NA` (with a message).
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' seu <- readRDS("my_seurat.rds")       # raw counts in the counts layer of the RNA assay
#' out <- run_pipeline(seu)              # the 4 bundled Hammond axes
#' seu <- add_rfind_scores(seu, out)     # columns RFind_DAM, RFind_Maturity, ...
#' FeaturePlot(seu, "RFind_DAM")
#' }
#'
#' @export
add_rfind_scores <- function(object, scores, prefix = "RFind_", overwrite = TRUE) {
  if (is.list(scores) && !is.data.frame(scores) && !is.null(scores$scores)) {
    scores <- scores$scores
  }
  scores <- as.matrix(scores)
  if (!is.numeric(scores) || is.null(rownames(scores)) || is.null(colnames(scores))) {
    stop("`scores` must be the result of a scoring function or a numeric matrix ",
         "with cell names as row names and panel names as column names.")
  }
  cols <- paste0(prefix, colnames(scores))
  if (anyDuplicated(cols)) {
    stop("Duplicated column names: ", paste(unique(cols[duplicated(cols)]), collapse = ", "))
  }
  is_seurat <- inherits(object, "Seurat")
  if (is_seurat) {
    if (!requireNamespace("SeuratObject", quietly = TRUE)) {
      stop("SeuratObject is needed to add scores to a Seurat object.")
    }
    cells <- colnames(object)
    present <- colnames(object[[]])
  } else if (is.data.frame(object)) {
    cells <- rownames(object)
    present <- colnames(object)
  } else {
    stop("`object` must be a Seurat object or a data frame with one row per cell.")
  }
  hit <- match(cells, rownames(scores))
  if (all(is.na(hit))) {
    stop("No cell of `object` is in `scores` (the row names of `scores` must be cell names).")
  }
  existing <- intersect(cols, present)
  if (length(existing) && !overwrite) {
    stop("Columns already present: ", paste(existing, collapse = ", "),
         " (use overwrite = TRUE to replace them).")
  }
  if (anyNA(hit)) {
    message(sum(is.na(hit)), " cell(s) of `object` are not in `scores`; their scores are NA.")
  }
  add <- scores[hit, , drop = FALSE]
  dimnames(add) <- list(cells, cols)
  add <- as.data.frame(add)
  if (is_seurat) return(SeuratObject::AddMetaData(object, metadata = add))
  object[cols] <- add
  object
}
