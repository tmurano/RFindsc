#' Build a Type A panel from a DEG table
#'
#' Constructs a bidirectional panel (paired UP/DOWN with log2FC) from a DEG
#' table with gene symbol and log2 fold change columns, ordered by fold change.
#' By default all genes are kept: [score_cells()] keeps the top
#' `config$top_n` per side (500) and then the genes of its pair universe, as in
#' the RFind-sc paper (0.3.8; before, the restriction came first). Set `top_n`
#' to cap the panel itself.
#'
#' The panel also records `eligible`: every gene of the table with a finite,
#' non-zero log2FC (before `padj_thresh` and `top_n` are applied), i.e. the
#' genes the study measured. [score_cells()] uses the genes of the dataset that
#' are also in `eligible` as the background of the comparison (the pair
#' universe), so genes the reference could never report do not count. Pass the
#' full DEG table (all tested genes, not only the significant ones) for this to
#' hold; set `panel$eligible <- NULL` to use all genes of the dataset instead.
#' Added in 0.3.5.
#'
#' @param deg A data frame with at least a gene symbol column and a log2FC
#'   column. NA values in either are dropped.
#' @param gene_col Character: name of the gene column. Default `"gene"`.
#' @param log2fc_col Character: name of the log2 fold-change column.
#'   Default `"log2FC"`.
#' @param top_n Integer or Inf: cap per side. Default Inf (no cap here;
#'   [score_cells()] keeps the top `config$top_n` per side, 500, and then the
#'   genes of the dataset, as in the RFind-sc paper; 1000 before 0.3.0).
#' @param label Character: panel label (stored as metadata).
#' @param source Character: free-form source citation. Default NA.
#' @param padj_col Character or NULL: optional adjusted p-value column.
#'   If non-NULL and `padj_thresh` is numeric, rows with padj > threshold are
#'   filtered out.
#' @param padj_thresh Numeric or NULL.
#'
#' @return A list with elements `up`, `down` (each a data.frame with columns
#'   `gene`, `diff`, `rank`), `eligible` (character vector; 0.3.5),
#'   `type = "deg"`, `label`, `source`.
#' @export
build_panel_from_deg <- function(deg,
                                 gene_col = "gene",
                                 log2fc_col = "log2FC",
                                 top_n = Inf,
                                 label = "user_panel",
                                 source = NA_character_,
                                 padj_col = NULL,
                                 padj_thresh = NULL) {
  stopifnot(is.data.frame(deg),
            gene_col %in% colnames(deg),
            log2fc_col %in% colnames(deg))
  d <- data.frame(
    gene   = as.character(deg[[gene_col]]),
    log2fc = suppressWarnings(as.numeric(deg[[log2fc_col]])),
    stringsAsFactors = FALSE
  )
  d <- d[!is.na(d$gene) & d$gene != "" & !is.na(d$log2fc), ]
  eligible <- unique(d$gene[is.finite(d$log2fc) & d$log2fc != 0])   # genes the study measured
  if (!is.null(padj_col) && !is.null(padj_thresh) &&
      padj_col %in% colnames(deg)) {
    padj <- suppressWarnings(as.numeric(deg[[padj_col]]))
    padj <- padj[!is.na(deg[[gene_col]]) & deg[[gene_col]] != "" &
                 !is.na(as.numeric(deg[[log2fc_col]]))]
    d <- d[is.na(padj) | padj <= padj_thresh, ]
  }
  up <- d[d$log2fc > 0, ]; up <- up[order(-up$log2fc), ]
  up <- utils::head(up, top_n)
  up_df <- data.frame(gene = up$gene, diff = up$log2fc,
                      rank = seq_len(nrow(up)), stringsAsFactors = FALSE)
  dn <- d[d$log2fc < 0, ]; dn <- dn[order(dn$log2fc), ]
  dn <- utils::head(dn, top_n)
  dn_df <- data.frame(gene = dn$gene, diff = dn$log2fc,
                      rank = seq_len(nrow(dn)), stringsAsFactors = FALSE)
  list(up = up_df, down = dn_df, eligible = eligible, type = "deg",
       label = label, source = source)
}
