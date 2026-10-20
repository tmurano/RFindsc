#' Build per-cell UP/DOWN signature from Pearson residuals (Delta = P2)
#'
#' For one cell, ranks genes by the Pearson residual of observed vs expected
#' counts, `Delta = (obs - exp) / sqrt(exp + 1)`, where `exp` is the count the
#' gene would receive under the baseline rate (dataset-wide, or the cell's own
#' group with `baseline_groups` in [score_cells()]) scaled to this cell's
#' library size. Returns (1) the top `n_cell` genes with the largest
#' positive residual (`up_idx`, over-expressed relative to baseline) and (2)
#' the top `n_cell` genes with the most negative residual (`down_idx`,
#' under-expressed / dropped out relative to baseline). Fewer than `n_cell`
#' genes are returned on a side when that many signed residuals are unavailable
#' (wrong-sign genes are never padded in).
#' Count differences within 64 machine epsilons times the larger of 1, observed
#' and expected counts are numerical zeroes. [score_cells()] also allows for
#' summation over the cells contributing to each baseline group.
#'
#' This residual formulation (P2) replaces the earlier log-ratio delta and its
#' dropout-only DOWN side: the DOWN signature is now every gene expressed below
#' expectation ranked by deficit, which is robust across cell types (a
#' baseline-high housekeeping gene that merely dropped out no longer dominates
#' the DOWN list).
#'
#' @param obs Numeric vector: observed counts for one cell, one entry per gene
#'   in the scoring universe.
#' @param expected Numeric vector: expected counts under the baseline for this
#'   cell (baseline rate x library size), same length and order as `obs`.
#' @param n_cell Max number of genes per side. Default 500.
#'
#' @return A list with `up_idx` and `down_idx` integer vectors (1-based indices
#'   into the universe).
#' @keywords internal
build_cell_signature <- function(obs, expected, n_cell = 500) {
  r <- rank_cell_residuals(obs, expected)
  list(up_idx = utils::head(r$up, n_cell), down_idx = utils::head(r$down, n_cell))
}

#' Full signed ranking of one cell's genes by the P2 residual
#'
#' All genes with a positive residual, largest first (`up`), and all genes with
#' a negative residual, most negative first (`down`). [build_cell_signature()]
#' keeps the first `n_cell` of each; when a panel has its own pair universe
#' (see [score_cells()]) the lists are first restricted to that universe and
#' then cut, so the ranking must not be cut before. Added in 0.3.5.
#'
#' @inheritParams build_cell_signature
#' @param baseline_size Number of cells contributing to the baseline. The
#'   count-scale zero tolerance is `8 * .Machine$double.eps * max(8, baseline_size) *
#'   pmax(1, abs(obs), abs(expected))`.
#' @return A list with integer vectors `up` and `down` (1-based indices into
#'   the universe).
#' @keywords internal
rank_cell_residuals <- function(obs, expected, baseline_size = 1L) {
  # Same count-scale bound as Python: allow for summation in a group's mean
  # plus the CPM conversions. This is a numerical bound, not an effect-size filter.
  diff <- obs - expected
  tolerance <- 8 * .Machine$double.eps * max(8, baseline_size) * pmax(1, abs(obs), abs(expected))
  diff[abs(diff) <= tolerance] <- 0
  resid <- diff / sqrt(expected + 1)
  pos <- which(resid > 0)
  neg <- which(resid < 0)
  list(up   = if (length(pos)) pos[order(resid[pos], decreasing = TRUE)] else integer(0),
       down = if (length(neg)) neg[order(resid[neg], decreasing = FALSE)] else integer(0))
}

#' Prepare a panel for scoring against a dataset: pair universe and lists
#'
#' The background of the Running Fisher comparison is the pair universe: the
#' genes of the dataset (`gene_names`) that the reference could contain. A
#' ranked reference built from a DEG table carries `eligible`, the genes its
#' source study reported with a finite, non-zero change (see
#' [build_panel_from_deg()]); its pair universe is `gene_names` intersected with
#' `eligible`. A panel without `eligible` (a membership gene set, which can
#' contain any gene, or a panel whose tested genes are unknown) uses all of
#' `gene_names`. The size of the universe counts distinct lower-case symbols
#' (e.g. `PISD` and `Pisd` count once; before 0.3.5 every row counted). A
#' ranked list is first cut to its top `top_n` genes per side and then
#' restricted to the pair universe (order kept), so that a reference keeps the
#' genes it ranks highest whatever the dataset measures; a gene set
#' (`ranked = FALSE`) is used whole, restricted to the pair universe. This is
#' the rule of the RFind-sc paper (<https://github.com/tmurano/RFind-sc>).
#' Added in 0.3.5; since 0.3.8 the cut comes before the restriction (before,
#' the lists were restricted first and then cut, and gene sets were cut too).
#'
#' @param b A panel in internal form ([to_b2()]).
#' @param gene_names Lower-case gene symbols of the scoring universe.
#' @param top_n Cap per side (`config$top_n`).
#' @return `b` with restricted `up` / `down`, `n_up` / `n_down`,
#'   `universe_size` and `in_u` (a logical vector over `gene_names`, or `NULL`
#'   when the pair universe is the whole scoring universe).
#' @keywords internal
prepare_panel <- function(b, gene_names, top_n) {
  if (is.null(b$eligible)) {
    b$in_u <- NULL
    u <- gene_names
    b$universe_size <- length(unique(gene_names))   # symbols differing only in case count once, as in the paper
  } else {
    b$in_u <- gene_names %in% b$eligible
    u <- gene_names[b$in_u]
    b$universe_size <- length(unique(u))
  }
  cut <- function(g) {                 # ranked: top_n of the reference first, then the pair universe (0.3.8)
    if (!isFALSE(b$ranked)) g <- utils::head(g, top_n)
    g[g %in% u]
  }
  b$up <- cut(b$up); b$down <- cut(b$down)
  b$n_up <- length(b$up); b$n_down <- length(b$down)
  b
}

#' Score one cell's ranking against one prepared panel
#'
#' Cuts the cell's UP / DOWN lists to `n_cell` genes inside the panel's pair
#' universe (all genes when `b$in_u` is `NULL`) and returns the four-term score
#' ([rfind_score_cell()]). Added in 0.3.5.
#'
#' @param ranking Output of [rank_cell_residuals()].
#' @param b A panel from [prepare_panel()].
#' @param gene_names Lower-case gene symbols of the scoring universe.
#' @param n_cell Genes per side of the cell's lists.
#' @param clip Score cap.
#' @return Numeric score.
#' @keywords internal
score_cell_ranking <- function(ranking, b, gene_names, n_cell, clip = 300) {
  up <- ranking$up; dn <- ranking$down
  if (!is.null(b$in_u)) { up <- up[b$in_u[up]]; dn <- dn[b$in_u[dn]] }
  rfind_score_cell(gene_names[utils::head(up, n_cell)], gene_names[utils::head(dn, n_cell)], b,
                   universe_size = b$universe_size, clip = clip)
}

#' Score all cells against all panels (streaming, parallel)
#'
#' Computes RFind-sc per-cell scores for every panel in `biosets`. Cells are
#' processed independently and in parallel; per-cell signatures are built once
#' (see [build_cell_signature()]) and reused across panels.
#'
#' @param counts A sparse or dense gene x cell counts matrix (raw counts).
#'   Library size is taken as `colSums(counts)`, so pass the full counts (do
#'   not pre-subset genes) for a faithful baseline. Counts must be finite and
#'   non-negative; gene row names must be present, non-missing and non-empty.
#' @param biosets A bioset list (see [default_biosets()]). A panel with
#'   `ranked = FALSE` (a membership gene set, e.g. the bundled `Interferon`
#'   group) is used whole; other panels are ranked lists (see [to_b2()] and
#'   [running_fisher()]). A panel with `eligible` (the genes its source study
#'   reported, as set by [build_panel_from_deg()] and in the bundled ranked
#'   references) is scored in its pair universe (see Details).
#' @param config Configuration from [rfindsc_config()].
#' @param return_signatures Logical. If TRUE, also return per-cell UP/DOWN
#'   index lists. Default FALSE.
#' @param libsize Optional numeric vector: total UMI per cell (one per column
#'   of `counts`). Default `colSums(counts)`. Pass the totals over all genes when
#'   `counts` has already been restricted to a gene universe (as in the RFind-sc
#'   paper, where the universe is chosen before scoring). Values must be finite
#'   and strictly positive.
#' @param baseline_groups Optional vector (character, factor or numeric) with
#'   one label per cell (column of `counts`), e.g. the cell type. Default `NULL`:
#'   one dataset-wide baseline, the mean CPM over all cells. When given, each
#'   cell's expected counts use the mean CPM over the cells of its own group
#'   instead, as the RFind-sc paper does for De Donno et al. (per cell type,
#'   both conditions pooled, so that cell identity does not dominate the
#'   signature). The scoring universe (`config$nonzero_rate`) is still chosen
#'   over all cells. Missing labels are not allowed. Same argument as
#'   `baseline_groups` in the Python package rfindsc. Added in 0.3.3.
#'
#' @details
#' Expected counts: `exp[g, j] = pi[g, k(j)] * libsize[j]`, where
#' `pi[g, k] = mean over the cells j of group k of counts[g, j] / libsize[j]`
#' (the mean CPM / 1e6). Without `baseline_groups` there is a single group.
#'
#' Pair universe (0.3.5): the background of each comparison is the scoring
#' universe intersected with the panel's `eligible` genes, the genes the
#' reference could contain; a panel without `eligible` uses the whole scoring
#' universe (see [prepare_panel()]). Each cell's UP / DOWN lists are
#' restricted to it before the top `config$n_cell` genes are taken; a ranked
#' panel keeps its top `config$top_n` genes per side and then those in the pair
#' universe (0.3.8), and a gene set is used whole, as in the RFind-sc paper.
#'
#' Parallelism and memory: cells are scored in chunks of up to 200 cells with
#' `parallel::mclapply()` on `config$n_cores` cores (one core on Windows by
#' default). Before the workers start, the counts of the scoring universe are
#' copied into a dense matrix when it takes at most 2 GB (genes x cells x 8
#' bytes; 8 GB before 0.3.9); larger inputs are made dense one chunk at a time,
#' which keeps the memory small and is as fast on several cores. The per-cell
#' signatures are returned from the workers only with `return_signatures = TRUE`
#' (0.3.9). Two environment variables change these settings: `RFIND_CHUNK_SIZE`
#' (largest chunk, in cells) and `RFIND_DENSE_MAX_GB` (largest dense copy, in GB;
#' a larger value can make one-core runs on large inputs faster on a machine
#' with enough memory). Unset, or set to a value that is not a number, each
#' gives its default (200 cells, 2 GB). Scores do not depend on them.
#'
#' @return A list:
#'   * `scores`: numeric matrix (cell x panel)
#'   * `universe`: character vector of genes in the scoring universe
#'   * `universe_size`: integer vector, the pair universe size per panel
#'     (0.3.5)
#'   * `signatures`: (optional) list of per-cell `up_idx` / `down_idx` (the top
#'     `n_cell` genes of the whole scoring universe)
#'
#' @export
score_cells <- function(counts, biosets, config = rfindsc_config(),
                        return_signatures = FALSE, libsize = NULL,
                        baseline_groups = NULL) {
  stopifnot(is.list(biosets), length(biosets) > 0)
  panel_names <- names(biosets)
  if (is.null(panel_names) || anyNA(panel_names) || any(!nzchar(trimws(panel_names)))) {
    stop("`biosets` must have unique non-empty names for each panel.")
  }
  axis_names <- panel_names[panel_names != "meta"]
  if (anyDuplicated(axis_names)) stop("`biosets` must have unique non-empty names for each panel.")
  stopifnot(length(axis_names) > 0)

  if (!(is.matrix(counts) || inherits(counts, "Matrix")) || ncol(counts) == 0L) {
    stop("counts must be a gene x cell matrix with at least one cell.")
  }
  genes <- rownames(counts)
  if (is.null(genes) || anyNA(genes) || any(!nzchar(trimws(genes)))) {
    stop("gene names (rownames of counts) are required and must not be missing or empty.")
  }
  values <- if (inherits(counts, "Matrix")) {
    if ("x" %in% methods::slotNames(counts)) counts@x else 1
  } else counts
  if (!(is.numeric(values) || is.logical(values)) || is.complex(values) ||
      any(!is.finite(values)) || any(values < 0)) {
    stop("counts must contain only finite, non-negative values.")
  }
  if (is.null(libsize)) libsize <- Matrix::colSums(counts)
  if (!is.numeric(libsize) || length(libsize) != ncol(counts) ||
      any(!is.finite(libsize)) || any(libsize <= 0)) {
    stop("`libsize` must contain one finite, positive library size per cell.")
  }
  group <- baseline_group_index(baseline_groups, ncol(counts))
  cpm <- Matrix::t(Matrix::t(counts) / libsize) * 1e6
  keep <- Matrix::rowMeans(cpm > 0) >= config$nonzero_rate
  counts_u <- counts[keep, , drop = FALSE]   # raw counts, scoring universe only
  cpm <- cpm[keep, , drop = FALSE]
  universe <- rownames(cpm)
  p_universe <- length(universe)
  # Baseline expected-count rate per gene (rows) and baseline group (columns):
  # mean CPM / 1e6 over the cells of the group (one group = all cells unless
  # `baseline_groups` is given). Multiplying by a cell's library size gives that
  # cell's expected count, matching the P2 delta used to build per-cell
  # signatures.
  expc_rate <- if (is.null(group)) {
    matrix(as.numeric(Matrix::rowMeans(cpm)) / 1e6, ncol = 1L)
  } else {
    vapply(seq_len(max(group)), function(k)
      as.numeric(Matrix::rowMeans(cpm[, group == k, drop = FALSE])) / 1e6,
      numeric(nrow(cpm)))
  }
  if (is.null(group)) group <- rep(1L, ncol(counts))
  baseline_sizes <- tabulate(group)
  gene_names <- tolower(universe)
  # Names are matched without regard to case, so two rows whose names differ only in case (e.g. PISD and Pisd) count as
  # the same gene twice in a cell's lists (0.3.9: a warning, once; the scores are computed as before)
  dup_names <- unique(gene_names[duplicated(gene_names)])
  if (length(dup_names)) {
    warning(sprintf(paste("%d gene name(s) occur in more than one row of the scoring universe when case is ignored",
                          "(e.g. %s); such a gene can be counted twice in a cell's lists. Merge or rename these rows",
                          "before scoring."), length(dup_names), paste(utils::head(dup_names, 3), collapse = ", ")),
            call. = FALSE)
  }

  # Each panel is scored in its pair universe (prepare_panel): the genes of the
  # dataset that the reference could contain (`eligible`), or all of them for a
  # panel without `eligible` (a gene set). A ranked panel keeps its top
  # config$top_n genes per side and then those in the pair universe (order kept;
  # 0.3.8), as in the RFind-sc paper: the hypergeometric null is the pair
  # universe, so genes that can never be in both lists are not counted. A gene set
  # (b$ranked = FALSE) has no order: it is used whole, with K = the number of
  # its genes in the universe (rfind_score_cell -> running_fisher(g2_ranked)).
  b2_list <- lapply(axis_names, function(nm)
    prepare_panel(to_b2(biosets[[nm]]), gene_names, config$top_n))
  names(b2_list) <- axis_names
  duplicated_panels <- axis_names[vapply(b2_list, function(b)
    anyDuplicated(b$up) > 0L || anyDuplicated(b$down) > 0L, logical(1))]
  if (length(duplicated_panels)) {
    warning(structure(class = c("rfindsc_duplicated_names", "warning", "condition"),
      list(message = paste0("Reference panel(s) name a gene more than once after case-insensitive matching: ",
                            paste(duplicated_panels, collapse = ", "),
                            ". Remove duplicates before scoring; reference ranks are otherwise unchanged."),
           call = NULL)))
  }

  n_cells <- ncol(counts_u)
  cell_names <- colnames(counts_u)
  libsize <- as.numeric(libsize)

  # Pre-densify the universe counts matrix in the parent process before
  # fork(). Sparse Matrix S4 slot dereferences (`@i`, `@p`, `@x`) trigger
  # refcount writes under R 4.5.1, and those writes inside mclapply children
  # cause fork-COW page faults that effectively serialize the workers. A dense
  # matrix's refcount sits on a single header page, so column reads from shared
  # data storage do not write to other pages and workers run in parallel
  # without COW contention. Memory cost is P x N x 8 bytes (~430 MB for 6k
  # cells, ~5.7 GB for 75k cells of 9.4k genes); above the cap the workers make
  # each chunk dense themselves. Measured in 0.3.9 (R 4.5.2, Apple M3, 16 GB;
  # 75,572 cells x 9,414 genes): per-chunk extraction was as fast on 6 cores and
  # used far less memory; with the 5.7 GB copy, 6 cores had been no faster than
  # one on an 8 GB machine (0.3.8). On one core the copy was about 10% faster. Hence the cap of 2 GB (8 GB before 0.3.9), configurable
  # via env var RFIND_DENSE_MAX_GB (in GB; an unset variable or a value that is
  # not a number gives the default; 0 turns the copy off). Scores do not depend
  # on it.
  mem_gb <- nrow(counts_u) * ncol(counts_u) * 8 / 1e9
  dense_max_gb <- {
    env_v <- suppressWarnings(as.numeric(Sys.getenv("RFIND_DENSE_MAX_GB")))
    if (!is.na(env_v)) env_v else 2
  }
  predensify <- mem_gb <= dense_max_gb
  if (predensify) {
    counts_for_workers <- as.matrix(counts_u)
  } else {
    counts_for_workers <- counts_u  # workers will densify per chunk
  }

  # Chunked dispatch: process up to chunk_size cells per task to amortize
  # mclapply dispatch overhead. Adaptive: smaller of CHUNK_SIZE_MAX (memory
  # / dispatch cap) and ceiling(n_cells / n_cores) so all cores get work
  # even on small inputs. Configurable via env var RFIND_CHUNK_SIZE (in cells;
  # default 200; values below 1 count as 1). As for RFIND_DENSE_MAX_GB, an
  # unset variable or a value that cannot be read as an integer gives the
  # default.
  chunk_size_max <- {
    env_v <- suppressWarnings(as.integer(Sys.getenv("RFIND_CHUNK_SIZE")))
    if (!is.na(env_v)) max(1L, env_v) else 200L
  }
  chunk_size <- min(chunk_size_max,
                    max(1L, as.integer(ceiling(n_cells / config$n_cores))))
  chunk_starts <- seq(1L, n_cells, by = chunk_size)
  chunk_idx_list <- lapply(chunk_starts, function(s)
    seq.int(s, min(s + chunk_size - 1L, n_cells)))

  result_chunks <- parallel::mclapply(chunk_idx_list, function(idx_chunk) {
    chunk_data <- if (predensify) {
      counts_for_workers[, idx_chunk, drop = FALSE]
    } else {
      as.matrix(counts_for_workers[, idx_chunk, drop = FALSE])
    }
    lapply(seq_along(idx_chunk), function(j) {
      obs <- chunk_data[, j]
      expected <- expc_rate[, group[idx_chunk[j]]] * libsize[idx_chunk[j]]
      ranking <- rank_cell_residuals(obs, expected, baseline_sizes[group[idx_chunk[j]]])
      scores <- vapply(axis_names, function(nm) {
        tryCatch(
          # the duplicated names were reported once above, not for every cell
          withCallingHandlers(
            score_cell_ranking(ranking, b2_list[[nm]], gene_names,
                               n_cell = config$n_cell, clip = config$clip),
            rfindsc_duplicated_names = function(w) invokeRestart("muffleWarning")),
          error = function(e) NA_real_
        )
      }, numeric(1))
      # the per-cell signatures are sent back from the workers only when asked for (0.3.9): for 75,572 cells they are
      # about 300 MB, and copying them back took about a third of the time on 6 cores
      if (isTRUE(return_signatures)) {
        list(scores = scores,
             up_idx = utils::head(ranking$up, config$n_cell),
             down_idx = utils::head(ranking$down, config$n_cell))
      } else {
        list(scores = scores)
      }
    })
  }, mc.cores = config$n_cores)

  results <- unlist(result_chunks, recursive = FALSE)

  score_mat <- do.call(rbind, lapply(results, `[[`, "scores"))
  rownames(score_mat) <- cell_names
  colnames(score_mat) <- axis_names

  out <- list(scores = score_mat, universe = universe,
              universe_size = vapply(b2_list, `[[`, 0L, "universe_size"))
  if (isTRUE(return_signatures)) {
    out$signatures <- list(
      up_idx   = lapply(results, `[[`, "up_idx"),
      down_idx = lapply(results, `[[`, "down_idx"),
      n_cell   = config$n_cell,
      universe_size = p_universe
    )
  }
  out
}

#' Baseline group of each cell (helper for `score_cells(baseline_groups = )`)
#'
#' @param baseline_groups `NULL` or one label per cell.
#' @param n_cells Number of cells.
#' @return `NULL` (one dataset-wide baseline) or an integer vector of group
#'   indices `1..K`, groups in sorted label order.
#' @keywords internal
baseline_group_index <- function(baseline_groups, n_cells) {
  if (is.null(baseline_groups)) return(NULL)
  labels <- as.character(baseline_groups)
  if (length(labels) != n_cells) {
    stop("`baseline_groups` must have one label per cell (column of `counts`).")
  }
  if (anyNA(labels)) stop("`baseline_groups` must not contain missing values.")
  match(labels, sort(unique(labels)))
}
