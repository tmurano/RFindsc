#' Running Fisher score of two ranked gene lists
#'
#' Computes `-log10` of the smallest hypergeometric upper-tail p-value,
#' `P(X >= M)`, over every pair of cut-offs: the top `G` genes of `g1` and the
#' top `K` genes of `g2` (`M` = genes shared by the two prefixes). This is the
#' kernel of RFind-sc (Kupershmidt et al. 2010 Running Fisher, G x K form). With
#' `M` fixed, `p` only grows with `G` and `K`, so the minimum lies at the pairs
#' of ranks of shared genes; only those corners are evaluated, which is
#' identical to the full G x K scan. `p` is computed on the log scale.
#'
#' A membership gene set (e.g. an MSigDB Hallmark set) has no order. With
#' `g2_ranked = FALSE` it is used whole, as a biogroup in the original Running
#' Fisher: `K` is fixed at the set size (`length(unique(g2))`) and only the
#' top-`G` cut-offs of the ranked list `g1` are scanned (at the ranks of the
#' shared genes), so the score does not depend on the order in which the set
#' is stored.
#'
#' Changed in 0.3.0: earlier versions scanned one list against the whole other
#' list in both directions and averaged the two p-values (geometric mean). The
#' current definition matches the RFind-sc paper. Changed in 0.3.4: `g2_ranked`
#' added; before, gene sets were scanned as if ranked in the order stored
#' (alphabetical for the bundled gene groups).
#'
#' Each list should name each gene once: a name given twice counts as two
#' genes and inflates the overlap `M`. Since 0.3.9 such input gives a warning
#' (class `rfindsc_duplicated_names`); the score is computed as before.
#'
#' @param g1 Character vector: ranked gene symbols (first list, most extreme
#'   first).
#' @param g2 Character vector: ranked gene symbols (second list), or a
#'   membership gene set when `g2_ranked = FALSE`.
#' @param universe_size Integer: total gene universe `P` used for the
#'   hypergeometric null.
#' @param clip Numeric: score cap (default 300).
#' @param g2_ranked Logical: `FALSE` when `g2` is an unordered gene set; then
#'   `K` is fixed at the set size and only the cut-offs of `g1` are scanned.
#'   Default `TRUE`. Added in 0.3.4.
#'
#' @return Numeric: \code{-log10(min p)}, clipped to `[0, clip]`. Returns 0 when
#'   either list is empty or they share no gene.
#'
#' @examples
#' set.seed(1)
#' pool <- paste0("g", 1:500)
#' g1 <- sample(pool, 50)
#' g2 <- c(g1[1:10], sample(setdiff(pool, g1), 40))
#' running_fisher(g1, g2, universe_size = 500)
#' # g2 as an unordered gene set: the result does not depend on its order
#' running_fisher(g1, g2, universe_size = 500, g2_ranked = FALSE)
#' running_fisher(g1, sort(g2), universe_size = 500, g2_ranked = FALSE)
#'
#' @export
running_fisher <- function(g1, g2, universe_size, clip = 300, g2_ranked = TRUE) {
  if (!length(g1) || !length(g2)) return(0)
  if (anyDuplicated(g1) || anyDuplicated(g2)) {
    warning(structure(class = c("rfindsc_duplicated_names", "warning", "condition"),
      list(message = paste("running_fisher(): a gene list names a gene more than once; each copy counts as a",
                           "separate gene and inflates the overlap. Remove the duplicates before matching."),
           call = NULL)))
  }
  r1 <- which(g1 %in% g2); m <- length(r1)                 # ranks in g1 of the shared genes (increasing)
  if (!m) return(0)
  if (!g2_ranked) {
    # Unordered gene set: K = set size; the smallest p for a given G = r1[j] has M = j.
    K <- length(unique(g2))
    lp <- stats::phyper(seq_len(m) - 1, K, universe_size - K, r1,
                        lower.tail = FALSE, log.p = TRUE)
  } else {
    r2 <- match(g1[r1], g2); K <- sort(r2)
    M <- crossprod(outer(r1, r1, "<="), outer(r2, K, "<="))   # overlap of top r1[i] of g1 and top K[j] of g2
    Kv <- rep(K, each = m)
    lp <- stats::phyper(as.vector(M) - 1, Kv, universe_size - Kv, rep(r1, times = m),
                        lower.tail = FALSE, log.p = TRUE)
  }
  s <- -min(lp) / log(10)
  if (!is.finite(s)) return(clip)
  min(s, clip)
}

#' Running Fisher score (former name)
#'
#' Kept for backward compatibility; identical to [running_fisher()] since 0.3.0
#' (the earlier two-scan geometric-mean definition was replaced).
#'
#' @inheritParams running_fisher
#' @return See [running_fisher()].
#' @export
bid_running_fisher <- function(g1, g2, universe_size, clip = 300, g2_ranked = TRUE) {
  running_fisher(g1, g2, universe_size, clip, g2_ranked = g2_ranked)
}

#' Four-term RFind-sc score for a single cell against a single panel
#'
#' Composes the signed RFind-sc score:
#' \deqn{S = RF(b1u, b2u) - RF(b1u, b2d) - RF(b1d, b2u) + RF(b1d, b2d)}
#' where `b1u/b1d` are the per-cell UP/DOWN signatures, `b2u/b2d` the panel
#' UP/DOWN gene lists and `RF` is [running_fisher()]. For a gene set without a
#' DOWN side the two terms with `b2d` are 0. Positive = cell aligned with the
#' panel's case direction, negative = aligned with control. When
#' `b2$ranked` is `FALSE` (a membership gene set) the panel lists are used
#' whole in all four terms (`g2_ranked = FALSE`).
#'
#' @param b1_up,b1_down Character vectors: per-cell UP/DOWN signatures
#'   (typically `n_cell` genes each, from [build_cell_signature()]).
#' @param b2 Panel object in internal form (list with `up`, `down`, `n_up`,
#'   `n_down`, and `ranked`: `FALSE` for a membership gene set; a missing
#'   `ranked` means ranked lists). Use [to_b2()] to convert a bioset entry.
#' @param universe_size Integer: total gene universe `P`.
#' @param clip Score magnitude cap (symmetric). Default 300.
#'
#' @return Numeric: signed RFind-sc score, clipped to `[-clip, clip]`.
#' @keywords internal
rfind_score_cell <- function(b1_up, b1_down, b2, universe_size, clip = 300) {
  rk <- !isFALSE(b2$ranked)                                 # b2$ranked = FALSE for a membership gene set
  s <- running_fisher(b1_up,   b2$up,   universe_size, clip = clip, g2_ranked = rk) -
       running_fisher(b1_up,   b2$down, universe_size, clip = clip, g2_ranked = rk) -
       running_fisher(b1_down, b2$up,   universe_size, clip = clip, g2_ranked = rk) +
       running_fisher(b1_down, b2$down, universe_size, clip = clip, g2_ranked = rk)
  max(-clip, min(s, clip))
}

#' Is a panel a ranked reference or an unordered gene set?
#'
#' @param panel A bioset panel entry.
#' @return `FALSE` when `panel$ranked` is `FALSE`, or when `panel` has no
#'   `ranked` field and its `type` is a membership gene set type:
#'   `"gene_group"` (bundled MSigDB Hallmark groups), `"msigdb_type_b"`
#'   ([msigdbr_to_panels()]) or `"msigdb_cgp_paired"` ([c2_cgp_auto_pair()]).
#'   `TRUE` otherwise.
#' @keywords internal
panel_is_ranked <- function(panel) {
  if (!is.null(panel$ranked)) return(!isFALSE(panel$ranked))
  unranked_types <- c("gene_group", "msigdb_type_b", "msigdb_cgp_paired")
  !(is.character(panel$type) && length(panel$type) == 1L && panel$type %in% unranked_types)
}

#' Convert a bioset panel entry into internal scoring form
#'
#' @param panel A list with `$up` and `$down` data.frames each containing a
#'   `gene` column (as in [default_biosets()]). An optional logical `ranked`
#'   marks a membership gene set (`ranked = FALSE`), whose lists have no order
#'   and are used whole; without it, panels of type `"gene_group"`,
#'   `"msigdb_type_b"` or `"msigdb_cgp_paired"` are gene sets and all other
#'   panels are ranked lists.
#'
#'   An optional character vector `eligible` (the genes the reference's
#'   source study reported; 0.3.5) sets the pair universe (see
#'   [prepare_panel()]).
#'
#' @return A list with lowercase gene vectors, sizes, `ranked` and `eligible`
#'   (lower case, or `NULL`), ready for [prepare_panel()] and
#'   [rfind_score_cell()].
#' @keywords internal
to_b2 <- function(panel) {
  up_sym <- tolower(trimws(panel$up$gene))
  dn_sym <- tolower(trimws(panel$down$gene))
  el <- panel$eligible
  if (!is.null(el)) { el <- unique(tolower(trimws(as.character(el)))); el <- el[!is.na(el) & nzchar(el)] }
  list(up = up_sym, down = dn_sym,
       n_up = length(up_sym), n_down = length(dn_sym),
       ranked = panel_is_ranked(panel), eligible = el)
}
