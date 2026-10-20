.default_n_cores <- function() {
  if (.Platform$OS.type == "windows") return(1L)
  n <- parallel::detectCores()
  if (length(n) != 1L || !is.finite(n)) return(1L)
  max(1L, n - 1L)
}

#' RFind-sc configuration defaults
#'
#' Default hyperparameters for RFind-sc scoring. Pass a modified copy to
#' [run_pipeline()] or [score_cells()] to override.
#'
#' @param top_n Cap on panel UP/DOWN gene list length. Default 500 (the
#'   RFind-sc paper uses the 500 most up- and 500 most down-regulated genes of
#'   a signed reference; 1000 before 0.3.0). The top `top_n` of a ranked
#'   reference are taken before its genes are restricted to those of the
#'   dataset (0.3.8); gene sets are used whole.
#' @param n_cell Cap on per-cell UP/DOWN signature length. Default 500.
#' @param nonzero_rate Minimum fraction of cells in which a gene must be
#'   non-zero to enter the scoring universe. Default 0.01.
#' @param clip Score magnitude cap (symmetric, applied per panel). Default 300.
#' @param n_cores Number of cores for `parallel::mclapply`.
#'   Default `max(1, parallel::detectCores() - 1)`, and 1 on Windows, where
#'   `mclapply` cannot fork and supports only one core. Scores do not depend
#'   on it. If core detection is unavailable, defaults to 1.
#'
#' @return A named list of settings.
#' @export
rfindsc_config <- function(top_n = 500,
                           n_cell = 500,
                           nonzero_rate = 0.01,
                           clip = 300,
                           n_cores = .default_n_cores()) {
  list(top_n = top_n,
       n_cell = n_cell,
       nonzero_rate = nonzero_rate,
       clip = clip,
       n_cores = n_cores)
}
