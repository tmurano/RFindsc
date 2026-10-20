## Per-group (e.g. per-cell-type) baseline: score_cells(baseline_groups = ).

sim_counts <- function(seed, n_genes = 300, sizes = c(A = 25, B = 20, C = 15)) {
  set.seed(seed)
  groups <- rep(names(sizes), sizes)
  rate <- rgamma(n_genes, shape = 0.6) *
    matrix(exp(rnorm(n_genes * length(sizes), 0, 0.8)), n_genes)[, match(groups, names(sizes))]
  lib <- round(exp(rnorm(length(groups), log(3000), 0.4)))
  mu <- sweep(rate, 2, colSums(rate), "/") * rep(lib, each = n_genes)
  counts <- matrix(rnbinom(length(mu), size = 2, mu = mu), n_genes)
  dimnames(counts) <- list(paste0("g", seq_len(n_genes)), paste0("c", seq_along(groups)))
  list(counts = Matrix::Matrix(counts, sparse = TRUE), groups = groups)
}

sim_panels <- function() list(
  P1 = list(up = data.frame(gene = paste0("g", 1:40)), down = data.frame(gene = paste0("g", 41:80)), type = "deg"),
  P2 = list(up = data.frame(gene = paste0("g", seq(300, 200, by = -2))),
            down = data.frame(gene = character(0)), type = "gene_set"))

## The per-cell-type baseline and the per-cell signature loop of analysis/dedonno/02_signature.R of the RFind-sc
## paper (https://github.com/tmurano/RFind-sc), with the inputs (cf = universe counts, cs = total UMI per cell,
## ct = cell type) and NCELL as arguments. Differences from that script: this copy computes in CPM (cf x 1e6 / cs;
## expected = baseline CPM x cs / 1e6) where the script uses shares (cf / cs; expected = baseline share x cs), the
## same values up to floating-point rounding; it does not require 20 cells per cell type (one simulated group has 15
## cells); and it keeps only the top-NCELL lists, not the deeper NDEEP ranking that the script also writes.
paper_signature <- function(cf, cs, ct, NCELL = 500L, CH = 4000L) {
  G <- nrow(cf); nC <- ncol(cf)
  cpm <- cf %*% Matrix::Diagonal(x = 1e6 / cs)
  cts <- sort(unique(ct)); K <- length(cts); ctidx <- match(ct, cts)
  baseline_ct <- matrix(0, G, K, dimnames = list(NULL, cts))
  for (ti in seq_len(K)) { cols <- which(ctidx == ti)
    baseline_ct[, ti] <- Matrix::rowMeans(cpm[, cols, drop = FALSE]) }
  upI <- matrix(NA_integer_, NCELL, nC); dnI <- matrix(NA_integer_, NCELL, nC)
  for (st in seq(1L, nC, by = CH)) {
    en <- min(st + CH - 1L, nC); idx <- st:en
    obs  <- as.matrix(cf[, idx])
    expc <- sweep(baseline_ct[, ctidx[idx], drop = FALSE], 2, cs[idx] / 1e6, "*")
    P2   <- (obs - expc) / sqrt(expc + 1)
    for (k in seq_along(idx)) { p <- P2[, k]
      pos <- which(p > 0); pos <- pos[order(p[pos], decreasing = TRUE)]
      neg <- which(p < 0); neg <- neg[order(p[neg], decreasing = FALSE)]
      nu <- min(NCELL, length(pos)); nd <- min(NCELL, length(neg))
      if (nu) upI[seq_len(nu), idx[k]] <- pos[seq_len(nu)]
      if (nd) dnI[seq_len(nd), idx[k]] <- neg[seq_len(nd)] }
  }
  list(up = upI, down = dnI)
}

test_that("baseline_groups = NULL is the dataset-wide baseline (unchanged) and one group equals it", {
  d <- sim_counts(1)
  cfg <- rfindsc_config(n_cell = 50, n_cores = 1, nonzero_rate = 0)
  a <- score_cells(d$counts, sim_panels(), cfg, return_signatures = TRUE)
  b <- score_cells(d$counts, sim_panels(), cfg, return_signatures = TRUE,
                   baseline_groups = rep("all", ncol(d$counts)))
  expect_identical(a, b)
  # the dataset-wide baseline written out: rate = mean CPM / 1e6 over all cells
  lib <- Matrix::colSums(d$counts)
  rate <- Matrix::rowMeans(Matrix::t(Matrix::t(d$counts) / lib) * 1e6) / 1e6
  for (j in c(1, 30, 60)) {
    sig <- build_cell_signature(as.numeric(d$counts[, j]), rate * lib[j], n_cell = 50)
    expect_identical(unname(a$signatures$up_idx[[j]]), unname(sig$up_idx))
    expect_identical(unname(a$signatures$down_idx[[j]]), unname(sig$down_idx))
  }
})

test_that("per-cell-type signatures equal the paper pipeline (analysis/dedonno/02_signature.R)", {
  d <- sim_counts(2)
  cs <- Matrix::colSums(d$counts) + 500          # total UMI over all genes (counts = the universe only)
  for (n_cell in c(40L, 500L)) {                 # with and without the cap (300 genes)
    ref <- paper_signature(d$counts, cs, d$groups, NCELL = n_cell)
    out <- score_cells(d$counts, sim_panels(), rfindsc_config(n_cell = n_cell, n_cores = 1, nonzero_rate = 0),
                       return_signatures = TRUE, libsize = cs, baseline_groups = d$groups)
    expect_identical(out$universe, rownames(d$counts))
    for (j in seq_len(ncol(d$counts))) {
      expect_identical(unname(out$signatures$up_idx[[j]]), ref$up[!is.na(ref$up[, j]), j])
      expect_identical(unname(out$signatures$down_idx[[j]]), ref$down[!is.na(ref$down[, j]), j])
    }
  }
  # and the grouped baseline changes the signatures (it is not silently ignored)
  flat <- score_cells(d$counts, sim_panels(), rfindsc_config(n_cell = 40, n_cores = 1, nonzero_rate = 0),
                      return_signatures = TRUE, libsize = cs)
  expect_false(identical(flat$signatures$up_idx, out$signatures$up_idx))
})

test_that("baseline_groups: labels, not their order or type, define the groups", {
  d <- sim_counts(3)
  cfg <- rfindsc_config(n_cell = 50, n_cores = 1, nonzero_rate = 0.05)
  a <- score_cells(d$counts, sim_panels(), cfg, baseline_groups = d$groups)
  relabel <- c(A = "z", B = "a", C = "m")[d$groups]
  expect_identical(score_cells(d$counts, sim_panels(), cfg, baseline_groups = unname(relabel)), a)
  expect_identical(score_cells(d$counts, sim_panels(), cfg, baseline_groups = factor(d$groups, c("C", "A", "B"))), a)
  expect_identical(score_cells(d$counts, sim_panels(), cfg, baseline_groups = match(d$groups, c("C", "A", "B"))), a)
  # the universe is still chosen over all cells
  expect_identical(a$universe, score_cells(d$counts, sim_panels(), cfg)$universe)
})

test_that("baseline_groups is passed through by run_pipeline(), score_panel_pairs() and score_genesets()", {
  d <- sim_counts(4)
  cfg <- rfindsc_config(n_cell = 50, n_cores = 1, nonzero_rate = 0)
  ref <- score_cells(d$counts, sim_panels(), cfg, baseline_groups = d$groups)$scores
  expect_identical(run_pipeline(d$counts, sim_panels(), cfg, baseline_groups = d$groups)$scores, ref)
  pp <- score_panel_pairs(d$counts, list(P1 = list(up = paste0("g", 1:40), down = paste0("g", 41:80))), cfg,
                          baseline_groups = d$groups)$scores
  expect_equal(pp[, "P1"], ref[, "P1"], tolerance = 0)
  gs <- score_genesets(d$counts, list(P2 = paste0("g", seq(300, 200, by = -2))), cfg, baseline_groups = d$groups)$scores
  sets <- sim_panels()["P2"]; sets$P2$ranked <- FALSE   # score_genesets() uses a set whole (0.3.4)
  gs_ref <- score_cells(d$counts, sets, cfg, baseline_groups = d$groups)$scores
  expect_equal(gs[, "P2"], gs_ref[, "P2"], tolerance = 0)
})

test_that("baseline_groups is validated", {
  d <- sim_counts(5, sizes = c(A = 5, B = 5))
  cfg <- rfindsc_config(n_cell = 20, n_cores = 1, nonzero_rate = 0)
  expect_error(score_cells(d$counts, sim_panels(), cfg, baseline_groups = d$groups[-1]), "one label per cell")
  expect_error(score_cells(d$counts, sim_panels(), cfg, baseline_groups = replace(d$groups, 3, NA)),
               "missing values")
})

test_that("per-group baseline scores equal the Python package rfindsc (shared fixture)", {
  fx <- function(...) test_path("fixtures", "py_grouped", ...)
  skip_if_not(file.exists(fx("scores.tsv")))
  con <- gzfile(fx("counts.mtx.gz")); counts <- methods::as(Matrix::readMM(con), "CsparseMatrix")
  genes <- readLines(fx("genes.txt"))
  cells <- utils::read.delim(fx("cells.tsv"), colClasses = "character")
  dimnames(counts) <- list(genes, cells$cell)
  pan <- utils::read.delim(fx("panels.tsv"), colClasses = "character")
  rk <- utils::read.delim(fx("panels_ranked.tsv"), colClasses = c("character", "logical"))
  panels <- lapply(split(pan, factor(pan$panel, unique(pan$panel))), function(p) {
    b <- list(up = data.frame(gene = p$gene[p$side == "up"]), down = data.frame(gene = p$gene[p$side == "down"]),
              ranked = rk$ranked[rk$panel == p$panel[1]])
    if (any(p$side == "eligible")) b$eligible <- p$gene[p$side == "eligible"]   # pair universe (0.3.5)
    b })
  py <- as.matrix(utils::read.delim(fx("scores.tsv"), row.names = 1, check.names = FALSE))
  # the fixture has one name that differs from another only in case (cspp1): score_cells() says so once (0.3.9)
  expect_warning(out <- score_cells(counts, panels, rfindsc_config(n_cell = 100, n_cores = 1), return_signatures = TRUE,
                                    baseline_groups = cells$group), "occur in more than one row")
  expect_identical(colnames(out$scores), colnames(py))
  expect_identical(rownames(out$scores), rownames(py))
  expect_lt(max(abs(out$scores - py)), 1e-8)
  expect_identical(out$universe, readLines(fx("universe.txt")))
  sig <- utils::read.delim(gzfile(fx("signatures.tsv.gz")), colClasses = "character")
  idx <- lapply(strsplit(sig$idx, ","), as.integer)
  expect_identical(lapply(out$signatures$up_idx, unname), unname(idx[sig$side == "up_idx"]))
  expect_identical(lapply(out$signatures$down_idx, unname), unname(idx[sig$side == "down_idx"]))
})
