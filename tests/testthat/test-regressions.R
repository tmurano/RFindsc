test_that("identical populations have zero scores and empty signatures", {
  values <- c(rep(5, 20), rep(3, 80))
  genes <- paste0("g", seq_len(100))
  for (nc in c(10L, 40L, 1000L)) for (sparse in c(FALSE, TRUE)) {
    counts <- matrix(values, 100, nc, dimnames = list(genes, paste0("c", seq_len(nc))))
    if (sparse) counts <- Matrix::Matrix(counts, sparse = TRUE)
    out <- score_genesets(counts, list(markers = genes[1:20]), rfindsc_config(n_cores = 1),
                          return_signatures = TRUE)
    expect_true(all(out$scores == 0))
    expect_true(all(lengths(out$signatures$up_idx) == 0L))
    expect_true(all(lengths(out$signatures$down_idx) == 0L))
  }
})

test_that("proportional populations and singleton baselines have zero scores", {
  genes <- paste0("g", seq_len(100))
  counts <- outer(rep(1:10, 10), c(1, 2, 3, 7, 10))
  dimnames(counts) <- list(genes, paste0("c", seq_len(5)))
  cfg <- rfindsc_config(n_cores = 1)
  expect_true(all(score_genesets(counts, list(A = genes[1:20]), cfg)$scores == 0))
  set.seed(21)
  counts <- matrix(rpois(100 * 40, 3), 100, dimnames = list(genes, paste0("c", seq_len(40))))
  out <- score_genesets(counts, list(A = genes[1:20]), cfg,
                        baseline_groups = seq_len(40), return_signatures = TRUE)
  expect_true(all(out$scores == 0))
  expect_true(all(lengths(out$signatures$up_idx) == 0L))
  expect_true(all(lengths(out$signatures$down_idx) == 0L))
})

test_that("the numerical zero guard preserves resolvable signed differences", {
  sig <- build_cell_signature(c(5 + 1e-14, 5 + 1e-8, 5 - 1e-8), c(5, 5, 5))
  expect_identical(sig$up_idx, 2L)
  expect_identical(sig$down_idx, 3L)
})

test_that("invalid counts are rejected even with external library sizes", {
  for (bad in c(NA_real_, Inf, -Inf, -1)) for (sparse in c(FALSE, TRUE)) {
    counts <- matrix(c(bad, 1, 1, 2), 2, dimnames = list(c("a", "b"), c("c1", "c2")))
    if (sparse) counts <- Matrix::Matrix(counts, sparse = TRUE)
    expect_error(score_genesets(counts, list(A = "a"), rfindsc_config(n_cores = 1), libsize = c(10, 10)),
                 "finite, non-negative")
  }
})

test_that("library sizes must be finite and positive", {
  counts <- matrix(1, 2, 2, dimnames = list(c("a", "b"), c("c1", "c2")))
  for (bad in c(NA_real_, Inf, 0, -1)) {
    expect_error(score_genesets(counts, list(A = "a"), rfindsc_config(n_cores = 1), libsize = c(bad, 10)),
                 "finite, positive")
  }
})

test_that("gene names are required before matching references", {
  counts <- matrix(1, 2, 2, dimnames = list(NULL, c("c1", "c2")))
  expect_error(score_genesets(counts, list(A = "a")), "gene names")
  for (bad in c(NA_character_, "", " ")) {
    rownames(counts) <- c("a", bad)
    expect_error(score_genesets(counts, list(A = "a")), "gene names")
  }
})

test_that("duplicate panel names are rejected and metadata may be combined", {
  counts <- matrix(c(10, 1, 1, 10), 2, dimnames = list(c("a", "b"), c("c1", "c2")))
  cfg <- rfindsc_config(n_cores = 1)
  expect_error(score_genesets(counts, list(A = "a", A = "b"), cfg), "unique non-empty names")
  expect_error(score_panel_pairs(counts, list(A = list(up = "a"), A = list(up = "b")), cfg),
               "unique non-empty names")
  panel <- list(up = data.frame(gene = "a"), down = data.frame(gene = character()))
  expect_error(score_cells(counts, list(A = panel, A = panel), cfg), "unique non-empty names")
  out <- score_cells(counts, list(A = panel, meta = list(), meta = list()), cfg)
  expect_identical(colnames(out$scores), "A")
})

test_that("reference duplicates warn once before workers and retain reference ranks", {
  counts <- matrix(c(10, 1, 5, 2, 1, 10, 5, 2), 4,
                   dimnames = list(c("a", "b", "c", "d"), c("c1", "c2")))
  for (cores in if (.Platform$OS.type == "windows") 1L else c(1L, 2L)) {
    warnings <- capture_warnings(out <- score_panel_pairs(counts,
      list(A = list(up = c("a", " A ", "b"))), rfindsc_config(n_cores = cores)))
    expect_length(grep("Reference panel", warnings), 1L)
    expect_equal(unname(out$scores[, "A"]), c(log10(3), -log10(3)), tolerance = 1e-12)
  }
})

test_that("default cores fall back to one when detection is unavailable", {
  local_mocked_bindings(detectCores = function(...) NA_integer_, .package = "parallel")
  expect_identical(rfindsc_config()$n_cores, 1L)
  counts <- matrix(c(10, 1, 1, 10), 2, dimnames = list(c("a", "b"), c("c1", "c2")))
  expect_equal(score_genesets(counts, list(A = "a"))$scores,
                score_genesets(counts, list(A = "a"), rfindsc_config(n_cores = 1))$scores)
})
