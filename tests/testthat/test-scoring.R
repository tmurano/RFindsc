test_that("build_cell_signature ranks by P2 residual (over/under expectation)", {
  # obs above expectation -> UP; obs below expectation (incl. dropout) -> DOWN.
  obs      <- c(g1 = 50, g2 = 0,  g3 = 5, g4 = 20, g5 = 0)
  expected <- c(g1 = 10, g2 = 8,  g3 = 5, g4 = 40, g5 = 1)
  resid <- (obs - expected) / sqrt(expected + 1)      # +,-,0,-,-
  sig <- build_cell_signature(obs, expected, n_cell = 500)
  # Only strictly-positive residuals in UP; strictly-negative in DOWN; 0 in neither
  expect_equal(unname(sig$up_idx), unname(which(resid > 0)))
  # DOWN ordered most-negative first
  expect_equal(unname(sig$down_idx),
               unname(which(resid < 0)[order(resid[resid < 0])]))
  expect_false(3L %in% c(sig$up_idx, sig$down_idx))   # g3 residual == 0
})

test_that("build_cell_signature caps each side at n_cell", {
  set.seed(7)
  obs      <- c(rpois(300, 20), rep(0, 200))          # 300 likely up, 200 dropouts
  expected <- rep(5, 500)
  sig <- build_cell_signature(obs, expected, n_cell = 50)
  expect_lte(length(sig$up_idx), 50)
  expect_lte(length(sig$down_idx), 50)
})

test_that("score_cells returns a cell x panel matrix on a toy dataset", {
  set.seed(42)
  n_genes <- 200; n_cells <- 40
  counts <- matrix(rpois(n_genes * n_cells, lambda = 2),
                   nrow = n_genes, ncol = n_cells)
  rownames(counts) <- paste0("g", seq_len(n_genes))
  colnames(counts) <- paste0("c", seq_len(n_cells))

  # Minimal 2-panel bioset
  bs <- list(
    P1 = list(up = data.frame(gene = paste0("g", 1:20),
                              diff = runif(20), rank = 1:20),
              down = data.frame(gene = paste0("g", 21:40),
                                diff = -runif(20), rank = 1:20),
              type = "deg"),
    P2 = list(up = data.frame(gene = paste0("g", 41:60),
                              diff = runif(20), rank = 1:20),
              down = data.frame(gene = character(0),
                                diff = numeric(0),
                                rank = integer(0)),
              type = "gene_set")
  )
  cfg <- rfindsc_config(n_cell = 50, n_cores = 1, nonzero_rate = 0)
  out <- score_cells(counts, bs, config = cfg)
  expect_true(is.matrix(out$scores))
  expect_equal(dim(out$scores), c(n_cells, 2))
  expect_equal(colnames(out$scores), c("P1", "P2"))
})

test_that("score_genesets accepts AUCell-style gene set list", {
  set.seed(43)
  n_genes <- 200; n_cells <- 30
  counts <- matrix(rpois(n_genes * n_cells, lambda = 2),
                   nrow = n_genes, ncol = n_cells)
  rownames(counts) <- paste0("g", seq_len(n_genes))
  colnames(counts) <- paste0("c", seq_len(n_cells))

  gene_sets <- list(
    set_a = paste0("g", 1:20),
    set_b = paste0("g", 50:75)
  )
  cfg <- rfindsc_config(n_cell = 50, n_cores = 1, nonzero_rate = 0)
  out <- score_genesets(counts, gene_sets, config = cfg)
  expect_true(is.matrix(out$scores))
  expect_equal(dim(out$scores), c(n_cells, 2))
  expect_equal(colnames(out$scores), c("set_a", "set_b"))
})

test_that("score_panel_pairs accepts character vector UP/DOWN", {
  set.seed(44)
  n_genes <- 200; n_cells <- 30
  counts <- matrix(rpois(n_genes * n_cells, lambda = 2),
                   nrow = n_genes, ncol = n_cells)
  rownames(counts) <- paste0("g", seq_len(n_genes))
  colnames(counts) <- paste0("c", seq_len(n_cells))

  panels <- list(
    pair1 = list(up = paste0("g", 1:20), down = paste0("g", 21:40)),
    pair2 = list(up = paste0("g", 41:60))   # down empty
  )
  cfg <- rfindsc_config(n_cell = 50, n_cores = 1, nonzero_rate = 0)
  out <- score_panel_pairs(counts, panels, config = cfg)
  expect_true(is.matrix(out$scores))
  expect_equal(dim(out$scores), c(n_cells, 2))
  expect_equal(colnames(out$scores), c("pair1", "pair2"))
})

test_that("score_genesets validates input", {
  set.seed(45)
  counts <- matrix(rpois(20 * 5, 1), nrow = 20)
  rownames(counts) <- paste0("g", seq_len(20))
  colnames(counts) <- paste0("c", seq_len(5))
  cfg <- rfindsc_config(n_cell = 5, n_cores = 1, nonzero_rate = 0)
  expect_error(score_genesets(counts, list(), config = cfg),
               "non-empty named list")
  expect_error(score_genesets(counts,
                               list(c("g1", "g2")), config = cfg),
               "unique non-empty names")
})

test_that("score_genesets uses each set whole: gene order does not matter", {
  set.seed(46)
  n_genes <- 200; n_cells <- 30
  counts <- matrix(rpois(n_genes * n_cells, lambda = 2),
                   nrow = n_genes, ncol = n_cells)
  counts[1:20, 1:10] <- counts[1:20, 1:10] + rpois(200, 6)
  rownames(counts) <- paste0("g", seq_len(n_genes))
  colnames(counts) <- paste0("c", seq_len(n_cells))
  cfg <- rfindsc_config(n_cell = 50, n_cores = 1, nonzero_rate = 0)
  sets <- list(set_a = paste0("g", c(1:20, 101:110)), set_b = paste0("g", 50:75))
  out <- score_genesets(counts, sets, config = cfg)$scores
  shuffled <- lapply(sets, function(g) rev(g))
  expect_identical(score_genesets(counts, shuffled, config = cfg)$scores, out)
  # the same sets as UP-only panels marked ranked = FALSE give the same scores
  pairs <- lapply(sets, function(g) list(up = sort(g), ranked = FALSE))
  expect_identical(score_panel_pairs(counts, pairs, config = cfg)$scores, out)
  # as ranked UP-only panels the order matters (the scores are those of the ranked kernel)
  ranked <- score_panel_pairs(counts, lapply(sets, function(g) list(up = g)), config = cfg)$scores
  expect_false(identical(ranked, out))
})

test_that("a panel with `eligible` is scored in its pair universe (0.3.5)", {
  set.seed(45)
  n_genes <- 300; n_cells <- 20
  counts <- matrix(rpois(n_genes * n_cells, lambda = 2), nrow = n_genes, ncol = n_cells,
                   dimnames = list(paste0("g", seq_len(n_genes)), paste0("c", seq_len(n_cells))))
  up <- paste0("g", 1:40); dn <- paste0("g", 41:80); el <- paste0("g", c(1:80, 151:250))
  bs <- list(
    A = list(up = data.frame(gene = up), down = data.frame(gene = dn), type = "deg", eligible = el),
    B = list(up = data.frame(gene = up), down = data.frame(gene = dn), type = "deg")
  )
  cfg <- rfindsc_config(n_cell = 30, n_cores = 1, nonzero_rate = 0)
  out <- score_cells(counts, bs, config = cfg)
  expect_equal(unname(out$universe_size), c(180L, 300L))
  # by hand: restrict each cell's ranking to the pair universe, then cut to n_cell
  gl <- tolower(rownames(counts)); ls <- colSums(counts)
  rate <- rowMeans(t(t(counts) / ls) * 1e6) / 1e6
  inu <- gl %in% tolower(el)
  b <- to_b2(bs$A)
  by_hand <- vapply(seq_len(n_cells), function(j) {
    r <- rank_cell_residuals(counts[, j], rate * ls[j])
    rfind_score_cell(gl[utils::head(r$up[inu[r$up]], 30)], gl[utils::head(r$down[inu[r$down]], 30)], b,
                     universe_size = 180)
  }, numeric(1))
  expect_equal(unname(out$scores[, "A"]), by_hand)
  # `eligible` covering every gene gives the same scores as no `eligible`
  bs$C <- bs$B; bs$C$eligible <- rownames(counts)
  out2 <- score_cells(counts, bs, config = cfg)
  expect_equal(out2$scores[, "C"], out2$scores[, "B"])
})

test_that("build_panel_from_deg records the measured genes as `eligible` (0.3.5)", {
  deg <- data.frame(gene = c("a", "b", "c", "d", "e"), log2FC = c(2, -1, 0, NA, 0.5),
                    padj = c(0.01, 0.2, 0.5, 0.1, 0.9))
  p <- build_panel_from_deg(deg, padj_col = "padj", padj_thresh = 0.05)
  expect_setequal(p$eligible, c("a", "b", "e"))   # finite, non-zero log2FC, before the padj filter
  expect_equal(p$up$gene, "a")
})

test_that("a ranked panel is cut to top_n before the pair universe; a gene set is used whole (0.3.8)", {
  gn <- paste0("g", 1:10)
  b <- list(up = c("g11", "g1", "g2", "g3"), down = c("g4", "g12", "g5"), ranked = TRUE, eligible = NULL)
  p <- prepare_panel(b, gn, top_n = 2)
  expect_identical(p$up, "g1")                      # the top 2 are g11 and g1; g11 is not in the dataset
  expect_identical(p$down, "g4")                    # the top 2 are g4 and g12
  s <- prepare_panel(utils::modifyList(b, list(ranked = FALSE)), gn, top_n = 2)
  expect_identical(s$up, c("g1", "g2", "g3"))       # a gene set is used whole
})

test_that("score_cells warns once about names that differ only in case and returns signatures only when asked (0.3.9)", {
  set.seed(3)
  counts <- matrix(rpois(60 * 30, lambda = 3), nrow = 60, ncol = 30,
                   dimnames = list(c("Gdup", "GDUP", paste0("g", 3:60)), paste0("c", 1:30)))
  bs <- list(P = list(up = data.frame(gene = paste0("g", 3:12), diff = 10:1, rank = 1:10),
                      down = data.frame(gene = paste0("g", 13:22), diff = -(10:1), rank = 1:10), type = "deg"))
  cfg <- rfindsc_config(n_cores = 1)
  w <- capture_warnings(out <- score_cells(counts, bs, config = cfg))
  expect_length(grep("occur in more than one row", w), 1L)
  expect_null(out$signatures)
  out2 <- suppressWarnings(score_cells(counts, bs, config = cfg, return_signatures = TRUE))
  expect_equal(out2$scores, out$scores)
  expect_length(out2$signatures$up_idx, ncol(counts))
})
