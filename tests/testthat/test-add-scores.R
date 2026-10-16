test_that("add_rfind_scores adds one column per panel, matched by cell name", {
  sc <- matrix(c(1, 2, 3, -1, -2, -3), nrow = 3,
               dimnames = list(c("c1", "c2", "c3"), c("A", "B")))
  md <- data.frame(x = 1:4, row.names = c("c3", "c1", "c2", "c4"))
  expect_message(out <- add_rfind_scores(md, list(scores = sc)), "1 cell")
  expect_equal(out$RFind_A, c(3, 1, 2, NA))
  expect_equal(out$RFind_B, c(-3, -1, -2, NA))
  expect_equal(out$x, 1:4)                                   # other columns kept
  expect_equal(rownames(out), rownames(md))

  again <- suppressMessages(add_rfind_scores(out, sc * 10))  # overwrite (default)
  expect_equal(again$RFind_A, c(30, 10, 20, NA))
  expect_error(add_rfind_scores(out, sc, overwrite = FALSE), "already present")
  expect_equal(names(suppressMessages(add_rfind_scores(md, sc, prefix = "idx_")))[2:3],
               c("idx_A", "idx_B"))
})

test_that("add_rfind_scores rejects scores that do not name the cells", {
  md <- data.frame(x = 1:2, row.names = c("c1", "c2"))
  expect_error(add_rfind_scores(md, matrix(1:4, 2)), "cell names")
  expect_error(add_rfind_scores(md, matrix(1, 1, 1, dimnames = list("z", "A"))), "No cell")
  expect_error(add_rfind_scores(list(), matrix(1, 1, 1, dimnames = list("c1", "A"))), "Seurat object")
})

test_that("add_rfind_scores writes the scores of a Seurat object into its metadata", {
  skip_if_not_installed("SeuratObject", "5.0.0")
  skip_if_not_installed("Seurat", "5.0.0")
  set.seed(42)
  n_genes <- 200; n_cells <- 40
  counts <- Matrix::Matrix(rpois(n_genes * n_cells, lambda = 2), nrow = n_genes, sparse = TRUE,
                           dimnames = list(paste0("g", seq_len(n_genes)), paste0("c", seq_len(n_cells))))
  seu <- suppressWarnings(SeuratObject::CreateSeuratObject(counts))
  out <- score_genesets(seu, list(setA = paste0("g", 1:20), setB = paste0("g", 21:40)))
  seu <- add_rfind_scores(seu, out)
  expect_true(all(c("RFind_setA", "RFind_setB") %in% colnames(seu[[]])))
  expect_equal(seu[[]][colnames(seu), "RFind_setA"], unname(out$scores[colnames(seu), "setA"]))
  expect_equal(seu[[]][colnames(seu), "RFind_setB"], unname(out$scores[colnames(seu), "setB"]))
})
