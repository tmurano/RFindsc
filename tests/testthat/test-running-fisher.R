test_that("bid_running_fisher returns 0 on empty input", {
  expect_equal(bid_running_fisher(character(0), letters[1:3], 100), 0)
  expect_equal(bid_running_fisher(letters[1:3], character(0), 100), 0)
})

test_that("bid_running_fisher is bounded by clip", {
  g <- paste0("g", 1:10)
  # Identical lists -> maximally significant; must not exceed clip
  s <- bid_running_fisher(g, g, universe_size = 1000, clip = 50)
  expect_gt(s, 0)
  expect_lte(s, 50)
  # Raising clip should not decrease the score
  s2 <- bid_running_fisher(g, g, universe_size = 1000, clip = 500)
  expect_gte(s2, s - 1e-8)
})

test_that("bid_running_fisher is symmetric up to fp noise", {
  set.seed(1)
  pool <- paste0("g", 1:500)
  g1 <- sample(pool, 50); g2 <- sample(pool, 60)
  s12 <- bid_running_fisher(g1, g2, 500)
  s21 <- bid_running_fisher(g2, g1, 500)
  expect_equal(s12, s21, tolerance = 1e-8)
})

test_that("running_fisher equals the full G x K scan (smallest p over all cut-offs)", {
  brute <- function(g1, g2, P) {
    best <- 0
    for (G in seq_along(g1)) for (K in seq_along(g2)) {
      M <- length(intersect(g1[seq_len(G)], g2[seq_len(K)]))
      if (M > 0) best <- max(best, -stats::phyper(M - 1, K, P - K, G, lower.tail = FALSE, log.p = TRUE) / log(10))
    }
    best
  }
  set.seed(2)
  pool <- paste0("g", 1:300)
  for (i in 1:20) {
    g1 <- sample(pool, sample(5:40, 1)); g2 <- c(sample(g1, 4), sample(pool, sample(5:40, 1)))
    g2 <- unique(g2)
    expect_equal(running_fisher(g1, g2, 300), min(brute(g1, g2, 300), 300), tolerance = 1e-8)
  }
})

test_that("an unordered gene set (g2_ranked = FALSE) equals the scan over the top-G cut-offs only", {
  brute_set <- function(g1, set, P) {
    K <- length(unique(set)); best <- 0
    for (G in seq_along(g1)) {
      M <- length(intersect(g1[seq_len(G)], set))
      if (M > 0) best <- max(best, -stats::phyper(M - 1, K, P - K, G, lower.tail = FALSE, log.p = TRUE) / log(10))
    }
    best
  }
  set.seed(4)
  pool <- paste0("g", 1:300)
  for (i in 1:20) {
    g1 <- sample(pool, sample(5:60, 1)); set <- unique(c(sample(g1, 4), sample(pool, sample(5:60, 1))))
    s <- running_fisher(g1, set, 300, g2_ranked = FALSE)
    expect_equal(s, min(brute_set(g1, set, 300), 300), tolerance = 1e-8)
    # the order in which the set is stored does not matter
    expect_identical(running_fisher(g1, sample(set), 300, g2_ranked = FALSE), s)
    expect_identical(running_fisher(g1, sort(set), 300, g2_ranked = FALSE), s)
  }
  expect_equal(running_fisher(character(0), letters, 100, g2_ranked = FALSE), 0)
  expect_equal(running_fisher(letters[1:5], LETTERS, 100, g2_ranked = FALSE), 0)
  expect_equal(running_fisher(letters[1:5], letters[1:5], 1000, clip = 5, g2_ranked = FALSE), 5)
})

test_that("g2_ranked = TRUE is the default ranked kernel", {
  set.seed(5)
  pool <- paste0("g", 1:400)
  for (i in 1:10) {
    g1 <- sample(pool, 80); g2 <- sample(pool, 60)
    expect_identical(running_fisher(g1, g2, 400, g2_ranked = TRUE), running_fisher(g1, g2, 400))
  }
})

test_that("rfind_score_cell uses a gene set whole in all four terms", {
  set.seed(6)
  pool <- paste0("g", 1:500)
  b1u <- sample(pool, 60); b1d <- sample(setdiff(pool, b1u), 60)
  up <- sample(pool, 40); dn <- sample(setdiff(pool, up), 30)
  rf <- function(a, b) running_fisher(a, b, 500, g2_ranked = FALSE)
  b2 <- list(up = up, down = dn, ranked = FALSE)
  expect_equal(rfind_score_cell(b1u, b1d, b2, 500),
               rf(b1u, up) - rf(b1u, dn) - rf(b1d, up) + rf(b1d, dn), tolerance = 1e-12)
  expect_identical(rfind_score_cell(b1u, b1d, list(up = rev(up), down = sort(dn), ranked = FALSE), 500),
                   rfind_score_cell(b1u, b1d, b2, 500))
  # without `ranked` (or ranked = TRUE) the lists are ranked
  rr <- function(a, b) running_fisher(a, b, 500)
  expect_equal(rfind_score_cell(b1u, b1d, list(up = up, down = dn), 500),
               rr(b1u, up) - rr(b1u, dn) - rr(b1d, up) + rr(b1d, dn), tolerance = 1e-12)
})

test_that("to_b2 marks membership gene sets as unranked", {
  side <- function(g) data.frame(gene = g, stringsAsFactors = FALSE)
  expect_true(to_b2(list(up = side("A"), down = side("B"), type = "deg"))$ranked)
  expect_true(to_b2(list(up = side("A"), down = side(character(0)), type = "gene_set"))$ranked)
  expect_true(to_b2(list(up = side("A"), down = side(character(0))))$ranked)
  for (tp in c("gene_group", "msigdb_type_b", "msigdb_cgp_paired"))
    expect_false(to_b2(list(up = side("A"), down = side(character(0)), type = tp))$ranked)
  # an explicit `ranked` field wins over the type
  expect_false(to_b2(list(up = side("A"), down = side(character(0)), type = "gene_set", ranked = FALSE))$ranked)
  expect_true(to_b2(list(up = side("A"), down = side(character(0)), type = "gene_group", ranked = TRUE))$ranked)
  # bundled references: gene groups unranked, signed and ranked UP-only references ranked
  bs <- c(default_biosets(), load_dedonno_biosets()[c("Glucocorticoid", "Seizure", "UPR", "NF-kB")])
  rk <- vapply(setdiff(names(bs), "meta"), function(a) to_b2(bs[[a]])$ranked, NA)
  expect_identical(rk, c(DAM = TRUE, Maturity = TRUE, Interferon = FALSE, G2M = FALSE,
                         Glucocorticoid = TRUE, Seizure = TRUE, UPR = FALSE, `NF-kB` = FALSE))
})

test_that("bid_running_fisher is an alias of running_fisher", {
  set.seed(3)
  pool <- paste0("g", 1:200); g1 <- sample(pool, 30); g2 <- sample(pool, 30)
  expect_identical(bid_running_fisher(g1, g2, 200), running_fisher(g1, g2, 200))
})

test_that("default_biosets loads the 4 Hammond axes of the paper", {
  bs <- default_biosets()
  expect_type(bs, "list")
  expected <- c("DAM", "Maturity", "Interferon", "G2M")
  expect_setequal(setdiff(names(bs), "meta"), expected)
  expect_identical(bs$meta$axes, expected)
  expect_true(is.data.frame(bs$DAM$up))
  expect_true("gene" %in% colnames(bs$DAM$up))
  expect_equal(nrow(bs$DAM$up), 500)
  # gene-group panels have an empty DOWN side
  expect_equal(nrow(bs$Interferon$down), 0)
})

test_that("load_dedonno_biosets loads the 4 De Donno axes of the paper", {
  dd <- load_dedonno_biosets()
  expect_setequal(setdiff(names(dd), "meta"), c("Glucocorticoid", "Seizure", "UPR", "NF-kB"))
  expect_equal(nrow(dd$Seizure$down), 0)
  # labels follow the axis names (0.3.2)
  expect_identical(vapply(dd[setdiff(names(dd), "meta")], `[[`, "", "label"),
                   c(Glucocorticoid = "Glucocorticoid", Seizure = "Seizure", UPR = "UPR", `NF-kB` = "NF-kB"))
})

test_that("build_panel_from_deg splits by sign and caps to top_n", {
  deg <- data.frame(
    gene = paste0("g", 1:30),
    log2FC = c(runif(15, 0.1, 3), -runif(15, 0.1, 3))
  )
  p <- build_panel_from_deg(deg, top_n = 10, label = "test")
  expect_equal(nrow(p$up), 10)
  expect_equal(nrow(p$down), 10)
  expect_true(all(p$up$diff > 0))
  expect_true(all(p$down$diff < 0))
  expect_equal(p$type, "deg")
})

test_that("running_fisher warns when a list names a gene twice and computes the score as before (0.3.9)", {
  g1 <- c("a", "a", "b", "c"); g2 <- c("a", "d", "b")
  expect_warning(s <- running_fisher(g1, g2, universe_size = 50), class = "rfindsc_duplicated_names")
  expect_true(is.finite(s))
  expect_silent(running_fisher(unique(g1), g2, universe_size = 50))
})
