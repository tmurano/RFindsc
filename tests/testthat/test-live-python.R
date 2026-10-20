test_that("current R and Python sources agree on shared inputs", {
  python_source <- Sys.getenv("RFINDSC_PY_SOURCE")
  skip_if(!nzchar(python_source), "set RFINDSC_PY_SOURCE to run the live cross-language test")
  python_source <- normalizePath(python_source, mustWork = TRUE)
  script <- file.path(python_source, "tools/check_cross_language.py")
  expect_true(file.exists(script))
  python <- Sys.getenv("RFINDSC_PYTHON", unset = Sys.which("python3"))
  report <- tempfile(fileext = ".json")
  output <- system2(python, c(shQuote(script), "--r-source", shQuote(normalizePath(test_path("..", ".."))),
                             "--report", shQuote(report)), stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  expect_true(is.null(status) || status == 0L, info = paste(output, collapse = "\n"))
  expect_true(file.exists(report))
  if (file.exists(report)) {
    result <- jsonlite::fromJSON(report)
    expect_true(result$passed)
    expect_lt(result$score_max_abs_diff, 1e-8)
    expect_lt(result$kernel_max_abs_diff, 1e-8)
    expect_true(all(result$cases$signatures_equal))
  }
})
