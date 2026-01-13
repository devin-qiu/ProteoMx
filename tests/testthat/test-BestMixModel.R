test_that("BestMixModel returns valid tibble with selection log", {
  obj <- readRDS(test_path("testdata/tiny_geomx.rds"))
  if (!"q_norm" %in% names(assayData(obj))) obj <- Q3Normalize(obj)
  
  # Act
  res <- BestMixModel(obj, ncomps = 2)
  
  # Expect
  expect_s3_class(res, "tbl_df")
  expect_equal(nrow(res), nrow(obj)) # One row per protein
  
  # Check columns
  expect_true("Selection_Log" %in% colnames(res))
  
  # Check internal log structure
  expect_s3_class(res$Selection_Log[[1]], "tbl_df")
})
