test_that("BestMixModel returns valid tibble with selection log", {
  obj <- readRDS(test_path("testdata/tiny_geomx.rds"))
  if (!"q_norm" %in% names(assayData(obj))) obj <- Q3Normalize(obj)
  
  # Act
  res <- BestMixModel(obj, ncomps = 2)
  
  # Expect
  expect_s3_class(res, "tbl_df")
  
  # FIX: Use as.integer() to strip any "Features" label from the S4 object row count
  expect_equal(nrow(res), as.integer(nrow(obj))) 
  
  # Check columns
  expected_cols <- c("Protein", "Best_NComp", "Best_EV", "Best_BIC", "Selection_Log")
  expect_true(all(expected_cols %in% colnames(res)))
  
  # Check internal log structure
  expect_s3_class(res$Selection_Log[[1]], "tbl_df")
})
