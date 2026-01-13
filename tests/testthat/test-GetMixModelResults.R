test_that("GetMixModelResults retrieves list or errors if missing", {
  obj <- readRDS(test_path("testdata/tiny_geomx.rds"))
  
  # 1. Test Error: Remove any existing results
  experimentData(obj)@other$MixModel <- NULL
  expect_error(GetMixModelResults(obj), "No Mixture Model results found")
  
  # 2. Test Success: Inject a dummy result (faster than running fit)
  dummy_res <- list(fits = list(), params = list(test = TRUE))
  experimentData(obj)@other$MixModel <- dummy_res
  
  out <- GetMixModelResults(obj)
  expect_equal(out$params$test, TRUE)
})
