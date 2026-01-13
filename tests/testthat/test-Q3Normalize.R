test_that("Q3Normalize adds q_norm assayData element", {
  # Load real data
  obj <- readRDS(test_path("testdata/tiny_geomx.rds"))
  
  # SIMULATE RAW DATA: Remove 'q_norm' if it exists so we can test adding it
  if ("q_norm" %in% names(assayData(obj))) {
    assayDataElement(obj, "q_norm") <- NULL
  }
  
  # Expectation: It shouldn't be there yet
  expect_false("q_norm" %in% names(assayData(obj)))
  
  # Act: Run normalization
  res <- Q3Normalize(obj)
  
  # Expectation: Now it should be there
  expect_true("q_norm" %in% names(assayData(res)))
  expect_s4_class(res, "NanoStringGeoMxSet")
})