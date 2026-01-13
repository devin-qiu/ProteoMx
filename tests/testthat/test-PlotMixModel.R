test_that("PlotMixModel produces a ggplot object", {
  obj <- readRDS(test_path("testdata/tiny_geomx.rds"))
  if (!"q_norm" %in% names(assayData(obj))) obj <- Q3Normalize(obj)
  obj <- MixModelFit(obj, ncomps = 2)
  
  # Act: Plot the negative control (since we know it exists)
  p <- PlotMixModel(obj, protein = "GAPDH", neg_ctrl = "Rt IgG2a")
  
  # Expect
  expect_s3_class(p, "ggplot")
})