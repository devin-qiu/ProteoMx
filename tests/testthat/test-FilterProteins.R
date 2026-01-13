test_that("FilterProteins subsets the object", {
  obj <- readRDS(test_path("testdata/tiny_geomx.rds"))
  
  # Setup: Must have q_norm and Fits
  if (!"q_norm" %in% names(assayData(obj))) obj <- Q3Normalize(obj)
  obj <- MixModelFit(obj, ncomps = 2) # Run fresh
  
  # Act: Filter
  # We use Rt IgG2a because we forced it to be in the subset
  res <- FilterProteins(obj, neg_ctrl = "Rt IgG2a")
  
  # Expect
  expect_s4_class(res, "NanoStringGeoMxSet")
  
  # The Negative Control itself MUST always be kept
  expect_true("Rt IgG2a" %in% rownames(res))
  
  # The resulting object should have <= rows than input
  expect_lte(nrow(res), nrow(obj))
})