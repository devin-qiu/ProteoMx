test_that("MixModelFit calculates models and updates metadata", {
  obj <- readRDS(test_path("testdata/tiny_geomx.rds"))
  
  # Ensure q_norm exists (in case your saved subset didn't have it)
  if (!"q_norm" %in% names(assayData(obj))) {
    obj <- Q3Normalize(obj)
  }
  
  # Run Fit (Fast on 10 proteins)
  res <- MixModelFit(obj, ncomps = 2)
  
  # Check Metadata
  meta <- experimentData(res)@other$MixModel
  expect_type(meta, "list")
  expect_named(meta, c("fits", "params", "timestamp"))
  
  # Check if we got results for "Rt IgG2a"
  fits_evT_n2 <- meta$fits[["evT_ncomp_2"]]
  expect_true("Rt IgG2a" %in% fits_evT_n2$Protein)
})
