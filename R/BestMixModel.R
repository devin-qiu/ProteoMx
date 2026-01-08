#' Select Best Mixture Model (Manual Tournament)
#'
#' Manually fits models for ncomp = 1 to ncomps and picks the winner (Lowest BIC).
#' This bypasses mixR::select() to avoid crashes/NAs.
#'
#' @param geomx_set A GeoMxSet object containing @assayData$q_norm.
#' @param ncomps Integer. Max components to test (default 3).
#' @return A tibble with columns: Protein, Best_NComp, Best_EV, Best_BIC.
#' @export
BestMixModel <- function(geomx_set, ncomps = 3) {
  
  require(mixR)
  require(dplyr)
  require(tibble)
  require(Biobase)
  
  # --- 1. Data Prep ---
  if (!"q_norm" %in% names(assayData(geomx_set))) stop("Error: 'q_norm' missing.")
  q3_mat <- assayDataElement(geomx_set, "q_norm")
  
  # Clean columns
  bad_cols <- apply(q3_mat, 2, function(x) all(is.na(x) | is.infinite(x)))
  if (sum(bad_cols) > 0) q3_mat <- q3_mat[, !bad_cols]
  
  log_mat <- log10(q3_mat + 1)
  all_proteins <- rownames(log_mat)
  
  # Output containers
  best_ncomp_sel <- rep(NA_integer_, length(all_proteins))
  best_ev_sel    <- rep(NA, length(all_proteins)) 
  best_bic_sel   <- rep(NA_real_, length(all_proteins))
  
  message("Running Manual Model Tournament (1 to ", ncomps, " components)...")
  pb <- txtProgressBar(min = 0, max = length(all_proteins), style = 3)
  
  # --- 2. Loop Over Proteins ---
  for (i in seq_along(all_proteins)) {
    setTxtProgressBar(pb, i)
    
    vals <- as.vector(as.numeric(log_mat[i, ]))
    vals <- vals[is.finite(vals)]
    
    if (length(vals) < 10) next 
    
    # --- 3. The Tournament (Fit 1..N and compare) ---
    best_bic_local <- Inf
    best_model_local <- list(ncomp = NA, ev = NA)
    
    # We test BOTH Equal Variance (TRUE) and Unequal Variance (FALSE)
    for (ev_flag in c(FALSE, TRUE)) {
      for (k in 1:ncomps) {
        
        # Use do.call on the *working* mixfit function
        fit <- tryCatch(
          do.call(mixR::mixfit, list(x = vals, ncomp = k, family = "norm", ev = ev_flag, init.method = "kmeans")),
          error = function(e) NULL
        )
        
        if (!is.null(fit) && !is.null(fit$bic) && !is.nan(fit$bic)) {
          if (fit$bic < best_bic_local) {
            best_bic_local <- fit$bic
            best_model_local <- list(ncomp = k, ev = ev_flag)
          }
        }
      }
    }
    
    # Store the winner for this protein
    if (best_bic_local < Inf) {
      best_ncomp_sel[i] <- best_model_local$ncomp
      best_ev_sel[i]    <- best_model_local$ev
      best_bic_sel[i]   <- best_bic_local
    }
  }
  close(pb)
  
  # --- 4. Format Output ---
  tibble(
    Protein    = all_proteins,
    Best_NComp = best_ncomp_sel,
    Best_EV    = best_ev_sel,
    Best_BIC   = best_bic_sel
  )
}
