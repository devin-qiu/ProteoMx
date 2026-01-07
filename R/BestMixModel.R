#' Select Best Mixture Model
#'
#' Runs mixR::select() on every protein to find the optimal number of components
#' and variance structure (BIC-based).
#'
#' @param geomx_set A GeoMxSet object.
#' @param ncomps Integer. Max components to test (default 3).
#'
#' @return A tibble with columns: Protein, Best_NComp, Best_EV, Best_BIC.
#' @export
BestMixModel <- function(geomx_set, ncomps = 3) {
  
  require(mixR)
  require(dplyr)
  require(tibble)
  require(utils)
  
  # --- Data Prep ---
  q3_mat <- tryCatch(
    geomx_set@assayData$q_norm,
    error = function(e) stop("Could not access @assayData$q_norm.")
  )
  bad_cols <- apply(q3_mat, 2, function(x) all(is.na(x) | is.infinite(x)))
  if (sum(bad_cols) > 0) q3_mat <- q3_mat[, !bad_cols]
  log_mat <- log10(q3_mat + 1)
  
  all_proteins <- rownames(log_mat)
  
  # Initialize vectors
  best_ncomp_sel <- rep(NA_integer_, length(all_proteins))
  best_ev_sel    <- rep(NA, length(all_proteins)) # Logical or Char
  best_bic_sel   <- rep(NA_real_, length(all_proteins))
  
  message("Running Model Selection (mixR::select) for ", length(all_proteins), " proteins...")
  pb <- txtProgressBar(min = 0, max = length(all_proteins), style = 3)
  
  for (i in seq_along(all_proteins)) {
    setTxtProgressBar(pb, i)
    
    vals <- as.numeric(log_mat[i, ])
    vals <- vals[is.finite(vals)]
    
    # Needs enough data points
    if (length(vals) < 10) next 
    
    sel <- tryCatch(
      mixR::select(vals, ncomp = 1:ncomps, family = "normal", 
                   mstep.method = "newton", init.method = "kmeans", verbose = FALSE),
      error = function(e) NULL
    )
    
    if (!is.null(sel)) {
      # The 'best' column contains "*" for the winner.
      idx <- which(sel$best == "*")
      
      # Fallback: if "*" is missing, take min BIC
      if (length(idx) == 0) idx <- which.min(sel$bic)
      
      # Extract
      if (length(idx) > 0) {
        best_ncomp_sel[i] <- sel$ncomp[idx]
        best_ev_sel[i]    <- sel$equal.var[idx]
        best_bic_sel[i]   <- sel$bic[idx]
      }
    }
  }
  close(pb)
  
  # Format output
  best_model_df <- tibble(
    Protein    = all_proteins,
    Best_NComp = best_ncomp_sel,
    Best_EV    = best_ev_sel,
    Best_BIC   = best_bic_sel
  )
  
  return(best_model_df)
}