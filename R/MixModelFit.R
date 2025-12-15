#' Batch Mixture Model Fit for All Proteins
#'
#' @param geomx_set A GeoMxSet object (containing @assayData$q_norm).
#' @param ncomps Integer. The maximum number of components to test (default 3).
#' @param neg_ctrl Character string. Name of the negative control target (default "Rt IgG2a").
#'
#' @return A list containing:
#'   \item{fits}{A named list where keys are protein names and values are the fit objects.}
#'   \item{data}{The cleaned, log10-transformed data matrix (for plotting).}
#'   \item{threshold}{The numeric background threshold used.}
#'   \item{failed}{A vector of protein names that failed to fit (e.g. low counts).}
#' @export
MixModelFit <- function(geomx_set, ncomps = 3, neg_ctrl = "Rt IgG2a") {
  
  require(mixR)
  require(dplyr)
  require(tidyr)
  require(tibble)
  require(utils) # For progress bar
  
  # --- 1. Data Extraction & Cleaning ---
  q3_mat <- tryCatch(
    geomx_set@assayData$q_norm,
    error = function(e) stop("Could not access @assayData$q_norm.")
  )
  
  # Remove "Bad AOIs" (Columns entirely NaN/Inf)
  bad_cols <- apply(q3_mat, 2, function(x) all(is.na(x) | is.infinite(x)))
  if (sum(bad_cols) > 0) {
    message("Cleaning: Removing ", sum(bad_cols), " AOI columns containing only NaN/Inf.")
    q3_mat <- q3_mat[, !bad_cols]
  }
  
  if (!neg_ctrl %in% rownames(q3_mat)) stop("Negative control '", neg_ctrl, "' not found.")
  
  # Prepare Data Matrix (Log10 scale)
  # We do this once for the whole matrix to save time
  log_mat <- log10(q3_mat + 1)
  
  # --- 2. Calculate Global Threshold ---
  neg_vals <- as.numeric(log_mat[neg_ctrl, ])
  neg_vals <- neg_vals[is.finite(neg_vals)]
  
  bg_threshold_log <- mean(neg_vals) + 1 * sd(neg_vals)
  message("Background Threshold (Log10): ", round(bg_threshold_log, 3))
  
  # --- 3. Loop Over All Proteins ---
  all_proteins <- rownames(log_mat)
  fit_results  <- list()
  failed_prots <- c()
  
  # Initialize Progress Bar
  pb <- txtProgressBar(min = 0, max = length(all_proteins), style = 3)
  
  for (i in seq_along(all_proteins)) {
    pname <- all_proteins[i]
    setTxtProgressBar(pb, i)
    
    vals <- as.numeric(log_mat[pname, ])
    vals <- vals[is.finite(vals)]
    
    # Skip if insufficient data
    if (length(vals) < 10 || length(unique(vals)) < 3) {
      failed_prots <- c(failed_prots, pname)
      next
    }
    
    # -- Fit all requested models for this protein --
    protein_fits <- list()
    summary_list <- list()
    
    ev_settings <- c(FALSE, TRUE)
    
    for (ev in ev_settings) {
      for (k in 1:ncomps) {
        
        # Unique key (e.g., "evF_ncomp_2")
        key <- paste0("ev", if(ev) "T" else "F", "_ncomp_", k)
        
        fit <- tryCatch(
          mixR::mixfit(vals, ncomp = k, family = "normal", ev = ev, init.method = "kmeans"),
          error = function(e) NULL
        )
        
        if (!is.null(fit)) {
          # Shift means for plotting alignment (mixR quirk)
          mix_mean <- sum(fit$pi * fit$mu)
          fit$mu   <- fit$mu + (mean(vals) - mix_mean)
          
          # Store fit
          protein_fits[[key]] <- fit
          
          # Store stats for selection
          summary_list[[key]] <- tibble(
            ncomp = k,
            equal_var = ev,
            bic = fit$BIC
          )
        }
      }
    }
    
    # Identify Best Model (min BIC) if any fits succeeded
    best_model_info <- NULL
    if (length(summary_list) > 0) {
      df <- bind_rows(summary_list) %>% arrange(bic)
      best_row <- df[1, ]
      best_model_info <- list(ncomp = best_row$ncomp, ev = best_row$equal_var, bic = best_row$bic)
    }
    
    # Save results for this protein
    if (length(protein_fits) > 0) {
      fit_results[[pname]] <- list(
        all_fits = protein_fits,
        best_model = best_model_info
      )
    } else {
      failed_prots <- c(failed_prots, pname)
    }
  }
  
  close(pb)
  message("\nProcessed ", length(all_proteins), " proteins. ", length(fit_results), " succeeded, ", length(failed_prots), " failed/skipped.")
  
  # Return the Massive Object
  return(list(
    fits = fit_results,
    data = log_mat,         # Store the whole matrix for plotting later
    threshold = bg_threshold_log,
    failed_proteins = failed_prots
  ))
}