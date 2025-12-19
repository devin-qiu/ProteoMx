#' Mixture Model Fit (Batch Engine)
#'
#' Fits mixture models for all proteins across all specified combinations of components 
#' and equal variance assumptions. Does NOT run model selection.
#'
#' @param geomx_set A GeoMxSet object containing @assayData$q_norm.
#' @param ncomps Integer. Maximum number of components to fit (default 3).
#' @param neg_ctrl Character. Name of the negative control probe (default "Rt IgG2a").
#'
#' @return A list containing:
#'   \item{fits}{A list of tibbles, keys are "evX_ncomp_Y". Matches original script structure.}
#'   \item{data}{The log10-transformed data matrix.}
#'   \item{threshold}{The background threshold used.}
#' @export
MixModelFit <- function(geomx_set, ncomps = 3, neg_ctrl = "Rt IgG2a") {
  
  require(mixR)
  require(dplyr)
  require(tidyr)
  require(tibble)
  require(utils)
  
  # --- 1. Data Prep ---
  q3_mat <- tryCatch(
    geomx_set@assayData$q_norm,
    error = function(e) stop("Could not access @assayData$q_norm.")
  )
  
  # Clean columns
  bad_cols <- apply(q3_mat, 2, function(x) all(is.na(x) | is.infinite(x)))
  if (sum(bad_cols) > 0) q3_mat <- q3_mat[, !bad_cols]
  
  log_mat <- log10(q3_mat + 1)
  
  if (!neg_ctrl %in% rownames(log_mat)) stop("Negative control '", neg_ctrl, "' not found.")
  
  # --- 2. Threshold ---
  neg_vals <- as.numeric(log_mat[neg_ctrl, ])
  neg_vals <- neg_vals[is.finite(neg_vals)]
  threshold <- mean(neg_vals) + 1 * sd(neg_vals)
  message("Threshold calculated: ", round(threshold, 3))
  
  # --- 3. Internal Helper ---
  fit_single_protein <- function(vals, n, ev) {
    vals <- vals[is.finite(vals)]
    if (length(vals) < 10) return(NULL)
    
    fit_result <- tryCatch(
      mixR::mixfit(vals, ncomp = n, family = "norm", ev = ev, init.method = "kmeans"),
      error = function(e) NULL
    )
    
    if (is.null(fit_result)) return(NULL)
    
    # Shift means to match data scale
    mix_mean <- sum(fit_result$pi * fit_result$mu)
    fit_result$mu <- fit_result$mu + (mean(vals) - mix_mean)
    return(fit_result)
  }
  
  # --- 4. Main Loop ---
  all_proteins <- rownames(log_mat)
  ev_values <- c(FALSE, TRUE)
  fit_results_all <- list()
  
  total_steps <- length(ev_values) * length(1:ncomps)
  step_count <- 0
  
  message("Starting Batch Fit...")
  
  for (ev_flag in ev_values) {
    for (n in 1:ncomps) {
      step_count <- step_count + 1
      message(paste0("[", step_count, "/", total_steps, "] Fitting ev=", ev_flag, ", ncomp=", n, "..."))
      
      protein_fits <- vector("list", length(all_proteins))
      keep_flags   <- logical(length(all_proteins))
      
      pb <- txtProgressBar(min = 0, max = length(all_proteins), style = 3)
      
      for (i in seq_along(all_proteins)) {
        setTxtProgressBar(pb, i)
        vals <- as.numeric(log_mat[i, ]) # Access by index is faster
        
        fit <- fit_single_protein(vals, n, ev_flag)
        
        if (!is.null(fit)) {
          protein_fits[[i]] <- fit
          keep_flags[i] <- max(fit$mu, na.rm = TRUE) > threshold
        } else {
          keep_flags[i] <- FALSE
        }
      }
      close(pb)
      
      # Construct the tibble for this bucket
      fits_tbl <- tibble(
        Protein    = all_proteins,
        fit_result = protein_fits,
        keep       = keep_flags
      )
      
      key <- paste0("ev", if (ev_flag) "T" else "F", "_ncomp_", n)
      fit_results_all[[key]] <- fits_tbl
    }
  }
  
  return(list(
    fits = fit_results_all,
    data = log_mat,
    threshold = threshold
  ))
}