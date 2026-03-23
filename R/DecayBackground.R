#' Decay Background Signal (Soft Thresholding)
#'
#' Applies an exponential decay penalty to expression values that fall below 
#' the negative control threshold. This gracefully reduces the impact of background noise
#' on downstream differential expression without destroying the natural variance 
#' structure by forcing a hard floor.
#'
#' @param geomx_set A \code{NanoStringGeoMxSet} object.
#' @param neg_ctrl Character. The negative control probe to calculate the threshold.
#' @param n_sd Numeric. Number of standard deviations for the threshold (default: 1).
#' @param decay_rate Numeric. The intensity of the penalty. Higher = stronger squashing (default: 2).
#'
#' @return A \code{NanoStringGeoMxSet} object with a new assay added ('exp_decayed') containing the penalized counts.
#' 
#' @importFrom Biobase assayDataElement assayDataElement<- assayData
#' @export
DecayBackground <- function(geomx_set, neg_ctrl, n_sd = 1, decay_rate = 2) {
  
  require(Biobase)
  
  if (!"q_norm" %in% names(assayData(geomx_set))) {
    stop("Input assay 'q_norm' not found. Please ensure your object has been Q3 normalized.")
  }
  
  # --- 1. Extract Data ---
  raw_mat <- assayDataElement(geomx_set, "q_norm")
  log_mat <- log10(raw_mat + 1)
  
  if (!neg_ctrl %in% rownames(log_mat)) {
    stop(paste("Negative control '", neg_ctrl, "' not found in dataset."))
  }
  
  # --- 2. Calculate Threshold ---
  neg_vals <- as.numeric(log_mat[neg_ctrl, ])
  neg_vals <- neg_vals[is.finite(neg_vals)]
  threshold <- mean(neg_vals, na.rm = TRUE) + (n_sd * sd(neg_vals, na.rm = TRUE))
  
  message("Calculated Soft Threshold (", neg_ctrl, " + ", n_sd, "SD): ", round(threshold, 3))
  
  # --- 3. Apply Exponential Decay (Soft Thresholding) ---
  decayed_log_mat <- log_mat
  
  # Create a mask to only target biological noise (values greater than 0 but below threshold)
  mask <- decayed_log_mat < threshold & decayed_log_mat > 0
  
  # Apply the decay formula: x_new = x * e^(decay * (x - threshold))
  decayed_log_mat[mask] <- decayed_log_mat[mask] * exp(decay_rate * (decayed_log_mat[mask] - threshold))
  
  # --- 4. Convert back to linear scale ---
  decayed_mat <- (10 ^ decayed_log_mat) - 1
  
  # Catch any microscopic negative values caused by floating-point math limits
  decayed_mat[decayed_mat < 0] <- 0
  
  # --- 5. Save securely as a new assay ---
  assayDataElement(geomx_set, "exp_decayed") <- decayed_mat
  
  message("Successfully applied exponential decay penalty.")
  message("Decayed data safely stored in new assay: 'exp_decayed'")
  
  return(geomx_set)
}
