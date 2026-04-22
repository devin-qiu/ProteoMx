#' Decay Background Signal with ECDF Confidence Weighting
#'
#' Applies an exponential decay penalty to expression values below the negative 
#' control threshold. The intensity of the penalty is dynamically scaled point-by-point 
#' using the Empirical Cumulative Distribution Function (ECDF). Data points close to 
#' the threshold are spared (high confidence), while points deep in the left tail 
#' are aggressively decayed (low confidence).
#'
#' @param geomx_set A \code{NanoStringGeoMxSet} object.
#' @param neg_ctrl Character. The negative control probe to calculate the threshold.
#' @param n_sd Numeric. Number of standard deviations for the threshold (default: 1).
#' @param base_decay Numeric. The baseline intensity of the penalty (default: 2.0).
#'
#' @return A \code{NanoStringGeoMxSet} object with a new assay added ('exp_decayed') containing the penalized counts.
#'
#' @details
#' Traditional hard thresholding (setting sub-threshold values to zero or NA) creates 
#' artificial zero-variance spikes that violate the assumptions of downstream linear 
#' modeling (e.g., limma). To address this, \code{DecayBackground} applies a continuous 
#' soft-thresholding transformation.
#' 
#' For values below the threshold (T), the data undergoes an exponential decay:
#' \preformatted{
#'   X_new = X * exp(r * k * (X - T))
#' }
#' Where 'r' is the \code{base_decay} rate and 'k' is a dynamic confidence weight derived 
#' from the ECDF (F) of that specific protein across all ROIs:
#' \preformatted{
#'   k = 1 - (F(X) / F(T))
#' }
#' This approach creates a "buffer zone": it protects borderline data points just below 
#' the threshold (where k approaches 0) while aggressively suppressing unambiguous 
#' background noise deeper in the tail (where k approaches 1). This scales organically 
#' to individual antibody kinetics and beautifully preserves true biological variance.
#' 
#' @importFrom Biobase assayDataElement assayDataElement<- assayData
#' @importFrom stats ecdf sd
#' @export
DecayBackground <- function(geomx_set, neg_ctrl, n_sd = 1, base_decay = 2.0) {
  
  require(Biobase)
  require(stats)
  
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
  message("Applying ECDF-weighted dynamic penalties...")
  
  # --- 3. Apply Dynamic ECDF Exponential Decay ---
  decayed_log_mat <- log_mat
  n_proteins <- nrow(log_mat)
  
  for (i in seq_len(n_proteins)) {
    
    # Extract all valid values to build the full distribution profile
    vals <- as.numeric(log_mat[i, ])
    vals <- vals[is.finite(vals)]
    
    # Skip if perfectly uniform to prevent ECDF failure
    if (length(unique(vals)) < 3) next
    
    # Generate the ECDF for this specific protein
    f_ecdf <- ecdf(vals)
    F_T <- f_ecdf(threshold)
    
    # Only apply penalty if there is data below the threshold
    if (F_T > 0) {
      mask <- decayed_log_mat[i, ] < threshold & decayed_log_mat[i, ] > 0
      x_vals <- decayed_log_mat[i, mask]
      
      # The Confidence Penalty Scalar: 1 - (ECDF(x) / ECDF(Threshold))
      # Near threshold -> k approaches 0 -> penalty approaches 1 (preserved)
      # Deep noise -> k approaches 1 -> full exponent applied (squashed)
      k_scalar <- 1 - (f_ecdf(x_vals) / F_T)
      
      decayed_log_mat[i, mask] <- x_vals * exp(base_decay * k_scalar * (x_vals - threshold))
    }
  }
  
  # --- 4. Convert back to linear scale ---
  decayed_mat <- (10 ^ decayed_log_mat) - 1
  decayed_mat[decayed_mat < 0] <- 0
  
  # --- 5. Save securely as a new assay ---
  assayDataElement(geomx_set, "exp_decayed") <- decayed_mat
  
  message("Successfully applied ECDF-weighted exponential decay.")
  message("Decayed data safely stored in new assay: 'exp_decayed'")
  
  return(geomx_set)
}