#' Decay Background Signal with Dynamic Antibody Scaling
#'
#' Applies an exponential decay penalty to expression values below the negative 
#' control threshold. The intensity of the penalty is dynamically scaled for each 
#' protein based on the geometric slope (k) of its left-tail noise distribution.
#' Narrow noise peaks receive a steeper penalty, while wide/sparse noise distributions 
#' receive a milder penalty.
#'
#' @param geomx_set A \code{NanoStringGeoMxSet} object.
#' @param neg_ctrl Character. The negative control probe to calculate the threshold.
#' @param n_sd Numeric. Number of standard deviations for the threshold (default: 1).
#' @param base_decay Numeric. The baseline intensity of the penalty before antibody-specific scaling (default: 2.0).
#'
#' @return A \code{NanoStringGeoMxSet} object with a new assay added ('exp_decayed') containing the penalized counts.
#' 
#' @importFrom Biobase assayDataElement assayDataElement<- assayData
#' @importFrom stats density sd
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
  message("Calculating dynamic antibody-specific penalties (k)...")
  
  # --- 3. Calculate Dynamic Slope (k) for each protein ---
  n_proteins <- nrow(log_mat)
  k_values <- numeric(n_proteins)
  names(k_values) <- rownames(log_mat)
  
  for (i in seq_len(n_proteins)) {
    vals <- as.numeric(log_mat[i, ])
    vals <- vals[is.finite(vals)]
    
    # Skip if perfectly uniform (e.g. all zeroes) to prevent density() crash
    if (length(unique(vals)) < 3) {
      k_values[i] <- 0.1 # Minimal penalty for near-empty data
      next
    }
    
    # Fit density curve
    dens <- density(vals, na.rm = TRUE)
    
    # Find the Y value at the Threshold (intersection point)
    # We find the X in the density object closest to our threshold
    idx_thresh <- which.min(abs(dens$x - threshold))
    x_thresh <- dens$x[idx_thresh]
    y_thresh <- dens$y[idx_thresh]
    
    # Find the Y value at the start of the left tail
    # We use the minimum observed value (or 0) as the anchor
    min_x_val <- max(0, min(vals))
    idx_start <- which.min(abs(dens$x - min_x_val))
    x_start <- dens$x[idx_start]
    y_start <- dens$y[idx_start]
    
    # Calculate Geometric Slope (delta Y / delta X)
    delta_x <- x_thresh - x_start
    delta_y <- y_thresh - y_start
    
    if (delta_x <= 0) {
      # If threshold is lower than all data points, there is no left tail.
      k_values[i] <- 0.1 
    } else {
      # Absolute value ensures k is positive.
      raw_slope <- abs(delta_y / delta_x)
      k_values[i] <- raw_slope
    }
  }
  
  # Normalize k values so the median antibody has a scalar of 1.0.
  # This prevents extreme slopes from pushing the math to infinity.
  median_k <- median(k_values[k_values > 0], na.rm = TRUE)
  if (is.na(median_k) || median_k == 0) median_k <- 1
  k_values_normalized <- k_values / median_k
  
  # --- 4. Apply Dynamic Exponential Decay ---
  decayed_log_mat <- log_mat
  
  for (i in seq_len(n_proteins)) {
    prot_name <- rownames(decayed_log_mat)[i]
    k_scalar <- k_values_normalized[prot_name]
    
    # Mask for this specific protein's noise
    mask <- decayed_log_mat[i, ] < threshold & decayed_log_mat[i, ] > 0
    
    # Apply dynamically scaled formula: x * e^(base_decay * k * (x - threshold))
    x_vals <- decayed_log_mat[i, mask]
    decayed_log_mat[i, mask] <- x_vals * exp(base_decay * k_scalar * (x_vals - threshold))
  }
  
  # --- 5. Convert back to linear scale ---
  decayed_mat <- (10 ^ decayed_log_mat) - 1
  decayed_mat[decayed_mat < 0] <- 0
  
  # --- 6. Save securely as a new assay ---
  assayDataElement(geomx_set, "exp_decayed") <- decayed_mat
  
  message("Successfully applied dynamically scaled exponential decay.")
  message("Decayed data safely stored in new assay: 'exp_decayed'")
  
  # Optional: Store the k_values in the experimentData for user transparency
  experimentData(geomx_set)@other$Antibody_Decay_Scalars <- k_values_normalized
  
  return(geomx_set)
}