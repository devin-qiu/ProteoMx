#' Filter Proteins by Best Mixture Model Detection
#'
#' Subsets the GeoMxSet object to retain only proteins that have a detectable signal
#' above the negative control background in their statistically optimal mixture model.
#'
#' @description
#' This function uses the "Best Model" selected by \code{\link{BestMixModel}} (via Bayes Factors).
#' It compares the signal of the best model against a dynamic negative control threshold.
#' Proteins are retained only if their optimal model's highest component mean exceeds the threshold.
#' Negative control probes are mathematically excluded but silently preserved in the output object
#' for downstream plotting functions.
#'
#' @param geomx_set A \code{NanoStringGeoMxSet} object processed by \code{\link{BestMixModel}}.
#' @param neg_ctrl Character. The negative control probe to use as the noise baseline.
#'   The threshold is calculated as \code{mean(neg_ctrl) + n_sd * SD(neg_ctrl)}.
#' @param n_sd Numeric. The number of standard deviations above the mean to set the threshold (default: 1).
#'
#' @return A subsetted \code{NanoStringGeoMxSet} object containing the detected proteins and negative controls.
#'
#' @importFrom Biobase assayDataElement experimentData
#' @importFrom dplyr filter
#' @export
FilterProteins <- function(geomx_set, neg_ctrl = "Rt IgG2a", n_sd = 1) {
  
  require(dplyr)
  require(Biobase)
  
  # --- 1. Validation Checks ---
  best_model_df <- experimentData(geomx_set)@other$Best_Model_Summary
  
  if (is.null(best_model_df)) {
    stop("BestMixModel summary not found in experimentData. Please run geomx_set <- BestMixModel(geomx_set) first.")
  }
  
  if (!"q_norm" %in% names(assayData(geomx_set))) stop("q_norm data missing. Run Q3Normalize() first.")
  q3_mat <- assayDataElement(geomx_set, "q_norm")
  
  # --- 2. Calculate Dynamic Threshold ---
  log_mat <- log10(q3_mat + 1)
  
  if (!neg_ctrl %in% rownames(log_mat)) {
    stop("Negative control '", neg_ctrl, "' not found in dataset.")
  }
  
  neg_vals <- as.numeric(log_mat[neg_ctrl, ])
  neg_vals <- neg_vals[is.finite(neg_vals)]
  
  threshold <- mean(neg_vals, na.rm = TRUE) + (n_sd * sd(neg_vals, na.rm = TRUE))
  message("Filtering Threshold (", neg_ctrl, " + ", n_sd, "SD): ", round(threshold, 3))
  
  # --- 3. Identify Proteins to Keep ---
  all_proteins <- rownames(log_mat)
  
  # Identify all negative controls in the dataset
  common_negs <- unique(c("Hmr IgG", "Ms IgG2b", "Ms IgG2a", "Rb IgG", "Rt IgG2a", "Ms IgG1", neg_ctrl))
  negs_in_data <- intersect(all_proteins, common_negs)
  
  passed_proteins <- c()
  
  # Iterate through the Best Model Summary table
  for (i in seq_len(nrow(best_model_df))) {
    prot_name <- best_model_df$Protein[i]
    signal_mean <- best_model_df$Best_Signal_Mean[i]
    
    if (is.na(signal_mean)) next
    
    if (signal_mean > threshold) {
      passed_proteins <- c(passed_proteins, prot_name)
    }
  }
  
  passed_proteins <- unique(passed_proteins)
  
  # FORCE negative controls out of the 'passed' list (even if they mathematically crossed the threshold)
  passed_proteins <- setdiff(passed_proteins, negs_in_data)
  
  # Failed proteins are everything else (which now strictly includes the negative controls)
  failed_proteins <- setdiff(all_proteins, passed_proteins)
  
  # The actual list we will use to subset the object (Biological Passes + Negative Controls)
  object_proteins <- unique(c(passed_proteins, negs_in_data))
  
  # --- 4. Report Results ---
  n_total <- length(all_proteins)
  n_kept  <- length(passed_proteins)
  n_excl  <- length(failed_proteins)
  
  message("\n==========================================")
  message(" FILTERING SUMMARY (Strict Bayes Logic)")
  message("==========================================")
  message("Total Proteins Input: ", n_total)
  message("Proteins Retained (Biological): ", n_kept)
  message("Proteins Excluded (Includes Neg Probes): ", n_excl)
  
  if (n_excl > 0) {
    message("\n--- EXCLUDED PROTEINS (Best Model Signal < Threshold) ---")
    cat(paste(failed_proteins, collapse = ", "), "\n")
  }
  
  if (n_kept > 0) {
    message("\n--- RETAINED PROTEINS ---")
    cat(paste(passed_proteins, collapse = ", "), "\n")
  }
  
  message("\n* Note: Negative control probes are listed as 'Excluded' but are silently preserved in the output object for plotting.")
  message("==========================================\n")
  
  # --- 5. Execute Subset ---
  geomx_subset <- geomx_set[object_proteins, ]
  
  return(geomx_subset)
}