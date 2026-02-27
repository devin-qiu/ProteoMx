#' Filter Proteins by Best Mixture Model Detection
#'
#' Subsets the GeoMxSet object to retain only proteins that have a detectable signal
#' above the negative control background in their statistically optimal mixture model.
#'
#' @description
#' This function uses the "Best Model" selected by \code{\link{BestMixModel}} (via Bayes Factors).
#' It compares the signal of the best model against a dynamic negative control threshold.
#' Proteins are retained only if their optimal model's highest component mean exceeds the threshold.
#'
#' @param geomx_set A \code{NanoStringGeoMxSet} object processed by \code{\link{BestMixModel}}.
#' @param neg_ctrl Character. The negative control probe to use as the noise baseline.
#'   The threshold is calculated as \code{mean(neg_ctrl) + n_sd * SD(neg_ctrl)}.
#' @param n_sd Numeric. The number of standard deviations above the mean to set the threshold (default: 1).
#'
#' @return A subsetted \code{NanoStringGeoMxSet} object containing only the detected proteins.
#'   It also prints the lists of kept and excluded proteins to the console.
#'
#' @details
#' \strong{Filtering Logic:}
#' \enumerate{
#'   \item Retrieves the \code{Best_Model_Summary} table from \code{experimentData(geomx_set)@other}.
#'   \item Calculates the background threshold using the raw data of the selected \code{neg_ctrl}.
#'   \item For each protein, retrieves the \code{Best_Signal_Mean} (from the optimal model).
#'   \item If \code{Best_Signal_Mean > Threshold}, the protein is retained.
#'   \item Negative control probes are always retained for reference.
#' }
#'
#' @importFrom Biobase assayDataElement experimentData
#' @importFrom dplyr filter
#' @export
#'
#' @seealso \code{\link{BestMixModel}}, \code{\link{PlotMixModel}}
#'
#' @examples
#' \dontrun{
#'   # 1. Run the tournament first
#'   geomx_set <- BestMixModel(geomx_set, ncomps = 3)
#'
#'   # 2. Filter using the best models and Rt IgG2a background
#'   geomx_filtered <- FilterProteins(geomx_set, neg_ctrl = "Rt IgG2a")
#' }
FilterProteins <- function(geomx_set, neg_ctrl = "Rt IgG2a", n_sd = 1) {
  
  require(dplyr)
  require(Biobase)
  
  # --- 1. Validation Checks ---
  # Retrieve the summary from the new experimentData slot
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
  
  # Use the raw data of the negative control to set the baseline
  neg_vals <- as.numeric(log_mat[neg_ctrl, ])
  neg_vals <- neg_vals[is.finite(neg_vals)]
  
  threshold <- mean(neg_vals, na.rm = TRUE) + (n_sd * sd(neg_vals, na.rm = TRUE))
  
  message("Filtering Threshold (", neg_ctrl, " + ", n_sd, "SD): ", round(threshold, 3))
  
  # --- 3. Identify Proteins to Keep ---
  all_proteins <- rownames(log_mat)
  kept_proteins <- c()
  
  # Always keep the negative controls themselves for QC plots
  common_negs <- c("Hmr IgG", "Ms IgG2b", "Rb IgG", "Rt IgG2a", "Ms IgG1")
  # Also keep the specific user-selected control if not in the common list
  common_negs <- unique(c(common_negs, neg_ctrl))
  
  kept_proteins <- intersect(all_proteins, common_negs)
  
  # Iterate through the Best Model Summary table
  # This uses the strictly selected model parameters
  for (i in seq_len(nrow(best_model_df))) {
    prot_name <- best_model_df$Protein[i]
    signal_mean <- best_model_df$Best_Signal_Mean[i]
    
    # Skip if signal is NA (e.g., model failed to converge)
    if (is.na(signal_mean)) next
    
    # STRICT DECISION:
    if (signal_mean > threshold) {
      kept_proteins <- c(kept_proteins, prot_name)
    }
  }
  
  # Unique list of survivors
  kept_proteins <- unique(kept_proteins)
  
  # Identify Excluded
  excluded_proteins <- setdiff(all_proteins, kept_proteins)
  
  # --- 4. Report Results ---
  n_total <- length(all_proteins)
  n_kept  <- length(kept_proteins)
  n_excl  <- length(excluded_proteins)
  
  message("\n==========================================")
  message(" FILTERING SUMMARY (Strict Bayes Logic)")
  message("==========================================")
  message("Total Proteins Input: ", n_total)
  message("Proteins Retained:    ", n_kept)
  message("Proteins Excluded:    ", n_excl)
  
  if (n_excl > 0) {
    message("\n--- EXCLUDED PROTEINS (Best Model Signal < Threshold) ---")
    cat(paste(excluded_proteins, collapse = ", "), "\n")
  } else {
    message("\n--- NO PROTEINS EXCLUDED ---")
  }
  
  if (n_kept > 0) {
    message("\n--- RETAINED PROTEINS ---")
    cat(paste(kept_proteins, collapse = ", "), "\n")
  }
  
  message("==========================================\n")
  
  # --- 5. Execute Subset ---
  geomx_subset <- geomx_set[kept_proteins, ]
  
  return(geomx_subset)
}