#' Filter Proteins by Mixture Model Detection
#'
#' Subsets the GeoMxSet object to retain only proteins that have a detectable signal 
#' above the negative control background in at least one fitted model.
#'
#' @param geomx_set A \code{NanoStringGeoMxSet} object processed by \code{\link{MixModelFit}}.
#' @param neg_ctrl Character. The negative control probe to use as the noise baseline.
#'   The threshold is calculated as \code{mean(neg_ctrl) + 1 * SD(neg_ctrl)}.
#'   
#'   \strong{Recommended Options:}
#'   \itemize{
#'     \item "Rt IgG2a" (Often the lowest background)
#'     \item "Hmr IgG"
#'     \item "Ms IgG2b"
#'     \item "Rb IgG"
#'     \item "Ms IgG1"
#'   }
#'
#' @return A subsetted \code{NanoStringGeoMxSet} object containing only the detected proteins.
#'   It also prints the lists of kept and excluded proteins to the console.
#'
#' @details
#' \strong{Filtering Logic:}
#' \enumerate{
#'   \item Recalculates the background threshold based on the user-selected \code{neg_ctrl}.
#'   \item Iterates through every protein and every model fit stored in \code{experimentData}.
#'   \item Checks if the \strong{highest component mean} (max mu) of the model exceeds the threshold.
#'   \item If a protein's signal is \strong{below} the threshold in \strong{ALL} tested models 
#'         (ncomp 1..N, ev T/F), it is considered "undetected" and removed.
#'   \item Negative control probes are always retained for reference.
#' }
#'
#' @importFrom Biobase assayDataElement experimentData
#' @importFrom dplyr filter
#' @export
#'
#' @examples
#' \dontrun{
#'   # Filter using the standard Rt IgG2a background
#'   geomx_filtered <- FilterProteins(geomx_set, neg_ctrl = "Rt IgG2a")
#' }
FilterProteins <- function(geomx_set, neg_ctrl = "Rt IgG2a") {
  
  require(dplyr)
  require(Biobase)
  
  # --- 1. Validation Checks ---
  mix_res <- experimentData(geomx_set)@other$MixModel
  if (is.null(mix_res)) stop("No MixModel results found. Please run MixModelFit() first.")
  
  if (!"q_norm" %in% names(assayData(geomx_set))) stop("q_norm data missing. Run Q3Normalize() first.")
  q3_mat <- assayDataElement(geomx_set, "q_norm")
  
  # --- 2. Calculate Dynamic Threshold ---
  log_mat <- log10(q3_mat + 1)
  
  if (!neg_ctrl %in% rownames(log_mat)) {
    stop("Negative control '", neg_ctrl, "' not found in dataset.")
  }
  
  neg_vals <- as.numeric(log_mat[neg_ctrl, ])
  neg_vals <- neg_vals[is.finite(neg_vals)]
  
  threshold <- mean(neg_vals) + 1 * sd(neg_vals)
  
  message("Filtering Threshold (", neg_ctrl, " + 1SD): ", round(threshold, 3))
  
  # --- 3. Identify Proteins to Keep ---
  all_proteins <- rownames(log_mat)
  kept_proteins <- c()
  
  # Always keep the negative controls themselves
  common_negs <- c("Hmr IgG", "Ms IgG2b", "Rb IgG", "Rt IgG2a", "Ms IgG1")
  kept_proteins <- c(kept_proteins, intersect(all_proteins, common_negs))
  
  # Iterate through every model bucket
  for (key in names(mix_res$fits)) {
    
    fit_tbl <- mix_res$fits[[key]]
    
    for (i in seq_len(nrow(fit_tbl))) {
      p_name <- fit_tbl$Protein[i]
      fit_obj <- fit_tbl$fit_result[[i]]
      
      if (!is.null(fit_obj)) {
        # If highest component > threshold, keep it
        if (max(fit_obj$mu, na.rm = TRUE) > threshold) {
          kept_proteins <- c(kept_proteins, p_name)
        }
      }
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
  message(" FILTERING SUMMARY")
  message("==========================================")
  message("Total Proteins Input: ", n_total)
  message("Proteins Retained:    ", n_kept)
  message("Proteins Excluded:    ", n_excl)
  
  if (n_excl > 0) {
    message("\n--- EXCLUDED PROTEINS (Signal < Threshold) ---")
    cat(paste(excluded_proteins, collapse = ", "), "\n")
  } else {
    message("\n--- NO PROTEINS EXCLUDED ---")
  }
  
  if (n_kept > 0) {
    message("\n--- RETAINED PROTEINS ---")
    # Print FULL list without truncation
    cat(paste(kept_proteins, collapse = ", "), "\n")
  }
  
  message("==========================================\n")
  
  # --- 5. Execute Subset ---
  geomx_subset <- geomx_set[kept_proteins, ]
  
  return(geomx_subset)
}