#' Check Ambiguous Mixture Models
#'
#' Identifies proteins flagged with an "Ambiguity_Note" during the BestMixModel 
#' tournament and extracts the parameters for both the chosen model and 
#' the alternative competing model for easy manual review.
#'
#' @description
#' This function should be run *after* \code{\link{FilterProteins}}. It cross-references
#' the flagged proteins with the currently surviving proteins. It then extracts the 
#' exact model parameters (NComp, EV, BIC) for both the chosen model and the 
#' best alternative model, saving this comparison table back into the object.
#'
#' @param geomx_set The subsetted \code{NanoStringGeoMxSet} object output by \code{FilterProteins}.
#'
#' @return The updated \code{NanoStringGeoMxSet} object. A review table is saved to 
#'   \code{experimentData(geomx_set)@other$Ambiguous_Models_Review}.
#'
#' @importFrom Biobase experimentData
#' @importFrom dplyr filter arrange slice bind_rows
#' @importFrom tibble tibble
#' @export
CheckBestModel <- function(geomx_set) {
  
  require(dplyr)
  require(tibble)
  require(Biobase)
  
  # --- 1. Validation Checks ---
  exp_data <- experimentData(geomx_set)
  sum_df <- exp_data@other$Best_Model_Summary
  
  if (is.null(sum_df)) {
    stop("Best_Model_Summary not found. Please ensure the object was processed by BestMixModel().")
  }
  
  # --- 2. Extract Ambiguous Proteins ---
  ambig_df <- sum_df %>% filter(Ambiguity_Note != "")
  ambig_proteins <- ambig_df$Protein
  n_ambig <- length(ambig_proteins)
  
  if (n_ambig == 0) {
    message("Great news! No proteins were flagged with ambiguous model fits.")
    exp_data@other$Ambiguous_Models_Review <- tibble()
    experimentData(geomx_set) <- exp_data
    return(geomx_set)
  }
  
  # --- 3. Cross-Reference & Extract Alternative Parameters ---
  survivors <- rownames(geomx_set)
  ambig_retained <- intersect(ambig_proteins, survivors)
  ambig_excluded <- setdiff(ambig_proteins, survivors)
  
  # Build a clean review table by extracting the BEST ALTERNATIVE from the log
  review_list <- lapply(seq_len(nrow(ambig_df)), function(i) {
    prot <- ambig_df$Protein[i]
    status <- ifelse(prot %in% survivors, "Retained", "Excluded")
    
    # The chosen model parameters
    c_ncomp <- ambig_df$Chosen_NComp[i]
    c_ev <- ambig_df$Chosen_EV[i]
    c_bic <- ambig_df$Chosen_BIC[i]
    
    # Dig into the Selection_Log to find the true alternative
    log_tbl <- ambig_df$Selection_Log[[i]]
    
    # FIX: Remove the chosen model from the log, then find the one with the lowest BIC
    alt_row <- log_tbl %>% 
      filter(!(ncomp == c_ncomp & ev == c_ev)) %>% 
      arrange(bic) %>% 
      slice(1)
    
    tibble(
      Protein = prot,
      Status = status,
      Chosen_NComp = c_ncomp,
      Chosen_EV = c_ev,
      Chosen_BIC = round(c_bic, 2),
      Alt_NComp = alt_row$ncomp,
      Alt_EV = alt_row$ev,
      Alt_BIC = round(alt_row$bic, 2)
    )
  })
  
  review_df <- bind_rows(review_list)
  
  # Store the review table into the filtered object
  exp_data@other$Ambiguous_Models_Review <- review_df
  experimentData(geomx_set) <- exp_data
  
  # --- 4. Print Summary ---
  message("\n==========================================")
  message(" AMBIGUOUS MODEL FIT SUMMARY")
  message("==========================================")
  message("Total Ambiguous Proteins: ", n_ambig)
  message("------------------------------------------")
  
  message("-> RETAINED (Passed filter despite ambiguity): ", length(ambig_retained))
  if (length(ambig_retained) > 0) {
    cat(paste(ambig_retained, collapse = ", "), "\n")
  }
  
  message("-> EXCLUDED (Failed filter, check if alternative model would have passed): ", length(ambig_excluded))
  if (length(ambig_excluded) > 0) {
    cat(paste(ambig_excluded, collapse = ", "), "\n")
  }
  message("==========================================\n")
  
  # --- 5. Return Updated Object ---
  return(geomx_set)
}