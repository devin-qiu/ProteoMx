#' Select Best Mixture Model (Manual Tournament + Selection Log)
#'
#' Robustly identifies the optimal mixture model for each protein by manually fitting
#' all combinations of components and variance structures.
#'
#' @description
#' This function performs a "Manual Tournament" to find the best model (lowest BIC).
#' Unlike \code{mixR::select}, this function handles errors gracefully: if a specific
#' model fails to converge, it is simply skipped, preventing the entire run from crashing.
#'
#' It returns a main summary of the mathematical winners, plus a detailed
#' \code{Selection_Log} for every protein to allow for manual inspection (e.g., checking the "Elbow Rule").
#'
#' @param geomx_set A \code{NanoStringGeoMxSet} object. Must be processed by \code{\link{Q3Normalize}} first.
#' @param ncomps Integer. The maximum number of components to test (default: 3).
#'   The function will test all values from 1 to \code{ncomps}.
#'
#' @return A \code{\link[tibble]{tibble}} with one row per protein and the following columns:
#' \itemize{
#'   \item \code{Protein}: The name of the protein.
#'   \item \code{Best_NComp}: The number of components in the model with the lowest BIC.
#'   \item \code{Best_EV}: The variance assumption (TRUE/FALSE) of the winning model.
#'   \item \code{Best_BIC}: The BIC score of the winning model.
#'   \item \code{Selection_Log}: A list-column. Each entry contains a nested tibble with the
#'         history of all successful fits for that protein, sorted by model simplicity
#'         (fewer components first, then equal variance).
#' }
#'
#' @details
#' \strong{The Tournament Logic:}
#' For each protein, the function iterates through:
#' \itemize{
#'   \item \code{ncomp}: 1 to \code{ncomps}
#'   \item \code{ev}: \code{FALSE} (unequal variance) and \code{TRUE} (equal variance)
#' }
#' The function uses \code{\link{do.call}} to robustly pass data to \code{mixR::mixfit},
#' bypassing common scoping errors found in the native package.
#'
#' \strong{Selection Log Sorting:}
#' The inner tables in \code{Selection_Log} are sorted by \strong{Parsimony} (simplicity),
#' not by BIC.
#' \enumerate{
#'   \item Lowest \code{ncomp} (1) to highest.
#'   \item For the same \code{ncomp}, \code{ev=TRUE} (simpler) comes before \code{ev=FALSE}.
#' }
#' This allows users to easily scan the log and apply the "Elbow Rule" (choosing a simpler model
#' if the BIC improvement of a more complex one is negligible).
#'
#' @importFrom mixR mixfit
#' @importFrom dplyr bind_rows arrange desc
#' @importFrom tibble tibble
#' @importFrom Biobase assayDataElement
#' @export
#'
#' @seealso \code{\link{MixModelFit}}, \code{\link{PlotMixModel}}
#'
#' @examples
#' \dontrun{
#'   # 1. Run the tournament (computationally intensive)
#'   best_res <- BestMixModel(geomx_set, ncomps = 3)
#'
#'   # 2. View the winners
#'   print(head(best_res))
#'
#'   # 3. Inspect the detailed log for the first protein to check for overfitting
#'   # (Look for marginal BIC improvements between ncomp=2 and ncomp=3)
#'   print(best_res$Selection_Log[[1]])
#' }
BestMixModel <- function(geomx_set, ncomps = 3) {
  
  require(mixR)
  require(dplyr)
  require(tibble)
  require(Biobase)
  
  # --- 1. Data Prep ---
  if (!"q_norm" %in% names(assayData(geomx_set))) stop("Error: 'q_norm' missing. Run Q3Normalize() first.")
  q3_mat <- assayDataElement(geomx_set, "q_norm")
  
  # Clean columns
  bad_cols <- apply(q3_mat, 2, function(x) all(is.na(x) | is.infinite(x)))
  if (sum(bad_cols) > 0) q3_mat <- q3_mat[, !bad_cols]
  
  log_mat <- log10(q3_mat + 1)
  all_proteins <- rownames(log_mat)
  
  # Output containers
  best_ncomp_sel <- rep(NA_integer_, length(all_proteins))
  best_ev_sel    <- rep(NA, length(all_proteins)) 
  best_bic_sel   <- rep(NA_real_, length(all_proteins))
  selection_logs <- vector("list", length(all_proteins))
  
  message("Running Manual Tournament (1 to ", ncomps, " components)...")
  pb <- txtProgressBar(min = 0, max = length(all_proteins), style = 3)
  
  # --- 2. Loop Over Proteins ---
  for (i in seq_along(all_proteins)) {
    setTxtProgressBar(pb, i)
    
    vals <- as.vector(as.numeric(log_mat[i, ]))
    vals <- vals[is.finite(vals)]
    
    if (length(vals) < 10) next 
    
    # Storage for this protein's tournament
    protein_tournament <- list()
    
    # Variables to track the mathematical winner (Lowest BIC)
    best_bic_local <- Inf
    best_model_local <- list(ncomp = NA, ev = NA)
    
    # --- 3. The Tournament Loop ---
    for (ev_flag in c(FALSE, TRUE)) {
      for (k in 1:ncomps) {
        
        # Run fit safely using do.call
        fit <- tryCatch(
          do.call(mixR::mixfit, list(x = vals, ncomp = k, family = "norm", ev = ev_flag, init.method = "kmeans")),
          error = function(e) NULL
        )
        
        # Record Result
        bic_val <- NA_real_
        if (!is.null(fit) && !is.null(fit$bic) && !is.nan(fit$bic)) {
          bic_val <- fit$bic
          
          # Check if this is the new mathematical winner
          if (bic_val < best_bic_local) {
            best_bic_local <- bic_val
            best_model_local <- list(ncomp = k, ev = ev_flag)
          }
        }
        
        protein_tournament[[length(protein_tournament) + 1]] <- tibble(
          ncomp = k,
          ev = ev_flag,
          bic = bic_val
        )
      }
    }
    
    # --- 4. Compile Results for this Protein ---
    
    # Store the mathematical winner in the main columns
    if (best_bic_local < Inf) {
      best_ncomp_sel[i] <- best_model_local$ncomp
      best_ev_sel[i]    <- best_model_local$ev
      best_bic_sel[i]   <- best_bic_local
    }
    
    # Store the Full Log, sorted by SIMPLICITY (Parsimony)
    # Order: Lower ncomp first. For ties, Equal Variance (TRUE) first.
    log_tbl <- bind_rows(protein_tournament) %>%
      filter(!is.na(bic)) %>%
      arrange(ncomp, desc(ev)) # desc(TRUE) = 1 puts TRUE before FALSE (0)
    
    selection_logs[[i]] <- log_tbl
  }
  close(pb)
  
  # --- 5. Format Final Output ---
  tibble(
    Protein       = all_proteins,
    Best_NComp    = best_ncomp_sel,
    Best_EV       = best_ev_sel,
    Best_BIC      = best_bic_sel,
    Selection_Log = selection_logs
  )
}
