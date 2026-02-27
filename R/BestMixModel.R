#' Select Best Mixture Model (Bayes Factor + Occam's Razor)
#'
#' Robustly identifies the optimal mixture model for each protein by manually fitting
#' all combinations of components and variance structures, and using Bayes Factors
#' to penalize unnecessary complexity (Occam's Razor).
#'
#' @description
#' This function performs a "Manual Tournament" to find the best model.
#' It converts BIC scores to Bayes Factors using \code{bayestestR}. If a simpler
#' model (fewer components) is statistically indistinguishable from a more complex model
#' (log10(Bayes Factor) < 0.5), the simpler model is chosen to prevent overfitting.
#'
#' @param geomx_set A \code{NanoStringGeoMxSet} object. Must be processed by \code{\link{Q3Normalize}} first.
#' @param ncomps Integer. The maximum number of components to test (default: 3).
#'
#' @return A \code{\link[tibble]{tibble}} with one row per protein and the following columns:
#' \itemize{
#'   \item \code{Protein}: The name of the protein.
#'   \item \code{Chosen_NComp}: The number of components in the winning model.
#'   \item \code{Chosen_EV}: The variance assumption of the winning model.
#'   \item \code{Chosen_BIC}: The BIC score of the winning model.
#'   \item \code{Best_Signal_Mean}: The highest component mean from the winning model (used for filtering).
#'   \item \code{Ambiguity_Note}: Flag if the top 2 models were statistically indistinguishable.
#'   \item \code{Selection_Log}: Nested tibble containing the history and Bayes Factors of all successful fits.
#' }
#'
#' @importFrom mixR mixfit
#' @importFrom dplyr bind_rows arrange desc mutate filter
#' @importFrom tibble tibble
#' @importFrom Biobase assayDataElement
#' @export
BestMixModel <- function(geomx_set, ncomps = 3) {
  
  require(mixR)
  require(dplyr)
  require(tibble)
  require(Biobase)
  
  if (!requireNamespace("bayestestR", quietly = TRUE)) {
    stop("Package 'bayestestR' is required for Bayes Factor calculations.")
  }
  
  # --- 1. Data Prep ---
  if (!"q_norm" %in% names(assayData(geomx_set))) stop("Error: 'q_norm' missing. Run Q3Normalize() first.")
  q3_mat <- assayDataElement(geomx_set, "q_norm")
  
  bad_cols <- apply(q3_mat, 2, function(x) all(is.na(x) | is.infinite(x)))
  if (sum(bad_cols) > 0) q3_mat <- q3_mat[, !bad_cols]
  
  log_mat <- log10(q3_mat + 1)
  all_proteins <- rownames(log_mat)
  
  # Output containers
  chosen_ncomp_sel <- rep(NA_integer_, length(all_proteins))
  chosen_ev_sel    <- rep(NA, length(all_proteins)) 
  chosen_bic_sel   <- rep(NA_real_, length(all_proteins))
  best_mean_sel    <- rep(NA_real_, length(all_proteins))
  note_sel         <- rep("", length(all_proteins))
  selection_logs   <- vector("list", length(all_proteins))
  
  message("Running Manual Tournament with Bayes Factor Selection (1 to ", ncomps, " components)...")
  pb <- txtProgressBar(min = 0, max = length(all_proteins), style = 3)
  
  # --- 2. Loop Over Proteins ---
  for (i in seq_along(all_proteins)) {
    setTxtProgressBar(pb, i)
    
    vals <- as.vector(as.numeric(log_mat[i, ]))
    vals <- vals[is.finite(vals)]
    
    if (length(vals) < 10) next 
    
    protein_tournament <- list()
    fitted_models <- list() 
    
    # --- 3. The Tournament Loop ---
    for (ev_flag in c(TRUE, FALSE)) { 
      for (k in 1:ncomps) {
        
        fit <- tryCatch(
          do.call(mixR::mixfit, list(x = vals, ncomp = k, family = "norm", ev = ev_flag, init.method = "kmeans")),
          error = function(e) NULL
        )
        
        bic_val <- NA_real_
        max_mean <- NA_real_
        
        if (!is.null(fit) && !is.null(fit$bic) && !is.nan(fit$bic)) {
          bic_val <- fit$bic
          max_mean <- max(fit$mu, na.rm = TRUE)
          
          model_id <- paste0("n", k, "_ev", ev_flag)
          fitted_models[[model_id]] <- fit
        }
        
        protein_tournament[[length(protein_tournament) + 1]] <- tibble(
          model_id = paste0("n", k, "_ev", ev_flag),
          ncomp = k,
          ev = ev_flag,
          bic = bic_val,
          max_mean = max_mean
        )
      }
    }
    
    # --- 4. Bayes Factor Selection Logic ---
    log_tbl <- bind_rows(protein_tournament) %>% 
      filter(!is.na(bic)) %>%
      filter(!(ncomp == 1 & ev == FALSE)) # FIX: Remove redundant 1-component model
    
    if (nrow(log_tbl) > 0) {
      
      # 4a. Find absolute minimum BIC
      min_bic <- min(log_tbl$bic)
      
      # 4b. Calculate Bayes Factors relative to the minimum BIC model
      log_tbl <- log_tbl %>%
        mutate(
          BF = bayestestR::bic_to_bf(bic, denominator = min_bic),
          log10_BF = abs(log10(BF))
        ) %>%
        arrange(ncomp, desc(ev))
      
      # 4c. Occam's Razor: Pick the SIMPLEST model where log10_BF < 0.5
      valid_models <- log_tbl %>% filter(log10_BF < 0.5)
      
      if (nrow(valid_models) > 0) {
        chosen_model <- valid_models[1, ]
        
        chosen_ncomp_sel[i] <- chosen_model$ncomp
        chosen_ev_sel[i]    <- chosen_model$ev
        chosen_bic_sel[i]   <- chosen_model$bic
        best_mean_sel[i]    <- chosen_model$max_mean
        
        # 4d. Check for Ambiguity
        bic_sorted <- log_tbl %>% arrange(bic)
        if (nrow(bic_sorted) > 1 && bic_sorted$log10_BF[2] < 0.5) {
          note_sel[i] <- "Top 2 models have similar fitting score. Please check manually and decide which model works the best."
        }
      }
    }
    
    selection_logs[[i]] <- log_tbl
  }
  close(pb)
  
  # --- 5. Format Final Output ---
  final_res <- tibble(
    Protein          = all_proteins,
    Chosen_NComp     = chosen_ncomp_sel,
    Chosen_EV        = chosen_ev_sel,
    Chosen_BIC       = chosen_bic_sel,
    Best_Signal_Mean = best_mean_sel,
    Ambiguity_Note   = note_sel,
    Selection_Log    = selection_logs
  )
  
  exp_data <- experimentData(geomx_set)
  exp_data@other$Best_Model_Summary <- final_res
  experimentData(geomx_set) <- exp_data
  
  return(geomx_set)
}