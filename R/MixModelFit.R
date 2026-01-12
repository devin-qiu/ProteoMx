#' Fit Mixture Models (Batch Engine)
#'
#' Fits Gaussian mixture models for all proteins in a GeoMxSet object.
#' Results are stored directly inside the object's metadata.
#'
#' @param geomx_set A \code{NanoStringGeoMxSet} object. Must be Q3 normalized first.
#' @param ncomps Integer. The maximum number of components to fit (default: 3).
#' @param neg_ctrl Character. The name of the negative control probe (default: "Rt IgG2a").
#'   Used to calculate the detection threshold (Mean + 1 SD).
#'
#' @return The input \code{geomx_set} object, updated with fit results stored in 
#'   \code{experimentData(object)@other$MixModel}.
#'
#' @details
#' This function iterates through every protein and fits mixture models for every combination of:
#' \itemize{
#'   \item Components: 1 to \code{ncomps}
#'   \item Variance: Equal (ev=TRUE) and Unequal (ev=FALSE)
#' }
#' 
#' It uses \code{do.call} to robustly handle scoping issues in the underlying \code{mixR} package.
#'
#' @importFrom mixR mixfit
#' @importFrom Biobase experimentData experimentData<- assayDataElement
#' @importFrom dplyr tibble
#' @export
#'
#' @seealso \code{\link{Q3Normalize}}, \code{\link{PlotMixModel}}
#'
#' @examples
#' \dontrun{
#'   # Fit up to 3 components using 'Rt IgG2a' as background
#'   geomx_set <- MixModelFit(geomx_set, ncomps = 3, neg_ctrl = "Rt IgG2a")
#' }
MixModelFit <- function(geomx_set, ncomps = 3, neg_ctrl = "Rt IgG2a") {
  
  require(mixR)
  require(dplyr)
  require(tibble)
  require(Biobase)
  
  # --- 1. Robust Data Access ---
  if (!"q_norm" %in% names(assayData(geomx_set))) {
    stop("Error: 'q_norm' data is NULL. Please run Q3Normalize() first.")
  }
  
  q3_mat <- assayDataElement(geomx_set, "q_norm")
  
  # Dimension check
  if (is.null(q3_mat) || nrow(q3_mat) == 0 || ncol(q3_mat) == 0) {
    stop("Error: 'q_norm' matrix is empty.")
  }
  
  # Clean columns
  bad_cols <- apply(q3_mat, 2, function(x) all(is.na(x) | is.infinite(x)))
  if (sum(bad_cols) > 0) q3_mat <- q3_mat[, !bad_cols]
  
  log_mat <- log10(q3_mat + 1)
  
  if (!neg_ctrl %in% rownames(log_mat)) {
    stop("Error: Negative control '", neg_ctrl, "' not found.")
  }
  
  # --- 2. Threshold ---
  neg_vals <- as.numeric(log_mat[neg_ctrl, ])
  neg_vals <- neg_vals[is.finite(neg_vals)]
  threshold <- mean(neg_vals) + 1 * sd(neg_vals)
  message("Threshold calculated: ", round(threshold, 3))
  
  # --- 3. Internal Helper (THE FIX IS HERE) ---
  fit_single_protein <- function(vals, n, ev) {
    # Ensure strict vector format
    vals <- as.vector(vals)
    vals <- vals[is.finite(vals)]
    
    if (length(vals) < 10) return(NULL)
    
    # WE USE do.call TO FORCE EVALUATION OF 'vals'
    # This prevents mixfit from looking for the variable name 'vals'
    args_list <- list(
      x = vals, 
      ncomp = n, 
      family = "norm", 
      ev = ev, 
      init.method = "kmeans"
    )
    
    fit_result <- tryCatch(
      do.call(mixR::mixfit, args_list),
      error = function(e) NULL
    )
    
    if (is.null(fit_result)) return(NULL)
    
    # Shift means
    mix_mean <- sum(fit_result$pi * fit_result$mu)
    fit_result$mu <- fit_result$mu + (mean(vals) - mix_mean)
    return(fit_result)
  }
  
  # --- 4. Main Batch Loop ---
  all_proteins <- rownames(log_mat)
  ev_values <- c(FALSE, TRUE)
  fit_results_all <- list()
  
  total_steps <- length(ev_values) * length(1:ncomps)
  step_count <- 0
  
  message("Starting Batch Fit for ", length(all_proteins), " proteins...")
  
  for (ev_flag in ev_values) {
    for (n in 1:ncomps) {
      step_count <- step_count + 1
      message(paste0("[", step_count, "/", total_steps, "] Fitting ev=", ev_flag, ", ncomp=", n, "..."))
      
      protein_fits <- vector("list", length(all_proteins))
      keep_flags   <- logical(length(all_proteins))
      
      pb <- txtProgressBar(min = 0, max = length(all_proteins), style = 3)
      
      for (i in seq_along(all_proteins)) {
        setTxtProgressBar(pb, i)
        vals <- as.numeric(log_mat[i, ]) 
        
        fit <- fit_single_protein(vals, n, ev_flag)
        
        if (!is.null(fit)) {
          protein_fits[[i]] <- fit
          keep_flags[i] <- max(fit$mu, na.rm = TRUE) > threshold
        } else {
          keep_flags[i] <- FALSE
        }
      }
      close(pb)
      
      fits_tbl <- tibble(
        Protein    = all_proteins,
        fit_result = protein_fits,
        keep       = keep_flags
      )
      
      key <- paste0("ev", if (ev_flag) "T" else "F", "_ncomp_", n)
      fit_results_all[[key]] <- fits_tbl
    }
  }
  
  # --- 5. Store Results ---
  res_package <- list(
    fits = fit_results_all,
    threshold = threshold,
    params = list(ncomps = ncomps, neg_ctrl = neg_ctrl),
    timestamp = Sys.time()
  )
  
  ed <- experimentData(geomx_set)
  ed@other$MixModel <- res_package
  experimentData(geomx_set) <- ed
  
  message("Done! Fits stored in experimentData(object)@other$MixModel")
  return(geomx_set)
}