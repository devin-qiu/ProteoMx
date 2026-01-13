#' Plot Mixture Model Fit
#'
#' Visualizes the fitted Gaussian mixture model for a specific protein, overlaid with a
#' dynamic negative control threshold to assess signal detection.
#'
#' @description
#' This function generates a density plot of the log10-transformed data for a selected protein.
#' It overlays the fitted mixture components (as filled areas) and the total model fit (as a black line).
#' A red dashed line indicates the background noise threshold, calculated dynamically based on
#' the selected negative control probe.
#'
#' @param geomx_set A \code{NanoStringGeoMxSet} object processed by \code{\link{MixModelFit}}.
#' @param protein Character. The name of the protein to plot.
#' @param ncomp Integer (Optional). The number of components to plot.
#'   If \code{NULL}, the function auto-selects the best model based on the lowest BIC score
#'   stored in the object.
#' @param ev Logical (Optional). Equal variance assumption.
#'   If \code{NULL}, the function auto-selects the best model based on BIC.
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
#' @return A \code{\link[ggplot2]{ggplot}} object showing:
#' \itemize{
#'   \item Histogram of observed data (density scale)
#'   \item Gaussian curves for individual components (filled areas)
#'   \item Total mixture curve (black line)
#'   \item Threshold line (red dashed) with annotation
#'   \item Subtitle indicating if the protein is considered "Detected" (max component mean > threshold)
#' }
#'
#' @details
#' The function retrieves fit parameters directly from \code{experimentData(object)@other$MixModel}.
#' It calculates the detection threshold on the fly using the provided \code{neg_ctrl}.
#' If \code{ncomp} and \code{ev} are not provided, it scans the stored results to find the
#' model combination (ncomp/ev) with the lowest BIC for that specific protein.
#'
#' @import ggplot2
#' @importFrom Biobase experimentData assayDataElement
#' @importFrom dplyr filter
#' @importFrom tibble tibble
#' @export
#'
#' @seealso \code{\link{MixModelFit}}, \code{\link{FilterProteins}}
#'
#' @examples
#' \dontrun{
#'   # Plot using auto-detected best parameters and default control
#'   PlotMixModel(geomx_set, "CD44")
#'
#'   # Plot using a specific control and model parameters
#'   PlotMixModel(geomx_set, "CD44", ncomp = 2, ev = FALSE, neg_ctrl = "Ms IgG1")
#' }
PlotMixModel <- function(geomx_set, protein, ncomp = NULL, ev = NULL, neg_ctrl = "Rt IgG2a") {
  
  require(ggplot2)
  require(dplyr)
  require(tibble)
  require(Biobase)
  
  # --- 1. Retrieve Analysis Results ---
  mix_res <- experimentData(geomx_set)@other$MixModel
  
  if (is.null(mix_res)) {
    stop("No mixture model results found. Please run: geomx_set <- MixModelFit(geomx_set) first.")
  }
  
  # --- 2. Calculate Threshold Dynamically ---
  # We access the data fresh to ensure we match the current object state
  q3_mat <- assayDataElement(geomx_set, "q_norm")
  log_mat <- log10(q3_mat + 1)
  
  if (!neg_ctrl %in% rownames(log_mat)) {
    stop("Negative control '", neg_ctrl, "' not found.")
  }
  
  neg_vals <- as.numeric(log_mat[neg_ctrl, ])
  neg_vals <- neg_vals[is.finite(neg_vals)]
  threshold <- mean(neg_vals) + 1 * sd(neg_vals)
  
  # --- 3. Auto-Select Best Model (if params are NULL) ---
  if (is.null(ncomp) || is.null(ev)) {
    best_bic <- Inf
    best_key <- NULL
    
    for (key in names(mix_res$fits)) {
      fit_row <- mix_res$fits[[key]] %>% filter(Protein == !!protein)
      if (nrow(fit_row) > 0) {
        fit_obj <- fit_row$fit_result[[1]]
        if (!is.null(fit_obj) && !is.na(fit_obj$bic)) {
          if (fit_obj$bic < best_bic) {
            best_bic <- fit_obj$bic
            best_key <- key
          }
        }
      }
    }
    
    if (is.null(best_key)) stop("Could not find a valid model fit for protein: ", protein)
    
    ev_str <- sub("ev([TF])_.*", "\\1", best_key)
    ev <- (ev_str == "T")
    ncomp <- as.integer(sub(".*_ncomp_", "", best_key))
    key <- best_key
  } else {
    key <- paste0("ev", if(ev) "T" else "F", "_ncomp_", ncomp)
  }
  
  # --- 4. Extract Specific Fit Data ---
  target_tbl <- mix_res$fits[[key]]
  if (is.null(target_tbl)) stop("No data found for parameter set: ", key)
  
  row_data <- target_tbl %>% filter(Protein == !!protein)
  if (nrow(row_data) == 0) stop("Protein not found in results.")
  
  fit <- row_data$fit_result[[1]]
  if (is.null(fit)) stop("Model fit is NULL for this protein/parameter combination.")
  
  # --- 5. Determine Detection Status ---
  # Logic: Detected if the highest component mean > threshold
  is_detected <- max(fit$mu, na.rm = TRUE) > threshold
  
  # --- 6. Prepare Plot Data ---
  vals <- as.numeric(log_mat[protein, ])
  vals <- vals[is.finite(vals)]
  
  x_vals <- seq(min(vals), max(vals), length.out = 1000)
  
  # Densities
  comp_mat <- vapply(seq_along(fit$mu), function(i) {
    fit$pi[i] * dnorm(x_vals, mean = fit$mu[i], sd = fit$sd[i])
  }, numeric(length(x_vals)))
  
  mix_y <- rowSums(comp_mat)
  mix_df  <- tibble(x = x_vals, density = mix_y)
  
  comp_df <- bind_rows(lapply(seq_along(fit$mu), function(i) {
    tibble(x = x_vals, density = comp_mat[, i], comp = paste0("Comp ", i))
  }))
  
  # --- 7. Plotting ---
  ggplot(tibble(val = vals), aes(x = val)) +
    geom_histogram(aes(y = after_stat(density)), bins = 40, fill = "grey90", color = "grey60") +
    geom_area(data = comp_df, aes(x = x, y = density, fill = comp), alpha = 0.4, position = "identity") +
    geom_line(data = mix_df, aes(x = x, y = density), size = 1, color = "black") +
    
    # Dynamic Threshold Line
    geom_vline(xintercept = threshold, linetype = "dashed", color = "red", size = 0.8) +
    annotate("text", x = threshold, y = Inf, 
             label = paste("Threshold:", round(threshold, 2), "(", neg_ctrl, ")"), 
             color = "red", angle = 90, vjust = -0.5, hjust = 1.1) +
    
    labs(
      title = paste0(protein),
      subtitle = paste0("Model: ncomp=", ncomp, ", ev=", ev, " | Detected? ", if(is_detected) "YES" else "NO"),
      x = "log10(Q3 Normalized Counts)", 
      y = "Density"
    ) +
    scale_fill_brewer(palette = "Set2") + 
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11, color = if(is_detected) "darkgreen" else "darkred"),
      legend.title = element_blank(),
      legend.position = "bottom"
    )
}