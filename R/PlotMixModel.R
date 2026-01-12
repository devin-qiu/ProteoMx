#' Plot Mixture Model Fit
#'
#' Visualizes the fitted mixture model for a specific protein. 
#'
#' @param geomx_set A \code{NanoStringGeoMxSet} object processed by \code{MixModelFit}.
#' @param protein Character. The name of the protein to plot.
#' @param ncomp Integer (Optional). Number of components to plot. If \code{NULL}, 
#'   the function auto-selects the best model based on the lowest BIC stored in the object.
#' @param ev Logical (Optional). Equal variance assumption. If \code{NULL}, 
#'   auto-selects the best model.
#'
#' @return A \code{ggplot} object showing:
#' \itemize{
#'   \item Histogram of observed data (density scale)
#'   \item Gaussian curves for individual components (filled areas)
#'   \item Total mixture curve (black line)
#'   \item Threshold line (red dashed)
#' }
#'
#' @details
#' This function retrieves fit parameters directly from \code{experimentData(object)@other$MixModel}.
#' It uses density-on-density plotting to ensure the curves align perfectly with the histogram.
#'
#' @import ggplot2
#' @importFrom Biobase experimentData assayDataElement
#' @export
#'
#' @examples
#' \dontrun{
#'   # Plot using auto-detected best parameters
#'   PlotMixModel(geomx_set, "CD44")
#'
#'   # Plot a specific model (e.g., 2 components, unequal variance)
#'   PlotMixModel(geomx_set, "CD44", ncomp = 2, ev = FALSE)
#' }
PlotMixModel <- function(geomx_set, protein, ncomp = NULL, ev = NULL) {
  
  require(ggplot2)
  require(dplyr)
  require(tibble)
  require(Biobase)
  
  # --- 1. Retrieve Analysis Results ---
  mix_res <- experimentData(geomx_set)@other$MixModel
  
  if (is.null(mix_res)) {
    stop("No mixture model results found. Please run: geomx_set <- MixModelFit(geomx_set) first.")
  }
  
  # --- 2. Auto-Select Best Model (if params are NULL) ---
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
  
  # --- 3. Extract Specific Fit Data ---
  target_tbl <- mix_res$fits[[key]]
  if (is.null(target_tbl)) stop("No data found for parameter set: ", key)
  
  row_data <- target_tbl %>% filter(Protein == !!protein)
  if (nrow(row_data) == 0) stop("Protein not found in results.")
  
  fit <- row_data$fit_result[[1]]
  is_kept <- row_data$keep
  threshold <- mix_res$threshold
  
  if (is.null(fit)) stop("Model fit is NULL for this protein/parameter combination.")
  
  # --- 4. Prepare Data ---
  # Re-extract raw data to ensure valid plot range
  q3_mat <- assayDataElement(geomx_set, "q_norm")
  vals <- as.numeric(log10(q3_mat[protein, ] + 1))
  vals <- vals[is.finite(vals)]
  
  # --- 5. Generate Curves (Density Scale) ---
  # Create a smooth grid for the curves
  x_vals <- seq(min(vals), max(vals), length.out = 1000)
  
  # Calculate densities. dnorm outputs DENSITY, which matches geom_histogram(after_stat(density))
  comp_mat <- vapply(seq_along(fit$mu), function(i) {
    fit$pi[i] * dnorm(x_vals, mean = fit$mu[i], sd = fit$sd[i])
  }, numeric(length(x_vals)))
  
  # Total mixture is the sum of components
  mix_y <- rowSums(comp_mat)
  
  mix_df  <- tibble(x = x_vals, density = mix_y)
  
  # Component DF for filled areas
  comp_df <- bind_rows(lapply(seq_along(fit$mu), function(i) {
    tibble(x = x_vals, density = comp_mat[, i], comp = paste0("Comp ", i))
  }))
  
  # --- 6. Plotting ---
  ggplot(tibble(val = vals), aes(x = val)) +
    # Histogram: Normalized to Density (Area=1)
    geom_histogram(aes(y = after_stat(density)), bins = 40, fill = "grey90", color = "grey60") +
    
    # Filled areas for components
    geom_area(data = comp_df, aes(x = x, y = density, fill = comp), alpha = 0.4, position = "identity") +
    
    # Black line for total fit
    geom_line(data = mix_df, aes(x = x, y = density), size = 1, color = "black") +
    
    # Threshold Line
    geom_vline(xintercept = threshold, linetype = "dashed", color = "red", size = 0.8) +
    
    # Annotations
    annotate("text", x = threshold, y = Inf, 
             label = paste("Threshold:", round(threshold, 2)), 
             color = "red", angle = 90, vjust = -0.5, hjust = 1.1) +
    
    # Styling
    labs(
      title = paste0(protein),
      subtitle = paste0("Model: ncomp=", ncomp, ", ev=", ev, " | Detected? ", if(is_kept) "YES" else "NO"),
      x = "log10(Q3 Normalized Counts)", 
      y = "Density"
    ) +
    scale_fill_brewer(palette = "Set2") + 
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11, color = if(is_kept) "darkgreen" else "darkred"),
      legend.title = element_blank(),
      legend.position = "bottom"
    )
}