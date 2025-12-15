#' Plot Mixture Model for a Specific Protein
#' 
#' @param mix_model_fit The list object returned by MixModelFit().
#' @param protein Required. The name of the protein to plot.
#' @param ncomp Integer. Number of components to visualize.
#' @param ev Logical. Equal variance assumption (TRUE or FALSE).
#' 
#' @return A ggplot2 object.
#' @export
PlotMixModel <- function(mix_model_fit, protein, ncomp, ev) {
  
  require(ggplot2)
  require(dplyr)
  require(tibble)
  
  # --- 1. Validation ---
  if (missing(protein)) stop("Argument 'protein' is required.")
  
  # Check if protein exists in the data
  if (!protein %in% rownames(mix_model_fit$data)) {
    stop("Protein '", protein, "' not found in the dataset.")
  }
  
  # Check if protein has a fit results
  if (is.null(mix_model_fit$fits[[protein]])) {
    stop("Protein '", protein, "' was excluded during fitting (likely insufficient data/low counts).")
  }
  
  # --- 2. Retrieve Specific Data ---
  # Get the raw values for this protein from the stored matrix
  vals <- as.numeric(mix_model_fit$data[protein, ])
  vals <- vals[is.finite(vals)]
  
  # Get the specific model fit object
  key <- paste0("ev", if(ev) "T" else "F", "_ncomp_", ncomp)
  fit <- mix_model_fit$fits[[protein]]$all_fits[[key]]
  
  if (is.null(fit)) {
    stop(paste("No fit found for", protein, "with parameters:", key))
  }
  
  thresh <- mix_model_fit$threshold
  
  # --- 3. Generate Densities (Plotting Logic) ---
  # Re-calculate density curves based on the stored parameters
  
  # Shift mean (alignment)
  mix_mean <- sum(fit$pi * fit$mu)
  fit$mu   <- fit$mu + (mean(vals) - mix_mean)
  
  x_vals <- seq(min(vals), max(vals), length.out = 1000)
  
  comp_mat <- vapply(
    seq_along(fit$mu),
    function(i) fit$pi[i] * dnorm(x_vals, mean = fit$mu[i], sd = fit$sd[i]),
    numeric(length(x_vals))
  )
  
  mix_y <- rowSums(comp_mat)
  
  # Data frames for ggplot
  mix_df  <- tibble(x = x_vals, density = mix_y)
  comp_df <- bind_rows(lapply(seq_along(fit$mu), function(i) {
    tibble(x = x_vals, density = comp_mat[, i], comp = paste0("Component ", i))
  }))
  obs_df  <- tibble(val = vals)
  
  # --- 4. Plot ---
  ggplot(obs_df, aes(x = val)) +
    geom_histogram(aes(y = after_stat(density)), bins = 70, fill = "grey80", color = "grey40") +
    geom_area(data = comp_df, aes(x = x, y = density, fill = comp), alpha = 0.45, position = "identity") +
    geom_line(data = mix_df, aes(x = x, y = density), linewidth = 1.1, color = "black") +
    geom_vline(xintercept = thresh, linetype = "dashed", color = "red") +
    annotate("text", x = thresh, y = max(mix_df$density) * 1.05, 
             label = "Neg Ctrl + 1SD", color = "red", angle = 90, vjust = -0.5, size = 3.5) +
    labs(
      title = paste0(protein, " (ev=", ev, ", ncomp=", ncomp, ")"),
      subtitle = paste0("BIC: ", round(fit$BIC, 2)),
      x = "log10(Q3-normalized count)", y = "Density", fill = "Components"
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold", hjust = 0.5))
}