#' Check Probe Performance and Dynamic Range
#'
#' Evaluates the dynamic range and performance of a specific antibody (protein)
#' by plotting its relative expression across all AOIs, ordered by intensity.
#'
#' @description
#' Generates an ordered scatter plot with a fitted LOESS curve. 
#' The Y-axis represents the relative expression (AOI Count / Max Count), 
#' and the X-axis represents the AOIs ranked from lowest to highest expression.
#' The points are color-coded based on a user-provided metadata variable to 
#' assess if antibody performance correlates with specific biological or technical factors.
#'
#' @param geomx_set A \code{NanoStringGeoMxSet} object.
#' @param protein Character. The name of the target protein (antibody) to evaluate.
#' @param plot_by Character. The column name in \code{phenoData(geomx_set)} used to color-code the scatter plot.
#' @param assay Character. The assay to use for expression values (default: "q_norm").
#'
#' @return A \code{\link[ggplot2]{ggplot}} object showing the performance curve.
#'
#' @import ggplot2
#' @importFrom Biobase assayDataElement pData
#' @importFrom dplyr arrange mutate row_number
#' @export
#'
#' @examples
#' \dontrun{
#'   # Check CD4 performance, color-coded by Tissue Type
#'   CheckProbePerformance(geomx_set, protein = "CD4", plot_by = "TissueType")
#'   
#'   # Check CD4 performance, color-coded by X coordinate
#'   CheckProbePerformance(geomx_set, protein = "CD4", plot_by = "X_Coordinate")
#' }
CheckProbePerformance <- function(geomx_set, protein, plot_by, assay = "q_norm") {
  
  require(ggplot2)
  require(Biobase)
  require(dplyr)
  
  # --- 1. Validate Inputs ---
  if (!assay %in% names(assayData(geomx_set))) {
    stop(paste("Error: Assay '", assay, "' not found in the object."))
  }
  
  if (!protein %in% rownames(geomx_set)) {
    stop(paste("Error: Protein '", protein, "' not found in the dataset."))
  }
  
  if (!plot_by %in% colnames(pData(geomx_set))) {
    stop(paste("Error: Variable '", plot_by, "' not found in phenoData of the object."))
  }
  
  # --- 2. Extract Data ---
  expr_mat <- assayDataElement(geomx_set, assay)
  expr_vals <- as.numeric(expr_mat[protein, ])
  
  pheno_df <- pData(geomx_set)
  
  # Combine into a plotting dataframe
  plot_df <- data.frame(
    AOI = colnames(expr_mat),
    Expression = expr_vals,
    ColorVar = pheno_df[[plot_by]],
    stringsAsFactors = FALSE
  )
  
  # Drop NAs safely
  plot_df <- plot_df[!is.na(plot_df$Expression), ]
  
  # --- 3. Process Data (Order and Relative Scale) ---
  max_val <- max(plot_df$Expression, na.rm = TRUE)
  
  if (max_val == 0) {
    warning("Maximum expression is 0. Data is entirely flat.")
    max_val <- 1 # Prevent division by zero mathematically
  }
  
  plot_df <- plot_df %>%
    arrange(Expression) %>%
    mutate(
      Rank = row_number(),
      Rel_Expression = Expression / max_val
    )
  
  # Check if the plot_by variable is continuous or categorical
  is_numeric_color <- is.numeric(plot_df$ColorVar)
  
  # --- 4. Generate Plot ---
  p <- ggplot(plot_df, aes(x = Rank, y = Rel_Expression)) +
    geom_point(aes(color = ColorVar), size = 2.5, alpha = 0.8) +
    
    # Add the trendline (LOESS provides a nice, smooth S-curve fit)
    geom_smooth(method = "loess", se = FALSE, color = "black", linetype = "dashed", size = 0.8) +
    
    labs(
      title = paste("Antibody Dynamic Range:", protein),
      subtitle = paste("Color-coded by:", plot_by, "| Assay:", assay),
      x = "AOIs (Ranked by Expression Level)",
      y = "Relative Expression (AOI Read / Max Read)",
      color = plot_by
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11, color = "grey30"),
      axis.title = element_text(face = "bold"),
      legend.position = "right"
    )
  
  # --- 5. Smart Color Scaling ---
  if (is_numeric_color) {
    # Continuous gradient for numeric data (like coordinates)
    p <- p + scale_color_viridis_c(option = "plasma")
  } else {
    # Distinct colors for categorical data (like Tissue Type or Sample ID)
    p <- p + scale_color_viridis_d(option = "turbo")
  }
  
  return(p)
}