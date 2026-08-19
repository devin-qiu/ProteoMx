#' Check Probe Performance and Dynamic Range Against Negative Probe(s)
#'
#' Evaluates the dynamic range and performance of a target protein/probe by
#' comparing its AOI-level ranked assay value against one or more negative
#' probes. Target and negative probes are ranked independently, while linkage
#' lines connect matched AOIs between the target and negative probe within each
#' panel.
#'
#' This version plots actual assay values, e.g. q_norm, rather than relative
#' expression normalized to each probe's maximum.
#'
#' @param geomx_set A \code{NanoStringGeoMxSet} object.
#' @param protein Character. The name of the target protein/probe to evaluate.
#' @param neg_probe Character vector or NULL. Negative probe(s) to compare against.
#'   If NULL, all common negative probes found in the object are used.
#' @param plot_by Character. Column name in \code{phenoData(geomx_set)} used for color coding.
#' @param assay Character. Assay to use for expression values. Default is "q_norm".
#' @param log2_y Logical. If TRUE, plot log2(assay value + pseudocount). Default FALSE.
#' @param pseudocount Numeric. Pseudocount used when log2_y = TRUE. Default 1.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object.
#'
#' @import ggplot2
#' @importFrom Biobase assayData assayDataElement pData
#' @importFrom dplyr arrange mutate row_number
#' @export
CheckProbePerformance <- function(
    geomx_set,
    protein,
    neg_probe = NULL,
    plot_by,
    assay = "q_norm",
    log2_y = FALSE,
    pseudocount = 1
) {
  
  require(ggplot2)
  require(Biobase)
  require(dplyr)
  
  if (!assay %in% names(assayData(geomx_set))) {
    stop(paste("Error: Assay '", assay, "' not found in the object."))
  }
  
  if (!protein %in% rownames(geomx_set)) {
    stop(paste("Error: Protein '", protein, "' not found in the dataset."))
  }
  
  if (!plot_by %in% colnames(pData(geomx_set))) {
    stop(paste("Error: Variable '", plot_by, "' not found in phenoData of the object."))
  }
  
  common_negs <- c(
    "Hmr IgG",
    "Ms IgG2b",
    "Ms IgG2a",
    "Rb IgG",
    "Rt IgG2a",
    "Ms IgG1"
  )
  
  if (is.null(neg_probe)) {
    neg_probe <- intersect(common_negs, rownames(geomx_set))
    
    if (length(neg_probe) == 0) {
      stop("No common negative probes were found in the dataset. Please specify neg_probe manually.")
    }
  } else {
    missing_negs <- setdiff(neg_probe, rownames(geomx_set))
    
    if (length(missing_negs) > 0) {
      stop(
        "The following negative probe(s) were not found in the dataset: ",
        paste(missing_negs, collapse = ", ")
      )
    }
  }
  
  expr_mat <- assayDataElement(geomx_set, assay)
  pheno_df <- pData(geomx_set)
  
  # Match phenotype rows to assay columns
  pheno_df <- pheno_df[colnames(expr_mat), , drop = FALSE]
  
  if (!all(rownames(pheno_df) == colnames(expr_mat))) {
    stop("phenoData rownames do not match assay matrix column names.")
  }
  
  # ------------------------------------------------------------
  # Target probe values
  # ------------------------------------------------------------
  
  target_values <- as.numeric(expr_mat[protein, ])
  
  if (log2_y) {
    target_values[target_values < 0] <- NA
    target_values <- log2(target_values + pseudocount)
  }
  
  target_df <- data.frame(
    AOI = colnames(expr_mat),
    Target_Value = target_values,
    ColorVar = pheno_df[[plot_by]],
    stringsAsFactors = FALSE
  )
  
  target_df <- target_df[is.finite(target_df$Target_Value), ]
  
  target_df <- target_df %>%
    arrange(Target_Value) %>%
    mutate(Target_Rank = row_number())
  
  # ------------------------------------------------------------
  # Negative probe values
  # ------------------------------------------------------------
  
  panel_list <- lapply(neg_probe, function(neg_name) {
    
    neg_values <- as.numeric(expr_mat[neg_name, ])
    
    if (log2_y) {
      neg_values[neg_values < 0] <- NA
      neg_values <- log2(neg_values + pseudocount)
    }
    
    neg_df <- data.frame(
      AOI = colnames(expr_mat),
      Neg_Value = neg_values,
      stringsAsFactors = FALSE
    )
    
    neg_df <- neg_df[is.finite(neg_df$Neg_Value), ]
    
    neg_df <- neg_df %>%
      arrange(Neg_Value) %>%
      mutate(Neg_Rank = row_number())
    
    merged_df <- merge(
      target_df,
      neg_df,
      by = "AOI",
      all = FALSE
    )
    
    merged_df$Neg_Probe <- neg_name
    
    merged_df
  })
  
  plot_df <- do.call(rbind, panel_list)
  
  if (nrow(plot_df) == 0) {
    stop("No matched AOIs remain after merging target and negative-probe values.")
  }
  
  # ------------------------------------------------------------
  # Plotting data
  # ------------------------------------------------------------
  
  target_points <- data.frame(
    AOI = plot_df$AOI,
    Rank = plot_df$Target_Rank,
    Value = plot_df$Target_Value,
    ColorVar = plot_df$ColorVar,
    Neg_Probe = plot_df$Neg_Probe,
    Probe_Type = protein,
    stringsAsFactors = FALSE
  )
  
  neg_points <- data.frame(
    AOI = plot_df$AOI,
    Rank = plot_df$Neg_Rank,
    Value = plot_df$Neg_Value,
    ColorVar = plot_df$ColorVar,
    Neg_Probe = plot_df$Neg_Probe,
    Probe_Type = plot_df$Neg_Probe,
    stringsAsFactors = FALSE
  )
  
  point_df <- rbind(target_points, neg_points)
  
  link_df <- rbind(
    data.frame(
      AOI = plot_df$AOI,
      Rank = plot_df$Target_Rank,
      Value = plot_df$Target_Value,
      Neg_Probe = plot_df$Neg_Probe,
      stringsAsFactors = FALSE
    ),
    data.frame(
      AOI = plot_df$AOI,
      Rank = plot_df$Neg_Rank,
      Value = plot_df$Neg_Value,
      Neg_Probe = plot_df$Neg_Probe,
      stringsAsFactors = FALSE
    )
  )
  
  is_numeric_color <- is.numeric(plot_df$ColorVar)
  
  y_label <- if (log2_y) {
    paste0("log2(", assay, " + ", pseudocount, ")")
  } else {
    assay
  }
  
  # ------------------------------------------------------------
  # Plot
  # ------------------------------------------------------------
  
  p <- ggplot() +
    geom_line(
      data = link_df,
      aes(
        x = Rank,
        y = Value,
        group = interaction(Neg_Probe, AOI)
      ),
      color = "grey70",
      alpha = 0.35,
      linewidth = 0.3
    ) +
    geom_point(
      data = target_points,
      aes(
        x = Rank,
        y = Value,
        color = ColorVar
      ),
      size = 2.5,
      alpha = 1,
      shape = 16
    ) +
    geom_point(
      data = neg_points,
      aes(
        x = Rank,
        y = Value,
        color = ColorVar
      ),
      size = 1.25,
      alpha = 0.5,
      shape = 17
    ) +
    geom_smooth(
      data = point_df,
      aes(
        x = Rank,
        y = Value,
        linetype = Probe_Type
      ),
      method = "loess",
      se = FALSE,
      color = "black",
      linewidth = 0.8
    ) +
    facet_wrap(~ Neg_Probe, nrow = 1) +
    labs(
      title = paste("Probe Performance:", protein, "vs Negative Probe(s)"),
      subtitle = paste(
        "Target and negative probes ranked independently | Color-coded by:",
        plot_by,
        "| Assay:",
        assay,
        "| Y-axis:",
        y_label
      ),
      x = "AOI Rank Within Each Probe",
      y = y_label,
      color = plot_by,
      linetype = "Probe"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11, color = "grey30"),
      axis.title = element_text(face = "bold"),
      legend.position = "right"
    )
  
  if (is_numeric_color) {
    p <- p + scale_color_viridis_c(option = "plasma")
  } else {
    p <- p + scale_color_viridis_d(option = "turbo")
  }
  
  return(p)
}