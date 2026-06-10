#' Check Probe Performance and Dynamic Range Against Negative Probe(s)
#'
#' Evaluates the dynamic range and performance of a target protein/probe by
#' comparing its AOI-level ranked relative signal against one or more negative
#' probes. Target and negative probes are ranked independently, while linkage
#' lines connect matched AOIs between the target and negative probe within each
#' panel.
#'
#' @param geomx_set A \code{NanoStringGeoMxSet} object.
#' @param protein Character. The name of the target protein/probe to evaluate.
#' @param neg_probe Character vector or NULL. Negative probe(s) to compare against.
#'   If NULL, all common negative probes found in the object are used.
#' @param plot_by Character. Column name in \code{phenoData(geomx_set)} used for color coding.
#' @param assay Character. Assay to use for expression values. Default is "q_norm".
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
    assay = "q_norm"
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
  
  target_df <- data.frame(
    AOI = colnames(expr_mat),
    Target_Expression = as.numeric(expr_mat[protein, ]),
    ColorVar = pheno_df[[plot_by]],
    stringsAsFactors = FALSE
  )
  
  target_df <- target_df[!is.na(target_df$Target_Expression), ]
  
  target_max <- max(target_df$Target_Expression, na.rm = TRUE)
  
  if (!is.finite(target_max) || target_max == 0) {
    warning("Maximum target-protein expression is 0 or non-finite. Using 1 to prevent division by zero.")
    target_max <- 1
  }
  
  target_df <- target_df %>%
    arrange(Target_Expression) %>%
    mutate(
      Target_Rank = row_number(),
      Target_Rel_Expression = Target_Expression / target_max
    )
  
  panel_list <- lapply(neg_probe, function(neg_name) {
    
    neg_df <- data.frame(
      AOI = colnames(expr_mat),
      Neg_Expression = as.numeric(expr_mat[neg_name, ]),
      stringsAsFactors = FALSE
    )
    
    neg_df <- neg_df[!is.na(neg_df$Neg_Expression), ]
    
    neg_max <- max(neg_df$Neg_Expression, na.rm = TRUE)
    
    if (!is.finite(neg_max) || neg_max == 0) {
      warning(
        "Maximum expression for negative probe '",
        neg_name,
        "' is 0 or non-finite. Using 1 to prevent division by zero."
      )
      neg_max <- 1
    }
    
    neg_df <- neg_df %>%
      arrange(Neg_Expression) %>%
      mutate(
        Neg_Rank = row_number(),
        Neg_Rel_Expression = Neg_Expression / neg_max
      )
    
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
  
  target_points <- data.frame(
    AOI = plot_df$AOI,
    Rank = plot_df$Target_Rank,
    Rel_Expression = plot_df$Target_Rel_Expression,
    ColorVar = plot_df$ColorVar,
    Neg_Probe = plot_df$Neg_Probe,
    Probe_Type = protein,
    stringsAsFactors = FALSE
  )
  
  neg_points <- data.frame(
    AOI = plot_df$AOI,
    Rank = plot_df$Neg_Rank,
    Rel_Expression = plot_df$Neg_Rel_Expression,
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
      Rel_Expression = plot_df$Target_Rel_Expression,
      Neg_Probe = plot_df$Neg_Probe,
      stringsAsFactors = FALSE
    ),
    data.frame(
      AOI = plot_df$AOI,
      Rank = plot_df$Neg_Rank,
      Rel_Expression = plot_df$Neg_Rel_Expression,
      Neg_Probe = plot_df$Neg_Probe,
      stringsAsFactors = FALSE
    )
  )
  
  is_numeric_color <- is.numeric(plot_df$ColorVar)
  
  p <- ggplot() +
    geom_line(
      data = link_df,
      aes(
        x = Rank,
        y = Rel_Expression,
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
        y = Rel_Expression,
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
        y = Rel_Expression,
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
        y = Rel_Expression,
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
        assay
      ),
      x = "AOI Rank Within Each Probe",
      y = "Relative Expression (Probe Read / Probe Max Read)",
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