#' Plot Relative Expression vs Discordance Score
#'
#' Creates scatter plot(s) comparing relative expression against AOI-normalized
#' discordance score for one or more proteins. Expression values are taken from
#' a selected assay, while discordance scores are taken from a stored discordance
#' matrix in \code{experimentData(geomx_set)@other$DiscordanceMatrix}.
#'
#' The stored discordance matrix is selected using \code{neg_ctrl}. For example,
#' if \code{neg_ctrl = "Rb IgG"}, the function will look for
#' \code{discordance_Rb_IgG}. Spaces in \code{neg_ctrl} are converted to
#' underscores.
#'
#' For each selected protein, expression values are scaled by that protein's
#' maximum expression across AOIs:
#'
#' Relative_Expression_j = Expression_j / max(Expression)
#'
#' Points with discordance scores less than or equal to the user-defined threshold
#' are colored red and plotted on top of the black points. If multiple proteins
#' are supplied, plots are combined into a multi-panel figure with 3 columns.
#'
#' @param geomx_set A \code{NanoStringGeoMxSet} object.
#' @param protein Character vector. One or more proteins to plot.
#' @param neg_ctrl Character. Negative control probe used to select the stored
#'   discordance matrix. Default is \code{"Rb IgG"}.
#' @param assay Character. Assay to use for expression values. Default is \code{"exprs"}.
#' @param threshold Numeric. Discordance-score cutoff used to color low-discordance
#'   points. Default is \code{0.01}.
#' @param ncol Integer. Number of plot columns when multiple proteins are supplied.
#'   Default is 3.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object if one protein is supplied, or
#'   a combined patchwork object if multiple proteins are supplied.
#'
#' @import ggplot2
#' @importFrom Biobase assayData assayDataElement experimentData
#' @export
PlotExpressionDiscordance <- function(
    geomx_set,
    protein,
    neg_ctrl = "Rb IgG",
    assay = "exprs",
    threshold = 0.0005,
    ncol = 3
) {
  
  require(Biobase)
  require(ggplot2)
  
  if (length(protein) > 1) {
    require(patchwork)
  }
  
  # --- 1. Validate assay ---
  if (!assay %in% names(assayData(geomx_set))) {
    stop(paste("Error: Assay '", assay, "' not found in the object."))
  }
  
  expr_mat <- assayDataElement(geomx_set, assay)
  
  missing_proteins <- setdiff(protein, rownames(expr_mat))
  if (length(missing_proteins) > 0) {
    stop(
      "The following protein(s) were not found in assay '", assay, "': ",
      paste(missing_proteins, collapse = ", ")
    )
  }
  
  # --- 2. Construct and validate discordance matrix name ---
  neg_ctrl_clean <- gsub("\\s+", "_", neg_ctrl)
  discordance_name <- paste0("discordance_", neg_ctrl_clean)
  
  discordance_list <- experimentData(geomx_set)@other$DiscordanceMatrix
  
  if (is.null(discordance_list)) {
    stop("No DiscordanceMatrix found in experimentData(geomx_set)@other.")
  }
  
  if (!discordance_name %in% names(discordance_list)) {
    stop(
      "Discordance matrix '", discordance_name, "' not found. Available matrices: ",
      paste(names(discordance_list), collapse = ", ")
    )
  }
  
  discordance_mat <- discordance_list[[discordance_name]]
  
  missing_discordance_proteins <- setdiff(protein, rownames(discordance_mat))
  if (length(missing_discordance_proteins) > 0) {
    stop(
      "The following protein(s) were not found in discordance matrix '",
      discordance_name,
      "': ",
      paste(missing_discordance_proteins, collapse = ", ")
    )
  }
  
  # --- 3. Internal single-protein plotting function ---
  plot_one_protein <- function(protein_name) {
    
    common_aois <- intersect(colnames(expr_mat), colnames(discordance_mat))
    
    if (length(common_aois) == 0) {
      stop("No shared AOIs found between expression matrix and discordance matrix.")
    }
    
    expr_vals <- as.numeric(expr_mat[protein_name, common_aois])
    discordance_vals <- as.numeric(discordance_mat[protein_name, common_aois])
    
    max_expr <- max(expr_vals, na.rm = TRUE)
    
    if (!is.finite(max_expr) || max_expr == 0) {
      stop(
        "Maximum expression for protein '",
        protein_name,
        "' is 0 or non-finite; cannot compute relative expression."
      )
    }
    
    rel_expr_vals <- expr_vals / max_expr
    
    plot_df <- data.frame(
      Protein = protein_name,
      AOI = common_aois,
      Relative_Expression = rel_expr_vals,
      Discordance = discordance_vals,
      Below_Threshold = discordance_vals <= threshold,
      stringsAsFactors = FALSE
    )
    
    plot_df <- plot_df[
      is.finite(plot_df$Relative_Expression) &
        is.finite(plot_df$Discordance),
    ]
    
    if (nrow(plot_df) == 0) {
      stop(
        "No finite relative-expression and discordance-score pairs remain for protein '",
        protein_name,
        "'."
      )
    }
    
    black_df <- plot_df[!plot_df$Below_Threshold, ]
    red_df <- plot_df[plot_df$Below_Threshold, ]
    
    ggplot() +
      geom_point(
        data = black_df,
        aes(x = Relative_Expression, y = Discordance),
        color = "black",
        alpha = 0.45,
        size = 1.8
      ) +
      geom_point(
        data = red_df,
        aes(x = Relative_Expression, y = Discordance),
        color = "red",
        alpha = 0.85,
        size = 1.8
      ) +
      geom_hline(
        yintercept = threshold,
        linetype = "dashed",
        color = "red",
        linewidth = 0.6
      ) +
      labs(
        title = protein_name,
        subtitle = paste0(
          "Neg ctrl: ", neg_ctrl,
          " | Assay: ", assay,
          " | Red: Discordance <= ", threshold
        ),
        x = paste0("Relative Expression (", assay, " / max ", assay, ")"),
        y = "Discordance Score"
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10, color = "grey30"),
        axis.title = element_text(face = "bold")
      )
  }
  
  # --- 4. Generate one or multiple plots ---
  plot_list <- lapply(protein, plot_one_protein)
  names(plot_list) <- protein
  
  if (length(plot_list) == 1) {
    return(plot_list[[1]])
  }
  
  combined_plot <- patchwork::wrap_plots(
    plot_list,
    ncol = ncol
  )
  
  return(combined_plot)
}