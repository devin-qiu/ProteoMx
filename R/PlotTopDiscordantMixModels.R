#' Plot Mixture Models for Top Discordant Proteins
#'
#' Plots mixture model fits for a vector of proteins and combines them into
#' a multi-panel figure with a user-defined number of plots per row.
#'
#' @param geomx_set A \code{NanoStringGeoMxSet} object processed by \code{MixModelFit}.
#' @param proteins Character vector. Protein names to plot, e.g. output from
#'   \code{FindTopDiscordantProteins}.
#' @param assay Character. Assay to use for plotting. Default is \code{"q_norm"}.
#' @param neg_ctrl Character vector. Negative control probe(s) used for thresholding.
#' @param ncol Integer. Number of plots per row. Default is 3.
#' @param ncomp Integer or NULL. Passed to \code{PlotMixModel}.
#' @param ev Logical or NULL. Passed to \code{PlotMixModel}.
#'
#' @return A combined patchwork ggplot object.
#'
#' @export
PlotTopDiscordantMixModels <- function(
    geomx_set,
    proteins,
    assay = "q_norm",
    neg_ctrl = "Rt IgG2a",
    ncol = 3,
    ncomp = NULL,
    ev = NULL
) {
  
  require(patchwork)
  
  if (missing(proteins) || is.null(proteins) || length(proteins) == 0) {
    stop("Please provide a non-empty character vector of proteins.")
  }
  
  proteins <- as.character(proteins)
  
  missing_proteins <- setdiff(proteins, rownames(geomx_set))
  
  if (length(missing_proteins) > 0) {
    stop(
      "The following protein(s) were not found in geomx_set: ",
      paste(missing_proteins, collapse = ", ")
    )
  }
  
  plot_list <- lapply(proteins, function(protein_name) {
    PlotMixModel(
      geomx_set = geomx_set,
      protein = protein_name,
      assay = assay,
      ncomp = ncomp,
      ev = ev,
      neg_ctrl = neg_ctrl
    )
  })
  
  names(plot_list) <- proteins
  
  combined_plot <- patchwork::wrap_plots(
    plot_list,
    ncol = ncol
  )
  
  return(combined_plot)
}