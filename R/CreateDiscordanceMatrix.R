#' Create AOI-normalized Protein Discordance Matrix
#'
#' Computes a protein-by-AOI discordance matrix by comparing the AOI rank of each
#' protein against the AOI rank of a selected negative probe. The resulting matrix
#' is stored in \code{experimentData(geomx_set)@other$DiscordanceMatrix} using the
#' name \code{discordance_<neg_probe>}. If the selected negative-probe name contains
#' spaces, spaces are replaced with underscores in the stored matrix name.
#'
#' For protein i and AOI j:
#'
#' d_ij = abs(rank_i_j - rank_neg_j)
#'
#' Then, within each AOI j:
#'
#' p_ij = d_ij / (sum_i d_ij + 1)
#'
#' The output value p_ij represents the fraction of AOI-level rank discordance
#' attributable to protein i, stabilized by adding 1 to the AOI-wise denominator.
#'
#' @param geomx_set A \code{NanoStringGeoMxSet} object.
#' @param neg_probe Character. The negative-control probe used as the rank reference.
#' @param assay Character. Assay to use for expression values. Default is \code{"q_norm"}.
#' @param proteins Character vector or NULL. Proteins to include. If NULL, all non-negative
#'   probes except the selected negative probe are used.
#'
#' @return The input \code{NanoStringGeoMxSet} object with the discordance matrix
#'   saved in \code{experimentData(geomx_set)@other$DiscordanceMatrix}.
#'
#' @importFrom Biobase assayData assayDataElement experimentData experimentData<-
#' @export
CreateDiscordanceMatrix <- function(
    geomx_set,
    neg_probe,
    assay = "q_norm",
    proteins = NULL
) {
  
  require(Biobase)
  
  if (!assay %in% names(assayData(geomx_set))) {
    stop(paste("Error: Assay '", assay, "' not found in the object."))
  }
  
  expr_mat <- assayDataElement(geomx_set, assay)
  
  if (!neg_probe %in% rownames(expr_mat)) {
    stop(paste("Error: Negative probe '", neg_probe, "' not found in the dataset."))
  }
  
  if (is.null(proteins)) {
    proteins <- rownames(expr_mat)
  } else {
    missing_proteins <- setdiff(proteins, rownames(expr_mat))
    
    if (length(missing_proteins) > 0) {
      stop(
        "The following protein(s) were not found in the dataset: ",
        paste(missing_proteins, collapse = ", ")
      )
    }
  }
  
  if (length(proteins) == 0) {
    stop("No valid proteins remain for discordance matrix calculation.")
  }
  
  target_mat <- expr_mat[proteins, , drop = FALSE]
  neg_vals <- as.numeric(expr_mat[neg_probe, ])
  
  neg_rank <- rank(neg_vals, ties.method = "average", na.last = "keep")
  
  protein_rank_mat <- t(apply(target_mat, 1, function(x) {
    rank(as.numeric(x), ties.method = "average", na.last = "keep")
  }))
  
  rownames(protein_rank_mat) <- proteins
  colnames(protein_rank_mat) <- colnames(expr_mat)
  
  rankdiff_mat <- sweep(
    protein_rank_mat,
    2,
    neg_rank,
    FUN = "-"
  )
  
  rankdiff_mat <- abs(rankdiff_mat)
  
  aoi_rankdiff_sum <- colSums(rankdiff_mat, na.rm = TRUE)
  
  discordance_mat <- sweep(
    rankdiff_mat,
    2,
    aoi_rankdiff_sum + 1,
    FUN = "/"
  )
  
  rownames(discordance_mat) <- proteins
  colnames(discordance_mat) <- colnames(expr_mat)
  
  neg_probe_clean <- gsub("\\s+", "_", neg_probe)
  discordance_name <- paste0("discordance_", neg_probe_clean)
  
  exp_data <- experimentData(geomx_set)
  
  if (is.null(exp_data@other$DiscordanceMatrix)) {
    exp_data@other$DiscordanceMatrix <- list()
  }
  
  exp_data@other$DiscordanceMatrix[[discordance_name]] <- discordance_mat
  
  experimentData(geomx_set) <- exp_data
  
  return(geomx_set)
}