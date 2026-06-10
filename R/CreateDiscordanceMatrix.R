#' Create AOI-normalized Protein Discordance Matrix
#'
#' Computes a protein-by-AOI discordance matrix by comparing the AOI rank of each
#' protein against the AOI rank of a selected negative probe.
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
#' @return A matrix with proteins as rows and AOIs as columns, containing AOI-normalized
#'   rank-discordance values.
#'
#' @importFrom Biobase assayData assayDataElement
#' @export
CreateDiscordanceMatrix <- function(
    geomx_set,
    neg_probe,
    assay = "q_norm",
    proteins = NULL
) {
  
  require(Biobase)
  
  # --- 1. Validate inputs ---
  if (!assay %in% names(assayData(geomx_set))) {
    stop(paste("Error: Assay '", assay, "' not found in the object."))
  }
  
  expr_mat <- assayDataElement(geomx_set, assay)
  
  if (!neg_probe %in% rownames(expr_mat)) {
    stop(paste("Error: Negative probe '", neg_probe, "' not found in the dataset."))
  }
  
  common_negs <- c(
    "Hmr IgG",
    "Ms IgG2b",
    "Ms IgG2a",
    "Rb IgG",
    "Rt IgG2a",
    "Ms IgG1"
  )
  
  if (is.null(proteins)) {
    proteins <- setdiff(rownames(expr_mat), common_negs)
  } else {
    missing_proteins <- setdiff(proteins, rownames(expr_mat))
    
    if (length(missing_proteins) > 0) {
      stop(
        "The following protein(s) were not found in the dataset: ",
        paste(missing_proteins, collapse = ", ")
      )
    }
    
    proteins <- setdiff(proteins, neg_probe)
  }
  
  if (length(proteins) == 0) {
    stop("No valid proteins remain for discordance matrix calculation.")
  }
  
  # --- 2. Subset expression matrix ---
  target_mat <- expr_mat[proteins, , drop = FALSE]
  neg_vals <- as.numeric(expr_mat[neg_probe, ])
  
  # --- 3. Rank AOIs for the negative probe ---
  # ties.method = "average" gives stable behavior when AOIs have identical values
  neg_rank <- rank(neg_vals, ties.method = "average", na.last = "keep")
  
  # --- 4. Rank AOIs independently for each protein ---
  protein_rank_mat <- t(apply(target_mat, 1, function(x) {
    rank(as.numeric(x), ties.method = "average", na.last = "keep")
  }))
  
  rownames(protein_rank_mat) <- proteins
  colnames(protein_rank_mat) <- colnames(expr_mat)
  
  # --- 5. Compute absolute rank difference d_ij ---
  rankdiff_mat <- sweep(
    protein_rank_mat,
    2,
    neg_rank,
    FUN = "-"
  )
  
  rankdiff_mat <- abs(rankdiff_mat)
  
  # --- 6. Compute AOI-wise denominator D_j + 1 ---
  aoi_rankdiff_sum <- colSums(rankdiff_mat, na.rm = TRUE)
  
  # --- 7. Normalize each protein's rank difference within each AOI ---
  discordance_mat <- sweep(
    rankdiff_mat,
    2,
    aoi_rankdiff_sum + 1,
    FUN = "/"
  )
  
  rownames(discordance_mat) <- proteins
  colnames(discordance_mat) <- colnames(expr_mat)
  
  return(discordance_mat)
}