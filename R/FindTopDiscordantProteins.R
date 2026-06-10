#' Find Top Discordant Proteins
#'
#' Identifies proteins with the highest total discordance score from a
#' protein-by-AOI discordance matrix.
#'
#' @param discordance_mat A numeric matrix returned by \code{CreateDiscordanceMatrix},
#'   with proteins as rows and AOIs as columns.
#' @param n Integer. Number of top discordant proteins to return. Default is 10.
#'
#' @return A character vector containing the names of the top \code{n} proteins
#'   with the highest row-summed discordance scores.
#'
#' @export
FindTopDiscordantProteins <- function(discordance_mat, n = 10) {
  
  if (!is.matrix(discordance_mat) && !is.data.frame(discordance_mat)) {
    stop("discordance_mat must be a matrix or data.frame.")
  }
  
  discordance_mat <- as.matrix(discordance_mat)
  
  if (is.null(rownames(discordance_mat))) {
    stop("discordance_mat must have protein names as rownames.")
  }
  
  if (!is.numeric(discordance_mat)) {
    stop("discordance_mat must be numeric.")
  }
  
  if (!is.numeric(n) || length(n) != 1 || n <= 0) {
    stop("n must be a single positive number.")
  }
  
  n <- floor(n)
  
  protein_scores <- rowSums(discordance_mat, na.rm = TRUE)
  
  top_proteins <- names(
    sort(protein_scores, decreasing = TRUE)
  )[seq_len(min(n, length(protein_scores)))]
  
  return(top_proteins)
}