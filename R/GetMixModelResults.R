#' Get Mixture Model Results
#'
#' A helper function to safely retrieve the mixture model analysis results 
#' stored inside a GeoMxSet object.
#'
#' @param geomx_set A \code{NanoStringGeoMxSet} object processed by \code{\link{MixModelFit}}.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{fits}: The list of model fits for all proteins.
#'   \item \code{threshold}: The calculated background threshold.
#'   \item \code{params}: The parameters used for the run (ncomps, negative control).
#'   \item \code{timestamp}: When the analysis was run.
#' }
#'
#' @details
#' This function accesses \code{experimentData(object)@other$MixModel}. 
#' It checks if the data exists and stops with a helpful error message if 
#' \code{MixModelFit} has not been run yet.
#'
#' @importFrom Biobase experimentData
#' @export
#'
#' @seealso \code{\link{MixModelFit}}
#'
#' @examples
#' \dontrun{
#'   # Retrieve the results list
#'   results <- GetMixModelResults(geomx_set)
#'   
#'   # Inspect the threshold used
#'   print(results$threshold)
#' }
GetMixModelResults <- function(geomx_set) {
  res <- experimentData(geomx_set)@other$MixModel
  if (is.null(res)) {
    stop("No Mixture Model results found. Please run MixModelFit() first.")
  }
  return(res)
}