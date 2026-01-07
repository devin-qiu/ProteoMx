#' Get Mixture Model Results
#'
#' Helper function to retrieve MixModel results from a GeoMxSet object.
#'
#' @param geomx_set A GeoMxSet object processed by MixModelFit.
#' @return The list of fits and thresholds.
#' @export
GetMixModelResults <- function(geomx_set) {
  res <- experimentData(geomx_set)@other$MixModel
  if (is.null(res)) {
    stop("No Mixture Model results found. Please run MixModelFit() first.")
  }
  return(res)
}