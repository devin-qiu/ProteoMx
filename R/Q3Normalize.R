#' This is a required preprocessing step before running mixture models.
#'
#' @param geomx_set A \code{NanoStringGeoMxSet} object.
#'
#' @return A \code{NanoStringGeoMxSet} object with the normalized matrix stored in 
#'   \code{assayDataElement(object, "q_norm")}.
#'
#' @details
#' This function checks if \code{q_norm} already exists. If not, it calls
#' \code{NanoStringNCTools::normalize} with \code{norm_method="quant"} and 
#' \code{desiredQuantile=0.75}.
#'
#' @importFrom NanoStringNCTools normalize
#' @importFrom Biobase assayData
#' @export
#'
#' @examples
#' \dontrun{
#'   # Basic usage
#'   geomx_set <- Q3Normalize(geomx_set)
#' }
Q3Normalize <- function(geomx_set) {
  
  require(NanoStringNCTools)
  require(GeomxTools)
  
  message("Performing Q3 (75th percentile) Normalization...")
  
  # Check if it already exists to avoid redundant computation (optional safety)
  if ("q_norm" %in% names(assayData(geomx_set))) {
    message("Note: 'q_norm' already exists. Overwriting...")
  }
  
  # Run standard normalization
  tryCatch({
    geomx_set <- NanoStringNCTools::normalize(
      geomx_set,
      norm_method = "quant", 
      desiredQuantile = 0.75,
      toElt = "q_norm"
    )
  }, error = function(e) {
    stop("Normalization failed. Please check input object structure.\nError: ", e$message)
  })
  
  message("Success! 'q_norm' added to assayData.")
  return(geomx_set)
}