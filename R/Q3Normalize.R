#' Q3 Normalization Helper
#'
#' Performs 75th percentile (Q3) normalization on a GeoMxSet object.
#' Stores the result in assayDataElement(object, "q_norm").
#'
#' @param geomx_set A GeoMxSet object.
#' @return A GeoMxSet object with the 'q_norm' assay data added.
#' @export
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