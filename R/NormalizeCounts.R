#' Normalize GeoMx Protein or RNA Counts
#'
#' This is a required preprocessing step before running mixture models or 
#' differential expression analysis. It allows for multiple normalization 
#' strategies commonly used in spatial omics.
#'
#' @param geomx_set A \code{NanoStringGeoMxSet} object.
#' @param method A character string specifying the normalization method. 
#'   Must be one of \code{"neg"} (Background/IgG normalization), 
#'   \code{"hk"} (Housekeeping normalization), or \code{"q3"} (75th percentile 
#'   quantile normalization). Defaults to \code{"neg"}.
#'
#' @return A \code{NanoStringGeoMxSet} object with the normalized matrix stored in 
#'   \code{assayDataElement(object, "neg_norm")}, \code{"hk_norm"}, or \code{"q_norm"} 
#'   depending on the chosen method.
#'
#' @details
#' This function checks if the specified normalized matrix already exists. If not, 
#' it calls \code{NanoStringNCTools::normalize} using the user-specified method. 
#' For spatial protein data, \code{"neg"} is generally recommended to account for 
#' non-specific antibody binding. 
#' 
#' Note: For \code{"neg"} and \code{"hk"} to work, the \code{featureData} slot of 
#' the object must have appropriate controls defined, and \code{featureType(object)} 
#' must be set to \code{"Target"}.
#'
#' @importFrom NanoStringNCTools normalize
#' @importFrom Biobase assayData
#' @export
#'
#' @examples
#' \dontrun{
#'   # Background (IgG) normalization (Default)
#'   geomx_set <- NormalizeCounts(geomx_set)
#'   
#'   # Housekeeping normalization
#'   geomx_set <- NormalizeCounts(geomx_set, method = "hk")
#'   
#'   # Q3 normalization
#'   geomx_set <- NormalizeCounts(geomx_set, method = "q3")
#' }
NormalizeCounts <- function(geomx_set, method = c("neg", "hk", "q3")) {
  
  # match.arg ensures the user inputs a valid option and sets the first item as default
  method <- match.arg(method)
  
  require(NanoStringNCTools)
  require(GeomxTools)
  
  # Configure variables based on the selected method
  if (method == "q3") {
    norm_method_str <- "quant"
    to_elt <- "q_norm"
    msg_str <- "Q3 (75th percentile)"
  } else if (method == "hk") {
    norm_method_str <- "hk"
    to_elt <- "hk_norm"
    msg_str <- "Housekeeping (HK)"
  } else if (method == "neg") {
    norm_method_str <- "neg"
    to_elt <- "neg_norm"
    msg_str <- "Background/IgG (Neg)"
  }
  
  message(sprintf("Performing %s Normalization...", msg_str))
  
  # Check if it already exists to avoid redundant computation
  if (to_elt %in% names(Biobase::assayData(geomx_set))) {
    message(sprintf("Note: '%s' already exists. Overwriting...", to_elt))
  }
  
  # Run standard normalization within a tryCatch block
  tryCatch({
    if (method == "q3") {
      geomx_set <- NanoStringNCTools::normalize(
        geomx_set,
        norm_method = norm_method_str, 
        desiredQuantile = 0.75,
        toElt = to_elt
      )
    } else {
      # "hk" and "neg" don't take the desiredQuantile argument
      geomx_set <- NanoStringNCTools::normalize(
        geomx_set,
        norm_method = norm_method_str,
        toElt = to_elt
      )
    }
  }, error = function(e) {
    stop("Normalization failed. Please check input object structure and ensure controls are defined.\nError: ", e$message)
  })
  
  message(sprintf("Success! '%s' added to assayData.", to_elt))
  return(geomx_set)
}