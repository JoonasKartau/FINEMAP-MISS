#' Inverse-variance weighting
#'
#'Performs a fixed-effects meta-analysis given a set of measurements and standard errors
#'across multiple datasets.
#'
#' @param betas Matrix of measurements. Each column contains the measurements
#'  from one study. Any missing values from a dataset should be set to 0.
#' @param ses Matrix of measurement standard errors. Each column contains the standard errors
#'  from one study. Any missing values from a dataset should be set to Inf.
#'
#' @return List containing the meta-analyzed measurements and standard errors.
#' @export
#'
#' @examples
#' #example meta-analysis of two datasets. Data point 3 is missing in dataset 1.
#' betas <- cbind(c(0.5, 0.1, 0), c(0.4, 0.15, 0.25))
#' ses <- cbind(c(0.11, 0.1, Inf), c(0.12, 0.11, 0.1))
#' meta_analysis <- IVW(betas, ses)
#'
IVW <- function(betas, ses){
  w <- ses^(-1)
  w <- w^2
  beta_meta <- rowSums(w*betas)/rowSums(w)
  se_meta <- rowSums(w)^(-1/2)
  return(list("beta meta" = beta_meta, "se meta" = se_meta))
}

