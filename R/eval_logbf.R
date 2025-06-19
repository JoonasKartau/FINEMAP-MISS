#' Evaluate log Bayes-factor
#'
#' Evaluate the log Bayes-factor of a configuration in FINEMAP-MISS. Computed using
#' the pre-computed vector \eqn{\boldsymbol{\widehat z}^{T} \boldsymbol{R_M}^{-1}\boldsymbol{R}_{\boldsymbol{\gamma}} }
#' and matrix \eqn{\mathbb{I}_p + \boldsymbol{R}_{\boldsymbol{\gamma}}^{T} \boldsymbol{R_M}^{-1} \boldsymbol{R}_{\boldsymbol{\gamma}} }
#'
#' @param z_RMi_R Pre-multiplied vector
#' @param I_tR_RMi_R Pre-multiplied matrix
#' @param configuration Index vector for the current configuration
#'
#' @return Numeric value of the log Bayes factor for the configuration
#'
.eval_logbf <- function(z_RMi_R, I_tR_RMi_R, configuration){
  if(length(configuration) > 1){
    return(-0.5*(determinant(I_tR_RMi_R[configuration, configuration], logarithm = T)$modulus - t(z_RMi_R[configuration, ]) %*% solve(I_tR_RMi_R[configuration, configuration]) %*% z_RMi_R[configuration, ]))
  } else {
    return(-0.5*(log(abs(I_tR_RMi_R[configuration, configuration])) - t(z_RMi_R[configuration, ]) %*% solve(I_tR_RMi_R[configuration, configuration]) %*% z_RMi_R[configuration, ]))
  }

}
