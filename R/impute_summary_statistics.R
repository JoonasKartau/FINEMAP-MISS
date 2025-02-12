#' Impute summary statistics
#'
#' @param R LD matrix
#' @param z z-scores
#' @param typed Index vector of observed variants
#' @param untyped Index vector of unobserved variants
#' @param shrink Shrinkage factor
#' @param return_z Return results as z-scores?
#' @param scale Scale z-scores after imputation?
#' @param n Sample size
#'
#' @return Matrix with columns: 1. Imputed z-scores (marginal effect if return_z == F), 2. Imputation quality, 3. (if return_z == F) standard error.
#' @export
#'
#' @examples
impute_summary_stats <- function(R, z, typed, untyped, shrink = 0.0001, return_z = TRUE, scale = T,
                   n = NULL){

  if((return_z == F)&(!is.numeric(n))){
    warning("n is not a numeric value")
  }
  R_shrunk <- (1-shrink)*R + shrink*diag(dim(R)[1])
  chol.R <- chol(R_shrunk[typed, typed])
  R_inv <- chol2inv(chol.R)

  imputed_z <- R_shrunk[untyped, typed] %*% R_inv %*% z
  imputed_z_var <- diag(R_shrunk[untyped, typed] %*% R_inv %*% R_shrunk[typed, untyped] )

  if(scale == T){
    if(return_z == T){
      return(cbind("z" = imputed_z/sqrt(imputed_z_var), "Quality" = imputed_z_var))
    } else {
      #return(cbind("beta" = imputed_z/sqrt(imputed_z_var)/sqrt(n), "Quality" = imputed_z_var, "se" = 1/sqrt(n)/sqrt(imputed_z_var)))
      return(cbind("beta" = imputed_z/sqrt(n), "Quality" = imputed_z_var, "se" = sqrt(imputed_z_var)/sqrt(n)))
    }
  } else {
    if(return_z == T){
      return(cbind("z" = imputed_z, "Quality" = imputed_z_var))
    } else {
      return(cbind("z" = imputed_z/sqrt(n), "Quality" = imputed_z_var, "se" = 1/sqrt(n)))
    }
  }
}
