#' Impute summary statistics
#'
#' @param R LD matrix
#' @param z z-scores
#' @param observed Index vector of observed variants
#' @param unobserved Index vector of unobserved variants
#' @param shrink Shrinkage factor
#' @param return_z Return results as z-scores?
#' @param scale Scale z-scores after imputation?
#' @param n Sample size
#'
#' @return Matrix with columns: 1. Imputed z-scores (marginal effect if return_z == F), 2. Imputation quality, 3. (if return_z == F) standard error.
#' @export
#'
#' @examples
#'#Loading toy data with one true causal variant.
#'data("toydata_FINEMAPMISS")
#'
#'betas <- toydata_FINEMAPMISS$betas
#'ses <- toydata_FINEMAPMISS$ses
#'MAF <- toydata_FINEMAPMISS$MAF
#'LD <- toydata_FINEMAPMISS$LD
#'
#'n <- toydata_FINEMAPMISS$study_sample_sizes
#'p <- dim(LD)[1]
#'
#'#Simulating missingness in 20% of variants (including the true causal variant)
#'missing_data <- unique(sort(c(toydata_FINEMAPMISS$causal_snp,
#'                               cbind(sample(1:p, round(p*0.2))))))
#'
#'#Which dataset are the variants missing from?
#'which_dataset_missing <- sample(1:2, length(missing_data), replace = TRUE)
#'
#'#Setting missing data to 0 or Inf.
#'betas[missing_data[which_dataset_missing == 1],1] <- 0
#'ses[missing_data[which_dataset_missing == 1],1] <- Inf
#'MAF[missing_data[which_dataset_missing == 1],1] <- 0
#'
#'#Setting variants sample sizes
#'variant_sample_sizes <- rep(sum(n), p)
#'variant_sample_sizes[missing_data] <- n[1]
#'
#'#Index of observed variants from study1
#'obs <- setdiff(1:p, missing_data[which_dataset_missing == 1])
#'
#'#Index of unobserved variants from study 1
#'unobs <- sort(missing_data[which_dataset_missing == 1])
#'
#'#z-scores from study 1
#'z_obs <- (betas/ses)[obs, 1]
#'
#'imputation <- impute_summary_stats(R = LD, z = z_obs,
#'                                     observed = obs,
#'                                     unobserved = unobs, n = n,
#'                                     return_z = FALSE, scale = TRUE)
#'
#'
impute_summary_stats <- function(R, z, observed, unobserved, shrink = 0.0001, return_z = TRUE, scale = TRUE,
                   n = NULL){

  if((return_z == F)&(!is.numeric(n))){
    warning("n is not a numeric value")
  }
  R_shrunk <- (1-shrink)*R + shrink*diag(dim(R)[1])
  chol.R <- chol(R_shrunk[observed, observed])
  R_inv <- chol2inv(chol.R)

  imputed_z <- R_shrunk[unobserved, observed, drop = FALSE] %*% R_inv %*% z
  imputed_z_var <- diag(R_shrunk[unobserved, observed, drop = FALSE] %*% R_inv %*% R_shrunk[observed, unobserved, drop = FALSE] )

  if(scale == T){
    if(return_z == T){
      return(cbind("z" = imputed_z/sqrt(imputed_z_var), "Quality" = imputed_z_var))
    } else {
      return(cbind("beta" = imputed_z/sqrt(n)/sqrt(imputed_z_var), "Quality" = imputed_z_var, "se" = 1/sqrt(n)))
    }
  } else {
    if(return_z == T){
      return(cbind("z" = imputed_z, "Quality" = imputed_z_var))
    } else {
      return(cbind("beta" = imputed_z/sqrt(n), "Quality" = imputed_z_var, "se" = 1/sqrt(n)))
    }
  }
}
