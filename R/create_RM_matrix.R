#' Create z-score covariance matrix
#'
#' Computes the covariance matrix \eqn{\boldsymbol{R_M}} for the z-scores from a GWAS meta-analysis
#' with possible missing information.
#'
#' @param ses A vector or matrix of GWAS marginal effect standard errors.
#' Takes the form of a vector if there is only a single dataset or the data has
#' been combined in advance. Otherwise, it is a matrix where each column contains
#' the standard errors for one dataset. For any unobserved variants in a study,
#' the standard errors should be set to Inf.
#' @param betas A vector or matrix of GWAS marginal effects.
#' Takes the form of a vector if there is only a single dataset or the data has
#' been combined in advance. Otherwise, it is a matrix where each column contains
#' the marginal effects for one dataset. For any unobserved variants in a study,
#' the marginal effects should be set to 0.
#' @param INFO A vector or matrix of variant imputation scores.
#' Takes the form of a vector if there is only a single dataset or the data has
#' been combined in advance. Otherwise, it is a matrix where each column contains
#' the marginal effects for one dataset. For any unobserved variants in a study,
#' the marginal effects should be set to 0.
#' @param R Reference LD matrix for the set of analyzed variants.
#' @param beta_shrink Binary parameter denoting whether beta shrinkage in the
#' marginal effect covariances is applied.
#' @param n_studies Number of studies. Should be equal to the number of columns in
#' ses abnd betas. If the data has been combined in advance, set to 1.
#' @param p Number of variants.
#' @param estimate_M Binary parameter denoting whether the sample overlap matrix
#' \eqn{\boldsymbol{M}} should be estimated rather than directly computed. Should
#' be used when the data has been combined in advance.
#' @param max_overlap Binary parameter denoting whether the maximum sample
#' overlap is assumed when estimating \eqn{\boldsymbol{M}}.
#'
#' @return z-score covariance matrix \eqn{\boldsymbol{R_M}}
#' @export
#'
#'
#'
.create_RM_matrix <- function(ses, betas = NULL, INFO, R, beta_shrink = T, n_studies, p,
                              estimate_M = FALSE, max_overlap = T){
  if(is.vector(ses) & is.vector(betas)){
    ses <- matrix(ses, ncol = 1)
    betas <- matrix(betas, ncol = 1)
  }


  w <- ses^(-1)
  w <- w^2
  W_list <- list()
  beta_list <- list()
  M <- matrix(0, p, p)


  if(estimate_M == FALSE){
    if(is.list(R)){
      for(ii in 1:n_studies){
        W_list[[ii]] <- sqrt(w[,ii]) %*% t(sqrt(w[,ii]))
        beta_list[[ii]] <- (1-cbind(betas[,ii]^2*INFO[,ii], rep(1,length(betas[,ii]))) %*% t(cbind(rep(1,length(betas[,ii])), betas[,ii]^2*INFO[,ii])) + (betas[,ii]*sqrt(INFO[,ii])) %*% t(betas[,ii]*sqrt(INFO[,ii]))*R[[ii]])/(sqrt(1 - betas[,ii]^2*INFO[,ii]) %*% t(sqrt(1 - betas[,ii]^2*INFO[,ii])))

      }
      for(ii in 1:n_studies){
        M <- M + W_list[[ii]]*beta_list[[ii]]*R[[ii]]
      }
    }
    if(is.matrix(R)){
      for(ii in 1:n_studies){
        W_list[[ii]] <- sqrt(w[,ii]) %*% t(sqrt(w[,ii]))
        beta_list[[ii]] <- (1-cbind(betas[,ii]^2*INFO[,ii], rep(1,length(betas[,ii]))) %*% t(cbind(rep(1,length(betas[,ii])), betas[,ii]^2*INFO[,ii])) + (betas[,ii]*sqrt(INFO[,ii])) %*% t(betas[,ii]*sqrt(INFO[,ii]))*R)/(sqrt(1 - betas[,ii]^2*INFO[,ii]) %*% t(sqrt(1 - betas[,ii]^2*INFO[,ii])))

      }
      for(ii in 1:n_studies){
        M <- M + W_list[[ii]]*beta_list[[ii]]*R
      }
    }
    return(M/(sqrt(rowSums(w)) %*% t(sqrt(rowSums(w)))))
  }

  if(estimate_M == TRUE){
    if(max_overlap == T){
      N <- sqrt(w) %*% t(sqrt(1/w))
      N <- ifelse(N > 1, 1/N, N)
      N <- N*(1-cbind(betas^2*INFO, rep(1,length(betas))) %*% t(cbind(rep(1,length(betas)), betas^2*INFO)) + ((betas*sqrt(INFO)) %*% t(betas*sqrt(INFO)))*R)/(sqrt(1 - betas^2*INFO) %*% t(sqrt(1 - betas^2*INFO)))
      return(N*R)
    }

    if(max_overlap == F){
      N <- ifelse((cbind(1, w) %*% t(cbind(w, 1)) - max(w)) < 0 , 0,  (cbind(1, w) %*% t(cbind(w, 1)) - max(w)))
      N <- N * (sqrt(1/w) %*% t(sqrt(1/w)))
      diag(N) <- 1
      N <- N*(1-cbind(betas^2*INFO, rep(1,length(betas))) %*% t(cbind(rep(1,length(betas)), betas^2*INFO)) + ((betas*sqrt(INFO)) %*% t(betas*sqrt(INFO)))*R)/(sqrt(1 - betas^2*INFO) %*% t(sqrt(1 - betas^2*INFO)))
      return(N*R)
    }
  }
}
