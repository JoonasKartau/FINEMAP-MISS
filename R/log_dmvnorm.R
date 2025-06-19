#' Log-dmvnorm
#' evaluate the log-density of a multivariate normal distribution.
#'
#' @param x vector, point to be evaluated.
#' @param mu vector, distribution mean
#' @param S matrix, distribution covariance
#' @param pivot should pivoting be used for cholseky decomposition?
#' @param cholesky is cholesky decomposition used for inversion? If not base R solve is used.
#'
#' @return log-density of the observed point x.
#'
.log_dmvnorm <- function(x,
                        mu = rep(0, length(x)),
                        S = diag(1, length(x)),
                        pivot = F,
                        cholesky = T){
  #returns log of density of MV-Normal(mean = mu, var = S) at x
  K = length(mu)
  stopifnot(all(dim(S) == K))
  stopifnot(length(x) == K)
  if(cholesky == T){
    chol.S = chol(S, pivot = pivot) #Cholesky decomposition
    log.det <- 2*sum(log(diag(chol.S)))
    return(-K/2*log(2*pi) - 0.5*(log.det + t(x - mu) %*% chol2inv(chol.S) %*% (x-mu)))
  } else {
      if(is.null(dim(S))){
        return(-K/2*log(2*pi) - 0.5*(t(x - mu) %*% (x-mu)))
      } else {
        return(-K/2*log(2*pi) - 0.5*(log(abs(det(S))) + t(x - mu) %*% solve(S) %*% (x-mu)))
      }
  }

}

