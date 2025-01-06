#' logsumexp
#'
#'Compute logsumexp for a vector of input values. Useful for computing
#'probabilities for large values that have been transformed to log scale.
#'
#' @param x Vector of values in natural log scale.
#'
#' @return Numeric value \eqn{\textrm{ln}(\sum_{i = 1}^n e^{x_i})}
#' @export
#'
#' @examples
#' log_bayes_factor <- c(1000, 1000.1, 1000.2) #vector of large values in natural log scale.
#' lse <- logsumexp(x = log_bayes_factor)
#' exp(log_bayes_factor - lse) #provides correct normalized values
#' exp(log_bayes_factor)/sum(exp(log_bayes_factor)) #numeric overflow
logsumexp <- function(x){
  return(max(x) + log(sum(exp(x - max(x)))))
}


