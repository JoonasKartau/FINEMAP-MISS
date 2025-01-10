#' Fine-mapping meta-analyses with missing information
#'
#' \code{FINEMAPMISS} performs Bayesian variable selection on a set of genetic
#' variants and a phenotype, using marginal effect estimates \eqn{\boldsymbol{\widehat \beta}},
#'  their standard errors \eqn{\boldsymbol{s}}, and  a reference LD panel \eqn{\boldsymbol{R}}
#'  as input. \cr
#'  \code{FINEMAPMISS} differs from previous methods, as it allows for the combination
#'  of multiple datasets across biobanks, while also functioning with missing data. \cr
#'  If given separate datasets, this function will automatically perform the meta-analysis
#'  and fine-map the data. \cr
#'  If the data have been combined into a single dataset in advance, \code{FINEMAPMISS}
#'  can still be run, but with reduced accuracy.
#'
#' @param ses A vector or matrix of GWAS marginal effect standard errors.
#' Takes the form of a vector if there is only a single dataset or the data has
#' been combined in advance. Otherwise, it is a matrix where each column contains
#' the standard errors for one dataset. For any unobserved variants in a study,
#' the standard errors should be set to \code{Inf}.
#' @param betas A vector or matrix of GWAS marginal effects.
#' Takes the form of a vector if there is only a single dataset or the data has
#' been combined in advance. Otherwise, it is a matrix where each column contains
#' the marginal effects for one dataset. For any unobserved variants in a study,
#' the marginal effects should be set to \code{0}.
#' @param INFO A vector or matrix of variant INFO scores.
#' Takes the form of a vector if there is only a single dataset or the data has
#' been combined in advance. Otherwise, it is a matrix where each column contains
#' the INFO scores for one dataset. For any unobserved variants in a study,
#' the marginal effects should be set to \code{0}.
#' @param R Reference LD matrix for the set of analyzed variants.
#' @param tau Prior standard error for the casual effects. Set by default to \code{0.05}.
#' @param adj Adjustment constant \eqn{\varepsilon}, shrinks the correlations in
#'  \eqn{\boldsymbol{R_M}} to make its inversion stable.
#' @param n_reps Maximum number of Stochastic Shotgun Search (SSS) iterations.
#' @param prob_threshold SSS termination threshold, based on probability mass added at each iteration.
#' @param max_causals Maximum number of causal variants.
#' @param cred_sizes Vector or numeric of size of credible sets to be evaluated. If no value
#' is provided, the rounded expected posterior for number of causal variants is used.
#' @param cred_eval Binary variable, should credible sets be computed?
#' @param meta_analyze Binary variable, should the data be meta-analyzed?
#' If input data (\code{ses, betas, freqs, INFO}) is in the form of a matrix, then
#' \code{meta_analyze} needs to be set to \code{TRUE}.
#' @param quiet Should the current progress of SSS be repressed?
#' @param supplied_matrices If the vector \eqn{\boldsymbol{\widehat z}^{T} \boldsymbol{R_M}^{-1}\boldsymbol{R}_{\boldsymbol{\gamma}} }
#' and matrix \eqn{\mathbb{I}_p + \boldsymbol{R}_{\boldsymbol{\gamma}}^{T} \boldsymbol{R_M}^{-1} \boldsymbol{R}_{\boldsymbol{\gamma}} }
#' have been computed in advance, they can be input as a list to avoid computing them again.
#' @param adj_match Should the \eqn{\boldsymbol{R}} be adjusted in a similar way to \eqn{\boldsymbol{R_M}} when computing log Bayes-factors?
#' @param RM If \eqn{\boldsymbol{R_M}} has been computed in advance, it can be provided as input.
#' @param estimate_M Should \eqn{\boldsymbol{M}} be estimated? This should be used if data has been meta-analyzed in advance.
#' @param beta_shrink Should beta shrinkage in the
#' marginal effect covariances be applied?
#' @param n_studies Number of studies. Should be equal to the number of columns in
#' ses abnd betas. If the data has been combined in advance, set to \code{1}.
#' @param variant_sample_sizes Vector of combined sample sizes per variant, from the meta-analysis.
#' @param max_overlap Binary parameter denoting whether the maximum sample
#' overlap is assumed when estimating \eqn{\boldsymbol{M}}.
#' @param freqs A vector or matrix of variant frequencies.
#' Takes the form of a vector if there is only a single dataset or the data has
#' been combined in advance. Otherwise, it is a matrix where each column contains
#' the variant frequencies for one dataset. For any unobserved variants in a study,
#' the marginal effects should be set to \code{0}.
#' @param scaled_data Has the data been scaled with allele frequencies in advance?
#' @param use_N Are the variant sample sizes or marginal effect standard errors used
#' for fine-mapping?
#'
#' @return A list of objects
#' \describe{
#'   \item{\code{pips}}{Numeric vector of posterior inclusion probabilities per variant.}
#'   \item{\code{stored_configs}}{A matrix of evaluated configurations and associated information.
#'    Columns: 1. causal configuration, 2. log Bayes-factor, 3. size of configuration, 4. Is this a unique configuration?}
#'   \item{\code{cred}}{List or matrix of credible sets.}
#'   \item{\code{post_k}}{Posterior probability distribution for the number of causal variants from \code{0:max_causals}}
#' }

#' @export
#'
#' @examples
FINEMAPMISS <- function(ses, betas, R, tau = 0.05, adj = 0.0001,
                             n_reps = 50,
                             prob_threshold = 0.001, max_causals = 5,
                             cred_sizes = NULL, cred_eval, meta_analyze = T,
                             quiet = T, supplied_matrices = NULL, adj_match = T,
                             RM = NULL,
                             estimate_M = FALSE,
                             beta_shrink = TRUE, n_studies,
                             variant_sample_sizes = NULL,
                             max_overlap = T, freqs,
                             scaled_data = F, use_N = F, INFO = NULL){


  # Checking for suitable input.
  ses <- as.matrix(ses)
  betas <- as.matrix(betas)

  if(!all(dim(ses) == dim(betas)) | !all(dim(ses) == dim(freqs))){
    stop("Error: dimensions of ses, betas, and freqs must be equal.")
  }
  if(!is.numeric(n_studies)){
    stop("Error: non-numeric input for n_studies")
  }

  if(dim(ses)[2] != n_studies | dim(betas)[2] != n_studies | dim(freqs)[2] != n_studies){
    stop("Error: Number of columns of ses, betas, freqs do not match n_studies parameter.")
  }

  if(dim(R)[1] != dim(R)[2]){
    stop("Error: R is not a square matrix")
  }
  if(dim(R)[1] != dim(ses)[1]){
    stop("Error: dimensions of R do not match data, (ses, betas, freqs)")
  }
  if(!is.numeric(R)){
    stop("Error: non-numeric input for R")
  }

  if(!is.null(RM)){
    if(dim(RM)[1] != dim(RM)[2]){
      stop("Error: RM is not a square matrix")
    }
    if(dim(RM)[1] != dim(ses)[1]){
      stop("Error: dimensions of RM do not match data, (ses, betas, freqs)")
    }
    if(!is.numeric(RM)){
      stop("Error: non-numeric input for RM")
    }
  }
  if(!is.numeric(betas)){
    stop("Error: non-numeric input for betas")
  }
  if(!is.numeric(ses)){
    stop("Error: non-numeric input for ses")
  }
  if(!is.numeric(freqs)){
    stop("Error: non-numeric input for freqs")
  }
  if(!is.numeric(tau)){
    stop("Error: non-numeric input for tau")
  }
  if(!is.numeric(adj)){
    stop("Error: non-numeric input for adj")
  }
  if(!is.numeric(n_reps)){
    stop("Error: non-numeric input for freqs")
  }
  if(!is.numeric(max_causals)){
    stop("Error: non-numeric input for freqs")
  }
  if(!is.numeric(prob_threshold)){
    stop("Error: non-numeric input for prob_threshold")
  }
  if(!is.null(cred_sizes)){
    if(!is.numeric(cred_sizes)){
      stop("Error: non-numeric input for cred_sizes")
    }
  }
  if(!is.null(variant_sample_sizes)){
    if(!is.numeric(variant_sample_sizes)){
      stop("Error: non-numeric input for variant_sample_sizes")
    }
  }



  #Scaling data
  if(scaled_data == F){
    betas <- betas*sqrt(2*freqs*(1-freqs))
    ses[which(ses != Inf)] <- (ses*sqrt(2*freqs*(1-freqs)))[which(ses != Inf)]
  }


  # Storing original vector of marginal effects.
  if(meta_analyze == T){
    meta_analysis <- IVW(betas = betas, ses = ses)
    beta_meta <- meta_analysis[[1]]
    se_meta <- meta_analysis[[2]]

    INFO_meta <- IVW(betas = INFO, ses = ses)[[1]]
  } else {
    beta_meta <- betas
    se_meta <- ses
  }


  # Creating sample size overlap matrix M
  p <- length(beta_meta)
  if(is.null(INFO)){
    INFO_meta <- rep(1,p)
  }

  print("Creating sample overlap matrix")

  if(is.null(RM)){
    RM <- .create_RM_matrix(ses = ses, betas = betas, R = R, beta_shrink = beta_shrink,
                   n_studies = n_studies, estimate_M = estimate_M, p = p,
                   max_overlap = max_overlap)
  }

  # Inverting RM with a diagonal adjustment.
  print("Inverting covariance matrix")
  if(is.null(supplied_matrices)){
    RMi <- solve((RM*(1 - adj) + diag(adj, p)))
    # Making the prior variance of a causal variant a vector of length p.
    tau_vec <- rep(tau[1], p)/se_meta
    z <- beta_meta/se_meta

    # Pre-multiplying matrices. We compute these only once and use specific
    #   cols and rows as needed, instead of computing them individually
    #   for each computation.
    print("Pre-multiplication matrices for fine-mapping")

    if(use_N == T){
      N <- sqrt(variant_sample_sizes/(1 - beta_meta^2))
    } else {
      N <- 1/se_meta/sqrt(INFO_meta)
    }
    if("adj_match" == T){
      RMi_R <-  tau[1]*RMi %*%  diag(N) %*% (R*(1 - adj) + diag(adj,p)) %*% diag(sqrt(INFO_meta))
      tR_RMi_R <-  tau[1]*( t(diag(N) %*% (R*(1 - adj) + diag(adj,p)) %*% diag(sqrt(INFO_meta)))  %*% RMi_R )
    } else {
      RMi_R <-  tau[1]*RMi %*% diag(N) %*% R %*% diag(sqrt(INFO_meta))
      tR_RMi_R <-  tau[1]*( t( diag(N) %*% R %*% diag(sqrt(INFO_meta))) %*% RMi_R )
    }
    z_RMi_R <- t(z %*% RMi_R)
    I_tR_RMi_R <- diag(1, p) + tR_RMi_R
  } else {
    z_RMi_R <- supplied_matrices[[1]]
    I_tR_RMi_R <- supplied_matrices[[2]]
  }


  # Comparison LD matrix. If a configuration has two SNPs with correlation
  #   above a threshold, then it is ignored.
  comp_R <- abs(R)
  diag(comp_R) <- 0

  # Setting an arbitrary initial configuration.
  init_config <- sample(1:p, 1)
  config <- init_config

  # Creating matrix to store configuration and their evaluated scores/LogBFs.
  print("Creating storage")
  max_stored_configs <- 1e6
  stored_configs <- stored_log_bf <- stored_unique_configs <- stored_config_sizes <- vector(length = 1e6)

  # This value keeps track of the number of evaluated configurations.
  tt <- 1

  #Creating log prior vector for number of causal variants k, (ranges from 0:p)
  log_prior <- -(0:p)*log(p) + (p-0:p)*log(1 - 1/p)
  log_prior[(max_causals + 2):(p+1)] <- -Inf
  log_prior[1] <- -Inf
  log_prior <- log_prior - logsumexp(log_prior)


  # Vector of initial configurations for SSS.
  init_configs <- vector()

  # Vector storing the size of each configuration.
  stored_init_config_sizes <- vector()

  #Newly discovered probability "mass" discovered per fine-mapping rep.
  new_mass <- 0

  print(c("Finemapping"))
  for(rep in 1:n_reps){
    if(quiet == F){
      print(paste(round(rep/n_reps*100), "% complete, config size: ", length(init_config), ", new mass: ", min(round(new_mass, 3), 1),", number of configs evaluated: ", tt, sep = ""))
    }

    # Storing the initial config.
    init_configs[length(init_configs) + 1] <- paste(init_config, collapse = ",")

    # Storing the size of the config
    n_causals <- length(init_config)
    stored_init_config_sizes[length(stored_init_config_sizes) + 1] <- n_causals

    # Which SNPs are not in the current initial config.
    non_causals <- setdiff(1:p, as.numeric(init_config))

    # Creating vectors to save log_bf values and configs sizes.
    log_bf <-vector()
    config.causals <- vector()

    # For loop. that goes over each config in the neighborhood of the initial config.
    log_bf <- vector()
    configs <- vector()
    config_sizes <- vector()


    for(kk in 1:(p + n_causals*p - n_causals^2)){
      if(kk < p+1){
        if(kk %in% init_config){
          config <- init_config[-which(init_config == kk)]
        } else {
          config <- sort.int(c(init_config, kk), method = "radix")
        }
      } else {
        causal_to_remove <- floor((kk - p - 1)/(p - n_causals)) + 1
        noncausal_to_add <- 1 + (kk - p - 1) %% (p - n_causals)
        config <- sort.int(c(init_config[-causal_to_remove], non_causals[noncausal_to_add]), method = "radix")
      }

      configs[kk] <- paste(config, collapse = ",")
      config_sizes[kk] <- length(config)
      if(any(comp_R[config,config] > 0.95)|(config_sizes[kk] == 0)){
        # The configuration is skipped if it contains SNPs in high LD.
        log_bf[kk] <- -Inf
      } else {
        # Here we evaluate the config if it does not contain highly correlated
        #   SNPs and is not the null config.
        #   The log prior is also added.
        log_bf[kk] <- .eval_logbf(z_RMi_R = z_RMi_R, I_tR_RMi_R = I_tR_RMi_R, configuration = config) + log_prior[length(config) + 1]

      }
    }
    unique_configs<- !(configs %in% stored_configs)


    if(tt + (p + n_causals*p - n_causals^2) < max_stored_configs){
      #storing data if it fits.
      stored_configs[tt:(tt + (p + n_causals*p - n_causals^2) - 1)] <- configs
      stored_log_bf[tt:(tt + (p + n_causals*p - n_causals^2) - 1)] <- log_bf
      stored_unique_configs[tt:(tt + (p + n_causals*p - n_causals^2) - 1)] <- unique_configs
      stored_config_sizes[tt:(tt + (p + n_causals*p - n_causals^2) - 1)] <- config_sizes
    } else {
      # If the newly evaluated configs do no fit in storage, then it is extended
      stored_configs <- c(stored_configs, vector(length = 1e6))
      stored_log_bf <- c(stored_log_bf, vector(length = 1e6))
      stored_unique_configs <- c(stored_unique_configs, vector(length = 1e6))
      stored_config_sizes <- c(stored_config_sizes, vector(length = 1e6))
      max_stored_configs <- max_stored_configs + 1e6

      stored_configs[tt:(tt + (p + n_causals*p - n_causals^2) - 1)] <- configs
      stored_log_bf[tt:(tt + (p + n_causals*p - n_causals^2) - 1)] <- log_bf
      stored_unique_configs[tt:(tt + (p + n_causals*p - n_causals^2) - 1)] <- unique_configs
      stored_config_sizes[tt:(tt + (p + n_causals*p - n_causals^2) - 1)] <- config_sizes
    }
    tt <- tt + (p + n_causals*p - n_causals^2)

    # Probability vector for sampling next initial config.
    prob <- rep(0, (p + n_causals*p - n_causals^2))
    prob <- exp(log_bf  - logsumexp(log_bf))


    new_mass <- sum(exp(log_bf[unique_configs] - logsumexp(stored_log_bf[stored_unique_configs])))
    if(new_mass < prob_threshold){
      break
    }

    # Whileloop here until a new initial config is selected.
    init_config <- .sample_new_config(init_config = init_config,
                                      init_configs = init_configs,
                                      prob = prob,
                                      n_causals = n_causals,
                                      p = p,
                                      non_causals = non_causals)
  }

  print("Calculating pips")
  # Creating vector of pips.

  lse <- logsumexp(stored_log_bf[stored_unique_configs])


  pips <- .calculate_pips(init_configs = init_configs,
                         p = p,
                         stored_init_config_sizes = stored_init_config_sizes,
                         stored_log_bf = stored_log_bf,
                         stored_unique_configs = stored_unique_configs,
                         lse = lse)


  # Computing posterior probability for number of causal variants 'k'
  post_k <- vector()
  for(ii in 1:max_causals){
    k_ind <- which((stored_config_sizes == ii)&(stored_unique_configs == T))
    post_k[ii] <- sum((exp(stored_log_bf[k_ind] - lse)))
  }

  # Evaluating credible sets
  if(cred_eval == T){
    print("Credible Sets")
    cred <- .create_credible_sets(cred_sizes = cred_sizes,
                                 post_k = post_k,
                                 max_causals = max_causals,
                                 stored_log_bf = stored_log_bf,
                                 stored_config_sizes = stored_config_sizes,
                                 stored_configs = stored_configs,
                                 log_prior = log_prior,
                                 z_RMi_R = z_RMi_R,
                                 I_tR_RMi_R = I_tR_RMi_R,
                                 comp_R = comp_R,
                                 p = p)
  } else {
    cred <- NULL
  }


  # Returning grouped pips, individual pips, snp groups and stored_configs.
    return(list("pips" = pips,
                "stored_configs" = cbind(stored_configs, stored_log_bf, stored_config_sizes, stored_unique_configs),
                "cred" = cred,
                "post_k" = post_k))
}

