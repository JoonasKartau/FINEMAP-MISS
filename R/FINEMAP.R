#' Fine-mapping meta-analyses with missing information
#'
#' \code{FINEMAP} performs Bayesian variable selection on a set of genetic
#' variants and a phenotype, using marginal effect estimates \eqn{\boldsymbol{\widehat \beta}},
#'  their standard errors \eqn{\boldsymbol{s}}, and  a reference LD panel \eqn{\boldsymbol{R}}
#'  as input. This is a recreation of the C++ implementation by Christian Benner and Matti Pirinen.
#'  Not recommended for meta-analyzed GWAS data.
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
#' @param R Reference LD matrix for the set of analyzed variants.
#' @param tau Prior standard error for the casual effects. Set by default to \code{0.05}.
#' @param n_reps Maximum number of Stochastic Shotgun Search (SSS) iterations.
#' @param prob_threshold SSS termination threshold, based on probability mass added at each iteration.
#' @param max_causals Maximum number of causal variants.
#' @param cred_sizes Vector or numeric of size of credible sets to be evaluated. If no value
#' is provided, the rounded expected posterior for number of causal variants is used.
#' @param cred_eval Binary variable, should credible sets be computed?
#' @param meta_analyze Binary variable, should the data be meta-analyzed? If data has been
#' combined in advance, should be set to \code{FALSE}.
#' @param quiet Should the current progress of SSS be repressed?
#' @param RM If \eqn{\boldsymbol{R_M}} has been computed in advance, it can be provided as input.
#' @param estimate_M Should \eqn{\boldsymbol{M}} be estimated? This should be used if data has been meta-analyzed in advance.
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
#' @param init_config Initial configuration, from which fine-mapping is started.
#' @param scaled_data Has the data been scaled with allele frequencies in advance?
#' @param use_N Are the variant sample sizes or marginal effect standard errors used
#' for fine-mapping?
#' @param INFO Vector of variant imputation INFO scores.
#' @param cholesky Is cholesky decomposition used for matrix inversion? If not base R \code{solve} is used.
#' @param rsid rsid vector (optional identifier).
#' @param chromosome chromosome vector (optional identifier).
#' @param position position vector (optional identifier).
#' @param allele1 reference allele vector (optional identifier).
#' @param allele2 alternate allele vector (optional identifier).
#'
#' @return A list of objects
#' \describe{
#'   \item{\code{summary_table}}{Data frame with entries:
#'   \itemize{
#'   \item{\code{rank}: The order in which the variants appear when sorted by PIP.}
#'   \item{\code{rsid}: Rsid (if supplied, NA otherwise).}
#'   \item{\code{chromosome}: Chromosome (if supplied, NA otherwise).}
#'   \item{\code{allele1}: Allele1 (if supplied, NA otherwise).}
#'   \item{\code{allele2}: Allele2 (if supplied, NA otherwise).}
#'   \item{\code{maf}: Minor allele frequency.}
#'   \item{\code{beta}: Marginal effect.}
#'   \item{\code{se}: Standard error.}
#'   \item{\code{z}: Z-score.}
#'   \item{\code{prob}: Posterior inclusion probability (PIP).}
#'   }}
#'   \item{\code{cred}}{List or matrix of credible sets.}
#'   \item{\code{post_k}}{Posterior probability distribution for the number of causal variants from \code{0:max_causals}}
#'   \item{\code{evaluated_configs}}{Optional output if \code{export_configs == TRUE}. Data frame of evaluated configurations and associated information.
#'   \itemize{
#'    \item{\code{configuration}: Causal configuration.}
#'    \item{\code{log_bf}: Causal configuration log Bayes-factor.}
#'    \item{\code{config_size}: Causal configuration size (how many causal variants?).}
#'    \item{\code{unique_config}: Is this a unique configuration?.}
#'   }}
#'  }

#' @export
#'
#' @examples
#'#Running FINEMAP with imputedd summary statistics
#'
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
#'betas[missing_data[which_dataset_missing == 2],2] <- 0
#'ses[missing_data[which_dataset_missing == 1],1] <- Inf
#'ses[missing_data[which_dataset_missing == 2],2] <- Inf
#'MAF[missing_data[which_dataset_missing == 1],1] <- 0
#'MAF[missing_data[which_dataset_missing == 2],2] <- 0
#'
#'#Setting variants sample sizes
#'variant_sample_sizes <- rep(sum(n), p)
#'variant_sample_sizes[missing_data] <- n[1]
#'
#'#Index of observed variants from each study
#'obs1 <- setdiff(1:p, missing_data[which_dataset_missing == 1])
#'obs2 <- setdiff(1:p, missing_data[which_dataset_missing == 2])
#'
#'#'#Index of unobserved variants from each study
#'unobs1 <- sort(missing_data[which_dataset_missing == 1])
#'unobs2 <- sort(missing_data[which_dataset_missing == 2])
#'
#'#z-scores from each study
#'z_obs1 <- (betas/ses)[obs1, 1]
#'z_obs2 <- (betas/ses)[obs2, 2]
#'
#'imputation1 <- impute_summary_stats(R = LD, z = z_obs1,
#'                                     observed = obs1,
#'                                     unobserved = unobs1, n = n[1],
#'                                     return_z = FALSE, scale = TRUE)
#'
#'
#'imputation2 <- impute_summary_stats(R = LD, z = z_obs2,
#'                                    observed = obs2,
#'                                    unobserved = unobs2, n = n[2],
#'                                    return_z = FALSE, scale = TRUE)
#'
#'#Copying over summary statistics
#'betas_imputed <- betas
#'ses_imputed <- ses
#'
#'#Using imputed values to replace missing observations
#'betas_imputed[unobs1,1] <- imputation1[,1]
#'ses_imputed[unobs1,1] <- imputation1[,3]
#'
#'betas_imputed[unobs2,2] <- imputation2[,1]
#'ses_imputed[unobs2,2] <- imputation2[,3]
#'
#'#Assuming all variants fully observed after imputation
#'variant_sample_sizes <- rep(sum(n),p)
#'
#'#Running FINEMAP
#'output_FM <- original_FINEMAP(ses = ses,
#'                          betas = betas,
#'                          R = LD,
#'                          n_studies = 2,
#'                          variant_sample_sizes = variant_sample_sizes,
#'                          freqs = MAF)
#'
#'
#'
#'
#'
original_FINEMAP <- function(betas,
                             ses,
                             R,
                             tau = 0.05,
                             n_reps = 50,
                             prob_threshold = 0.001,
                             max_causals = 5,
                             cred_sizes = NULL,
                             meta_analyze = T,
                             quiet = T,
                             RM = NULL,
                             estimate_M = FALSE,
                             n_studies,
                             variant_sample_sizes = NULL,
                             max_overlap = T,
                             freqs,
                             scaled_data = F,
                             use_N = F,
                             INFO = NULL,
                             cholesky = T,
                             rsid = NULL,
                             position = NULL,
                             allele1 = NULL,
                             allele2 = NULL,
                             chromosome = NULL,
                             export_configs = FALSE,
                             init_config = NULL){

  #Changing input into matrix form.


  ses <- as.matrix(ses)
  betas <- as.matrix(betas)
  freqs <- as.matrix(freqs)

  p <- dim(ses)[1]

  if(is.null(INFO)){
    INFO <- matrix(1, nrow = p, ncol = n_studies)
  } else {
    INFO <- as.matrix(INFO)
  }

  # Checking for suitable input.
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

  if(estimate_M == T & n_studies > 1){
    stop("Error: estimate_M can only be used if n_studies == 1")
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
    MAF_meta <- IVW(betas = freqs, ses = ses)[[1]]


  } else {
    beta_meta <- betas
    se_meta <- ses
    if(is.null(INFO)){
      INFO_meta <- rep(1,p)
    } else {
      INFO_meta <- INFO
    }
    MAF_meta <- freqs
  }
  z <- beta_meta/se_meta

  # Creating sample size overlap matrix M


  if(use_N == T){
    se_meta <- sqrt((1 - INFO_meta*beta_meta^2)/variant_sample_sizes)
    INFO_meta <- rep(1,p)
  }

  I2 <- diag(INFO_meta)
  Si_Ii <- diag(1/se_meta/sqrt(INFO_meta))


  print("Creating sample overlap matrix")

  if(is.null(RM)){
    RM <- .create_RM_matrix(ses = ses, betas = betas, R = R,
                            n_studies = n_studies, estimate_M = estimate_M, p = p,
                            max_overlap = max_overlap, INFO = INFO)
  }

  # Comparison LD matrix. If a configuration has two SNPs with correlation
  #   above a threshold, then it is ignored.
  comp_R <- abs(R)
  diag(comp_R) <- 0

  # Setting an arbitrary initial configuration if not defined by user.
  if(is.null(init_config)){
    init_config <- sample(1:p, 1)
  }
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
        log_bf[kk] <- .log_dmvnorm(x = z[config], S = RM[config,config] + tau^2*(Si_Ii[config,config] %*% R[config,config] %*% I2[config,config] %*% R[config, config] %*% Si_Ii[config,config])) -
          .log_dmvnorm(x = z[config], S = RM[config, config]) + log_prior[length(config) + 1]

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
  names(post_k) <- c(1:max_causals)
  # Evaluating credible sets

    print("Credible Sets")
    cred <- .create_credible_sets_FM(cred_sizes = cred_sizes,
                                  post_k = post_k,
                                  max_causals = max_causals,
                                  stored_log_bf = stored_log_bf,
                                  stored_config_sizes = stored_config_sizes,
                                  stored_configs = stored_configs,
                                  log_prior = log_prior,
                                  RM = RM,
                                  comp_R = comp_R,
                                  p = p,
                                  Si_Ii = Si_Ii,
                                  I2 = I2,
                                  R = R,
                                  z = z,
                                  tau = tau)



  # Returning grouped pips, individual pips, snp groups and stored_configs.
  summary_table <- data.frame("rank" = order(pips, decreasing = T),
                              "rsid" = if(is.null(rsid)){rep(NA, p)} else {rsid},
                              "chromosome" = if(is.null(chromosome)){rep(NA, p)} else {chromosome},
                              "position" = if(is.null(position)){rep(NA, p)} else {position},
                              "allele1" = if(is.null(allele1)){rep(NA, p)} else {allele1},
                              "allele2" = if(is.null(allele2)){rep(NA, p)} else {allele2},
                              "maf" = if(meta_analyze){MAF_meta} else {freqs},
                              "beta" = beta_meta,
                              "se" = se_meta,
                              "z" = z,
                              "prob" = pips)

  if(!is.matrix(cred)){
    cred <- matrix(cred, ncol = length(cred))
  }

  cred_cols <- ncol(cred)
  colnames(cred) <- 1:ncol(cred)
  for(ii in 1:(cred_cols/2)){
    colnames(cred)[(2*ii - 1):(2*ii)] <- c(paste0("CS", ii, "_variants"), paste0("CS", ii, "_prob"))
  }

  PENC <- sum(post_k*(1:max_causals))


  if(export_configs == T){
    evaluated_configs <- data.frame("configuration" = stored_configs,
                                    "log_bf" = stored_log_bf,
                                    "config_size" = stored_config_sizes,
                                    "unique_config" = stored_unique_configs)
    evaluated_configs <- evaluated_configs[which(stored_unique_configs == TRUE), ]

    return(list("summary_table" = summary_table,
                "credible_sets" = cred,
                "post_prob_n_causal_variants" = post_k,
                "post_expected_n_causal_variants" = PENC,
                "evaluated_configs" = evaluated_configs))
  } else {
    return(list("summary_table" = summary_table,
                "credible_sets" = cred,
                "post_prob_n_causal_variants" = post_k,
                "post_expected_n_causal_variants" = PENC))
  }
}

