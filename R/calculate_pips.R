#' FINEMAP-MISS PIPs
#'
#'Computes the postrior inclusion probabilities (PIPs), given the evaluated
#'causal configurations in FINEMAP-MISS.
#'
#' @param init_configs String vector of initial configurations for each iteration of
#' Shotgun Stochastic Search. Stored as a collapsed string.
#' @param p Number of variants.
#' @param stored_init_config_sizes  Numeric vector of sizes of initial configurations for each iteration of
#' Shotgun Stochastic Search.
#' @param stored_log_bf Numeric vector of log Bayes-factors of evaluated configurations.
#' @param stored_unique_configs Binary Vector of denoting uniqueness of evaluated configuration.
#' @param lse logsumexp of evaluated configuration log Bayes-factors
#'
#' @return Numeric vector of PIPs for each variant.
#'
.calculate_pips <- function(init_configs,
                            p,
                            stored_init_config_sizes,
                            stored_log_bf,
                            stored_unique_configs,
                            lse){
  pips <- rep(0,p)

  # Loop for each SNP to compute its pip.
  for(snp in 1:p){

    # Integer tallying our current position in the stored_configs matrix.
    sum_configs <- 0

    # Creating object to store index of configurations containing variant 'snp'.
    config_ind <- NULL

    # Since we have the list of initial configs, and we always evaluate neighbors in a
    #   specific. order, we can construct the index of configs that include any given SNP.
    for(ii in 1:length(init_configs)){
      if(snp %in% as.numeric(unlist(strsplit(init_configs[ii], split = ",")))){
        t <- (p+1):(p + p*stored_init_config_sizes[ii] - stored_init_config_sizes[ii]^2)

        # When the SNP is in the initial config, it appears in all but one of the
        #   following p configs.
        ind_in_first_p_configs <- setdiff(1:p, snp)

        # When the SNP is in the initial config, after the p first configs in the
        #   neighborhood, the SNPs are swapped with ones not in the initial config.
        #   one at a time. This index gives the configs in which the OTHER SNPs
        #   in the initial config are getting swapped.
        ind_in_configs_after_p <- p + which((as.numeric(unlist(strsplit(init_configs[ii], split = ","))))[(floor((t - p - 1)/(p - stored_init_config_sizes[ii])) + 1)] != snp)
        config_ind <- c(config_ind,  sum_configs + c(ind_in_first_p_configs, ind_in_configs_after_p))
      } else {
        # If the SNP is not in the initial config, it appears in only one of the
        #   first p configs of the neighborhood.
        ind_in_first_p_configs <- snp

        # The order of the current SNP in relation to those in the initial config.
        snp_ord <- sum(as.numeric(unlist(strsplit(init_configs[ii], split = ","))) < snp)

        # If the SNP is not in the initial config, after the p first configs in the
        #   in the neighborhood, the current SNP is included when it is swapped for
        #   one of the SNPs in the initial config.
        ind_in_configs_after_p <- p + snp - snp_ord + rep(p-stored_init_config_sizes[ii], stored_init_config_sizes[ii])*(0:(stored_init_config_sizes[ii]-1))
        config_ind <- c(config_ind, sum_configs +  c(ind_in_first_p_configs,  ind_in_configs_after_p))
      }
      # Tallying the number of configs that have been looped over.
      sum_configs <- sum_configs + p + p*stored_init_config_sizes[ii] - stored_init_config_sizes[ii]^2
    }

    # Calculating pip using logsumexp 'lse'.
    pips[snp] <- sum((exp(stored_log_bf[config_ind[stored_unique_configs[config_ind]]] - lse)))
  }
  return(pips)
}
