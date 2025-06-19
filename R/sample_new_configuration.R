#' Sample new configuration, Stochastic Shotgun Search (SSS)
#'
#' Given the fully evaluated neighborhood of configurations in a iteration in SSS,
#' this function samples the configuration that will be used as the initial
#' configuration in the next SSS iteration.
#'
#' @param init_config Numeric vector of the current initial configuration.
#' @param init_configs String vector of all initial configurations, collapsed into strings.
#' @param prob Probability vector for sampling the next initial configuration.
#' @param n_causals Number of causal configurations in init_config, 'length(init_config)'
#' @param p Number of variants.
#' @param non_causals Set of variants not in 'init_config'.
#'
#' @return Initial configuration for the next iteration of SSS.
#'
.sample_new_config <- function(init_config,
                               init_configs,
                               prob,
                               n_causals,
                               p,
                               non_causals){
  while(paste(init_config, collapse = ",") %in% init_configs ){
    if(all(prob == 0)){
      break
    }
    # Saving init_config
    init_config_save <- init_config

    # Sampling new initial config.
    new_ind <- sample(1:(p + n_causals*p - n_causals^2), size = 1, prob = prob)


    # Note on configs 1:
    # After init_config is sampled, for the first p configurations
    #   evaluated in the neighborhood a SNP is added if it is not in the
    #   init_config, otherwise it is removed.
    if(new_ind %in% init_config){
      init_config <- init_config[-which(init_config == new_ind)]
    }
    if((!(new_ind %in% init_config))&(new_ind < p+1)){
      init_config <- sort(c(init_config, new_ind))
    }
    # Note on configs 2:
    # For the remaining configs, one at a time, each SNP in the initial.config
    #   is replaced with a SNP not currently in the initial config.
    if(p < new_ind){
      causal_to_remove <- floor((new_ind - p - 1)/(p - n_causals)) + 1
      noncausal_to_add <- 1 + (new_ind - p - 1) %% (p - n_causals)
      init_config <- sort(c(init_config[-causal_to_remove], non_causals[noncausal_to_add]))
    }

    if(paste(init_config, collapse = ",") %in% init_configs ){
      prob[new_ind] <- 0
      init_config <- init_config_save
    }
  }
  return(init_config)
}
