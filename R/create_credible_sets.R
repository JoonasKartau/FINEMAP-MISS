#' FINEMAP-MISS credible sets
#'
#' Computes credible sets from the fine-mapping output of FINEMAP-MISS.
#'
#'
#' @param cred_sizes Vector or numeric of size of credible sets to be created.
#' If NULL, then the rounded expectation of the posterior number of causal variants
#' is used.
#' @param post_k Vector of posterior probability of number of casual variants,
#' from 0:max_causals
#' @param max_causals Numeric, maximum number of causal variants
#' @param stored_log_bf Vector of log Bayes-factors computed during fine-mapping.
#' @param stored_config_sizes Vector of size of configurations
#' evaluated during fine-mapping.
#' @param stored_configs Vector of evaluated configurations, stored as a collapsed string.
#' @param log_prior Natural log transformed prior for number of causal variants.
#' @param z_RMi_R Pre-multiplied vector for fine-mapping.
#' @param I_tR_RMi_R Pre-multiplied matrix for fine-mapping.
#' @param comp_R LD matrix, with diagonal set to 0.
#' @param p Number of variants.
#'
#' @return List of credible sets as matrices for each 'cred_size' in 'cred_sizes'. If 'length(cred_sizes) == 1', then the
#' single credible set matrix is returned.
#'
.create_credible_sets <- function(cred_sizes,
                                  post_k,
                                  max_causals,
                                  stored_log_bf,
                                  stored_config_sizes,
                                  stored_configs,
                                  log_prior,
                                  z_RMi_R,
                                  I_tR_RMi_R,
                                  comp_R,
                                  p){
  cred_list <- list()
  # If size of credible sets is not specified, expected posterior number of
  #   causal variants is used.
  if(is.null(cred_sizes)){
    cred_sizes <- round(sum(post_k*1:max_causals))
  }

  for(cred_size in cred_sizes){
    # Finding which configuration of size cred_size has the max log_bf.
    max_ind <- which.max(stored_log_bf[which(stored_config_sizes == cred_size)])
    max_config <- as.numeric(unlist(strsplit((stored_configs)[which(stored_config_sizes == cred_size)][max_ind], split = ",")))

    for(replaced_variant in 1:cred_size){
      # Leave out one SNP from the maximal config to create the conditional SNP set.
      conds <- max_config[-replaced_variant]
      non_conds <- setdiff(1:p, conds)

      #Keeping track of which SNPs are being conditioned on.
      non_cond_ind <- c(1:p)[non_conds]
      target <- max_config[replaced_variant]

      # Vector of log_abfs
      log_abf <- vector()

      #Special case if cred_size == 1
      if(cred_size == 1){
        non_cond_ind <- 1:p
        for(ii in 1:p){
          log_abf[ii] <- .eval_logbf(z_RMi_R = z_RMi_R, I_tR_RMi_R = I_tR_RMi_R, configuration = ii) + log_prior[2]
        }
        lse_abf <- logsumexp(log_abf)
        post <- sapply(1:p, function(x){
          sum((exp(as.numeric(log_abf[x]) - lse_abf)))
        })
      } else {
        # General case when cred_size > 1
        for(ii in non_cond_ind){
          config <- sort(c(ii,conds))
          if(ii %in% conds){
            #Skip any SNP in the condition set.
            log_abf[ii] <- -Inf

          } else {
            if(any(comp_R[ii, c(ii,conds)] > 0.95)){
              # Skip any SNP that is highly correlated with the condition set.
              log_abf[ii] <- -Inf
            } else {
              # Evaluating log_bfs
              log_abf[ii] <- .eval_logbf(z_RMi_R = z_RMi_R, I_tR_RMi_R = I_tR_RMi_R, configuration = config) + log_prior[length(config) + 1]
            }
          }
        }
        # Removing the skipped conditional SNPs.
        log_abf <- log_abf[-conds]

        # Computing posterior probability.
        lse_abf <- logsumexp(log_abf)
        post <- sapply(1:length(non_cond_ind), function(x){
          sum((exp(as.numeric(log_abf[x]) - lse_abf)))
        })
      }

      #Initializing credible set matrix
      if(replaced_variant == 1){
        cred <- matrix(NA, nrow = p, ncol = cred_size*2)
      }

      # Selecting SNPs until the summed probability is greater than the chosen threshold.
      prob <- 0
      post_orig <- post
      ii <- 1
      while(prob < 0.95){
        prob <- prob + post[which.max(post)]
        cred[ii, c(replaced_variant*2-1, replaced_variant*2)] <- c(non_cond_ind[which.max(post)], post[which.max(post)])
        post[which.max(post)] <- 0
        ii <- ii + 1
      }

      # Removing any rows containing only NAs.
      if(replaced_variant == cred_size){
        cred <- cred[1:(min(which(apply(cred, 1, function(x){all(is.na(x))})))-1),]
      }

    }
    if(length(cred_sizes) > 1){
      cred_list[[which(cred_sizes %in% cred_size)]]
    }
  }
  if(length(cred_sizes) > 1){
    return(cred_list)
  } else {
    return(cred)
  }

}
