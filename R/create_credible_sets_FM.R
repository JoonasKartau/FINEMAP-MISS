#' FINEMAP credible sets
#'
#' Computes credible sets from the fine-mapping output of the FINEMAP model.
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
#' @param comp_R LD matrix, with diagonal set to 0.
#' @param RM z-score covariance matrix
#' @param R LD matrix.
#' @param Si_Ii Diagonal matrix of standard error and square root of INFO score inverses.
#' @param I2 Diagonal matrix of INFO scores.
#' @param z vector of z-scores.
#' @param tau prior standard deviation.
#' @param p Number of variants.
#'
#' @return List of credible sets as matrices for each 'cred_size' in 'cred_sizes'. If 'length(cred_sizes) == 1', then the
#' single credible set matrix is returned.
#'
.create_credible_sets_FM <- function(cred_sizes,
                                  post_k,
                                  max_causals,
                                  stored_log_bf,
                                  stored_config_sizes,
                                  stored_configs,
                                  log_prior,
                                  comp_R,
                                  RM,
                                  R,
                                  p,
                                  Si_Ii,
                                  I2,
                                  z,
                                  tau){
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
          log_abf[ii] <- .log_dmvnorm(x = z[ii], S = RM[ii,ii] + tau^2*(Si_Ii[ii,ii] %*% R[ii,ii] %*% I2[ii,ii] %*% R[ii,ii] %*% Si_Ii[ii,ii])) -
            .log_dmvnorm(x = z[ii], S = RM[ii, ii]) + log_prior[2]
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
              log_abf[ii] <- .log_dmvnorm(x = z[config], S = RM[config,config] + tau^2*(Si_Ii[config,config] %*% R[config,config] %*% I2[config,config] %*% R[config, config] %*% Si_Ii[config,config])) -
                .log_dmvnorm(x = z[config], S = RM[config, config]) + log_prior[length(config) + 1]
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
      cred_list[[which(cred_sizes %in% cred_size)]] <- cred
    }
  }
  if(length(cred_sizes) > 1){
    return(cred_list)
  } else {
    return(cred)
  }

}
