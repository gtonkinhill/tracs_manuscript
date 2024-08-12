library(tidyverse)


# lambda <- 10
# trans_snp_time <- rep(1, 1000) #rexp(100, rate = 2)*10
# trans_snp_dist <- map_dbl(trans_snp_time, ~ rpois(n = 1, lambda = .x*lambda))
# 
# unrelated_snp_dist <- rnbinom(1000, mu = 1000, size = 1)
# 
# mix_snp_dist <- map_dfr(1:1000, ~{
#   
#   tt <- rexp(1, rate = 2)*10
#   if (runif(1)<0.2){
#     dd <- rpois(n = 1, lambda = tt*lambda)
#   } else {
#     dd <- rnbinom(1, mu = 1000, size = 1)
#   }
#   tibble(snp_dist=dd, time_dist=tt)
# 
# })
# 
# res <- calculate_snp_cutoff(mix_snp_dist$snp_dist, unrelated_snp_dist)

log_sum_exp <- function(log_a, log_b) {
  # Ensure log_a is the max
  if (log_a < log_b) {
    tmp <- log_a
    log_a <- log_b
    log_b <- tmp
  }
  # Return the sum in log space
  return(log_a + log(1 + exp(log_b - log_a)))
}



mixture_snp_cutoff <- function(trans_snp_dist, unrelated_snp_dist, 
                                 upper.tail=0.95, max_false_positive=0.05){
  
  if ((length(trans_snp_dist) >= 20) && (length(unrelated_snp_dist) >= 20)){
    nb_fit <- MASS::fitdistr(x=unrelated_snp_dist, densfun = "negative binomial")
    
    llk <- function(params, x){
      k <- params[[1]]
      lambda <- params[[2]]
      
      -sum(map_dbl(x, ~ {log_sum_exp(log(k) + dpois(x = .x, 
                                                        lambda =  lambda, 
                                                        log = TRUE),
                                         log(1-k) + dnbinom(x = .x, 
                                                            size = nb_fit$estimate['size'], 
                                                            mu = nb_fit$estimate['mu'], 
                                                            log = TRUE))}))
    }
    
    start_params <- c(0.5, 1) # Initial guesses
    result <- optim(par = start_params, fn = llk, x = trans_snp_dist, 
                    method = "L-BFGS-B", lower = c(0, 1e-10), upper = c(1, Inf))
    
    snp_threshold <- qpois(upper.tail, lambda = result$par[[2]])
  } else {
    warning("Insufficient data points to fit distributions!")
    return(
      tibble(
        snp_threshold=NA,
        lambda=NA,
        k=NA,
        estimated_fp=NA
      )
    )
  }
  
  if ((sum(unrelated_snp_dist<snp_threshold)/length(unrelated_snp_dist)) > max_false_positive){
    warning(paste0("Inferred SNP threshold may have a false positive rate above ",
                              max_false_positive, "!"))
  }
  
  estimated_fp <- pnbinom(q = snp_threshold,  
                          size = nb_fit$estimate['size'], 
                          mu = nb_fit$estimate['mu'], 
                          lower.tail = TRUE)
  
  return(
    tibble(
      snp_threshold=snp_threshold,
      lambda=result$par[[2]],
      k=result$par[[1]],
      estimated_fp=estimated_fp
    )
  )
  
}


mixture_snp_cutoff_time <- function(trans_snp_dist, trans_time_dist, unrelated_snp_dist, 
                                 max_time=NA, max_false_positive=0.05){
  
  if ((length(trans_snp_dist) >= 30) && (length(unrelated_snp_dist) >= 30)){
    nb_fit <- MASS::fitdistr(x=unrelated_snp_dist, densfun = "negative binomial")
    
    llk <- function(params, x, t){
      k <- params[[1]]
      lambda <- params[[2]]
      
      -sum(map2_dbl(x, t, ~ {log_sum_exp(log(k) + dpois(x = .x, 
                                                        lambda =  lambda*.y, 
                                                        log = TRUE),
                                         log(1-k) + dnbinom(x = .x, 
                                                            size = nb_fit$estimate['size'], 
                                                            mu = nb_fit$estimate['mu'], 
                                                            log = TRUE))}))
    }
    
    start_params <- c(0.5, 1) # Initial guesses
    result <- optim(par = start_params, fn = llk, x = trans_snp_dist, t=trans_time_dist, 
                    method = "L-BFGS-B", lower = c(0, 1e-10), upper = c(1, Inf))
    
    if (is.na(max_time)) {
      max_time <- 2*max(trans_time_dist)
    }
    
    snp_threshold <- result$par[[2]]*max_time
  } else {
    warning("Insufficient data points to fit distributions!")
    return(
      tibble(
        snp_threshold=NA,
        lambda=NA,
        k=NA,
      )
    )
  }
  
  if ((sum(unrelated_snp_dist<snp_threshold)/length(unrelated_snp_dist)) > max_false_positive){
    warning(paste0("Inferred SNP threshold may have a false positive rate above ",
                              max_false_positive, "!"))
  }
  
  return(
    tibble(
      snp_threshold=snp_threshold,
      lambda=result$par[[2]],
      k=result$par[[1]],
    )
  )
  
}


youden_snp_cutoff <- function(trans_snp_dist, unrelated_snp_dist, 
                              max_false_positive=0.05){
  
  if ((length(trans_snp_dist) >= 10) && (length(unrelated_snp_dist) >= 10)){
  
  youden_index <- map_dbl(trans_snp_dist, ~{
    tp <- sum(trans_snp_dist <= .x)
    tn <- sum(unrelated_snp_dist > .x)
    fn <- sum(trans_snp_dist > .x)
    fp <- sum(unrelated_snp_dist <= .x)
    
    return(tp/(tp+fn) + tn/(tn+fp) - 1)
  })
  youden_snp_threshold <- trans_snp_dist[which.max(youden_index)]
  
  return(
    tibble(
      snp_threshold=youden_snp_threshold,
      J=max(youden_index),
      estimated_fp=sum(unrelated_snp_dist<youden_snp_threshold)/length(unrelated_snp_dist)
    )
  )
  
  } else {
    warning("Insufficient data points to call cut-point!")
    return(
      tibble(
        snp_threshold=NA,
        J=NA,
        estimated_fp=NA
      )
    )
  }
}
