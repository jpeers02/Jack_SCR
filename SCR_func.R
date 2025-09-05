
library(spatstat)
library(secr)
library(acre)
library(tidyverse)
load("C:/Users/Jack/Downloads/test-data.RData")

# add data set to github

traps.1 <- test.data$traps

## Simulating a single session & occasion data set

sim_scr_data <- function(D, lambda0 = NULL, g0 = NULL, 
                         sigma, xlim, ylim, 
                         traps, detfn = c('hn','hhn')) {
  detfn <- match.arg(detfn)
  a_m2 <- diff(xlim) * diff(ylim)
  a_ha <- a_m2 / 10000
  N <- rpois(1, D * a_ha)
  
  activity_centers <- cbind(runif(N, xlim[1], xlim[2]),
                            runif(N, ylim[1], ylim[2]))
  
  dists <- crossdist(activity_centers[, 1], activity_centers[, 2],
                     traps[, 1], traps[, 2])
  if(detfn == "hn"){
    p <- g0 * exp(-dists^2 / (2 * sigma^2))
  } else if(detfn == "hhn"){
    lamd <- lambda0 * exp(-dists^2 / (2 * sigma^2))
    p <- 1 - exp(-lamd)
  }
  
  
  capt <- matrix(rbinom(length(p), 1, p), nrow = N)
  capt_filt <- capt[rowSums(capt) > 0, ]
  list(capt = capt_filt, N = N, activity_centers = activity_centers)
}

## Adapting the data set into  format required for secr and acre

# Adapt to allow for multi-session

sim_adapt_data <- function(D,lambda0 = NULL, g0 = NULL, sigma, xlim, ylim, traps, session_id = 1, detfn) {
  
  sim <- sim_scr_data(D, lambda0, g0, sigma, xlim, ylim, traps, detfn)
  
  capt <- sim$capt
  capthist_df <- data.frame()
  for (i in 1:nrow(capt)) {
    for (j in 1:ncol(capt)) {
      if (capt[i, j] > 0) {
        capthist_df <- rbind(capthist_df, 
                             data.frame(session = session_id,
                                        ID = paste0("ind", i),
                                        occasion = 1,
                                        stringsAsFactors = FALSE, trap = j))
      }
    }
  }
  
  return(list(sim_data = sim, capthist_df = capthist_df))
}

# Start simulating multi-session data - list of traps 

# Buffer = 5 * Sigma

# add hazard-half normal detection function simulation - good for acoustics - allow secr & acre to use.

model.sim <- function(n_sims, D, g0 = NULL, lambda0 = NULL, sigma, 
                      xlim, ylim, traps, buffer, detfn, 
                      method = c('secr', 'acre', 'both')) {
  
  method <- match.arg(method)
  results <- vector("list", n_sims)
  
  for (i in 1:n_sims) {
    
    sim_data <- sim_adapt_data(D = D, g0 = g0, lambda0 = lambda0, sigma = sigma, 
                               xlim = xlim, ylim = ylim, traps = traps, 
                               session_id = 1, detfn = detfn)
    
    sim_results <- list()
    

    if (method %in% c("secr", "both")) {
      secr.traps <- as.data.frame(traps)
      colnames(secr.traps) <- c("x", "y")
      secr_traps <- read.traps(data = secr.traps, detector = "proximity")
      
      capthist <- sim_data$capthist_df
      capt_hist <- make.capthist(capthist, secr_traps, fmt = "trapID")
      
      mask <- make.mask(secr_traps, buffer = buffer)
      fit_secr <- secr.fit(capt_hist, mask = mask)
      
      pred <- summary(fit_secr)$predicted
      df_secr <- data.frame(
        method    = "secr",
        parameter = rownames(pred),
        estimate  = pred[, "estimate"],
        se        = pred[, "SE.estimate"]
      )
      sim_results[["secr"]] <- df_secr
    }
    

    if (method %in% c("acre", "both")) {
      captu <- sim_data$capthist_df
      data_acre <- read.acre(captu, traps, control.mask = list(buffer = buffer))
      fit_acre <- fit.acre(data_acre)
      
      coef_acre <- summary(fit_acre)$coefs
      
      se_acre <- summary(fit_acre)$coefs_se |> as.vector()

      df_acre <- data.frame(
        method =  "acre",
        parameter = names(coef_acre), 
        estimate  = as.vector(coef_acre),
        se        = se_acre
      )
      
      sim_results[["acre"]] <- df_acre
    }
    
    if(detfn == "hn") {
      true_vals <- c(D = D, g0 = g0, sigma = sigma)
    } else if(detfn == "hhn") {
      true_vals <- c(D = D, lambda0 = lambda0, sigma = sigma)
    }    
    results[[i]] <- bind_rows(sim_results, .id = "fit_method") %>%
      mutate(sim = i,
         true = true_vals[parameter],
         pct_diff = 100 * ((estimate - true)/true))
  }
  
  all_results <- bind_rows(results)
  
  return(list(tidy = all_results, results = results))
}

sim_scr_data_multi <- function(D, g0, sigma, xlim, ylim, traps_list) {

  a_m2 <- diff(xlim) * diff(ylim)
  a_ha <- a_m2 / 10000
  
  N <- rpois(1, D * a_ha)
  
  activity_centers <- cbind(runif(N, xlim[1], xlim[2]),
                            runif(N, ylim[1], ylim[2]))
  
  capt_sessions <- vector("list", length(traps_list))
  
  for (s in seq_along(traps_list)) {
    traps <- traps_list[[s]]
    n_traps <- nrow(traps)
    
    capt_array <- array(0, dim = c(N, n_traps, n_occasions))
    
    dists <- crossdist(activity_centers[,1], activity_centers[,2],
                       traps[,1], traps[,2])
    
      p <- g0 * exp(-dists^2 / (2 * sigma^2))
      capt_array[,,i] <- matrix(rbinom(N * n_traps, 1, p), nrow = N)
 
    
    detected <- apply(capt_array, 1, sum) > 0
    capt_array <- capt_array[detected,,]
    
    capt_sessions[[s]] <- list(capt = capt_array,
                               traps = traps,
                               detected = sum(detected))
  }
  
  list(N_total = N,
       activity_centers = activity_centers,
       sessions = capt_sessions)
}

# Session - statistically independent capture data 
# Occasion - dependent - animals recognised between time & space - focus on single occasions for now.
# Acre - single occasion.
# Keep traps very separate to be independent. 
