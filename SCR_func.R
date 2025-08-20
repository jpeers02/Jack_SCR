
library(spatstat)
library(secr)
library(acre)
library(tidyverse)
load("C:/Users/Jack/Downloads/test-data.RData")

# add data set to github

traps.1 <- test.data$traps

## Simulating a single session & occasion data set

sim_scr_data <- function(D, g0, sigma, xlim, ylim, traps) {
  a_m2 <- diff(xlim) * diff(ylim)
  a_ha <- a_m2 / 10000
  N <- rpois(1, D * a_ha)
  
  activity_centers <- cbind(runif(N, xlim[1], xlim[2]),
                            runif(N, ylim[1], ylim[2]))
  
  dists <- crossdist(activity_centers[, 1], activity_centers[, 2],
                     traps[, 1], traps[, 2])
  
  p <- g0 * exp(-dists^2 / (2 * sigma^2))
  capt <- matrix(rbinom(length(p), 1, p), nrow = N)
  
  capt_filt <- capt[rowSums(capt) > 0, ]
  list(capt = capt_filt, N = N, activity_centers = activity_centers)
}

## Adapting the data set into  format required for secr and acre

sim_adapt_data <- function(D, g0, sigma, xlim, ylim, traps, session_id = 1) {
  
  sim <- sim_scr_data(D, g0, sigma, xlim, ylim, traps)
  
  capt <- sim$capt
  capthist_df <- data.frame()
  for (i in 1:nrow(capt)) {
    for (j in 1:ncol(capt)) {
      if (capt[i, j] > 0) {
        capthist_df <- rbind(capthist_df, 
                             data.frame(session = session_id,
                                        ID = paste0("ind", i),
                                        trap = j,
                                        stringsAsFactors = FALSE, occasion = 1))
      }
    }
  }
  
  return(list(sim_data = sim, capthist_df = capthist_df))
}

## Simulating data sets and fitting models to those data sets

acre_model.sim <- function(n_sims, D, g0, sigma, 
                           xlim, ylim, traps, buffer) {
  results <- vector("list", n_sims)
  
  for (i in 1:n_sims) {
    sim_data <- sim_adapt_data(D = D, g0 = g0, sigma = sigma, 
                               xlim = xlim, ylim = ylim, traps = traps, session_id = 1)
    captu <- sim_data$capthist_df
    df <- read.acre(captu, traps, control.mask = list(buffer = buffer))
    fit <- fit.acre(df)
    
    results[[i]] <- summary(fit)$coefs
  }
  
  all_results <- bind_rows(results, .id = "sim") %>% mutate(D = as.numeric(D), 
                                                            g0 = as.numeric(g0), 
                                                            sigma = as.numeric(sigma))
  all_results_tidy <- all_results %>%
    pivot_longer(cols = c("g0", "D", "sigma"), 
                 names_to = "parameter", 
                 values_to = "estimate")
  
  return(list(tidy = all_results_tidy, wide = all_results, results = results))

}

# Secr Simulation

# Start by simulating multi-session data - list of traps 

# Buffer = 5 * Sigma

# add hazard-half normal detection function simulation - good for acoustics.

# plots of % difference from truth

