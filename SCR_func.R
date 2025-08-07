
library(spatstat)

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


