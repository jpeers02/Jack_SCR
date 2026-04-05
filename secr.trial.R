
source("~/SCR_func.R")

secr.data <- sim_adapt_data(D = 10, g0 = 0.9, sigma = 80, xlim = c(-500, 900), 
                           ylim = c(-500, 900), traps = traps.1, detfn = 'hn')

secr.traps <- as.data.frame(traps.1)

colnames(secr.traps) <- c("x", "y")

secr.traps

secr_traps <- read.traps(data = secr.traps, detector = "proximity")

capthist <- secr.data$capthist_df

capt_hist <- make.capthist(capthist, secr_traps, fmt = "trapID")

mask <- make.mask(secr_traps, buffer = 400)

fit <- secr.fit(capt_hist, 
                mask = mask)

vals <- summary(fit)$predicted
params <- rownames(summary(fit)$predicted)

df <- data.frame(
  parameter = vals[, 'estimate'],
  estimate = rownames(vals), 
  se = vals[, 'SE.estimate']
  
)
params['estimate']
plot(capt_hist, varycol = FALSE)

