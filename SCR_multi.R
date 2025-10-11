source('~/multisession.test.R')
traps <- traps_list

sim_data <- sim_adapt_data_multi(3, sigma = 70, traps_list =  traps, g0 = 0.8, xlim = c(-900, 1200), ylim = c(-900, 1200))
buffer = 350
captu <- sim_data$capthist_df
data_acre <- read.acre(captu, traps, control.mask = list(buffer = buffer))
fit_acre <- fit.acre(data_acre)
data_acre$traps

secr_traps_list <- lapply(traps, function(df) {
  colnames(df) <- c("x", "y")
  read.traps(data = df, detector = "proximity")
})

capt_hist <- make.capthist(captu, secr_traps_list, fmt = "trapID")

mask <- make.mask(secr_traps_list, buffer = 350)
fit_secr <- secr.fit(capt_hist, mask = mask)

class(fit_secr)
fit_secr


trial.1 <- model.sim(n_sims = 50, D = 5, g0 = 0.8, sigma = 80, traps = traps_list, buffer = 400, detfn = 'hn', method = 'both', xlim = c(-900, 1200), ylim = c(-900, 1200))
