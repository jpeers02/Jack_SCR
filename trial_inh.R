
ndf <- sim_adapt_data_multi(D = seq(5, 20, length = 9), 
                   g0 = 0.8, sigma =80, xlim = c(-900, 1200), 
                   ylim = c(-900, 1200), traps_list = traps_list, detfn = 'hn')

secr.traps <-  lapply(traps_list, function(df) {
  colnames(df) <- c("x", "y")
  read.traps(data = df, detector = "proximity")
})
secr.traps
capthist <- ndf$capthist_df
capt_hist <- make.capthist(capthist, secr.traps, fmt = "trapID")

mask <- make.mask(secr.traps, buffer = 400)
fit_secr <- secr.fit(
  capt_hist,
  mask = mask,
  model = list(D ~ session)
)

rd <- rbinom(9, 1, 0.5)
est.d <- sapply(summary(fit_secr)$predicted, function(x) x$estimate[1]) 
r.d <- seq(5, 20, length = 9)
plot(r.d, est.d, col = 'blue')
abline(0, 1, col = 'red')
esa(fit_secr)
captu <- ndf$capthist_df

data_acre <- read.acre(captu, traps_list, control.mask = list(buffer = 400))
fit_acre <- fit.acre(data_acre)
