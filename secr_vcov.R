
df1 <- sim_adapt_data_multi(beta0 = log(3), betas = c(0.5, 0.5), covs_session = covs_session, traps_list = traps_list, g0 = 0.8,
                     alpha0 = 4.05, alphas = c(-0.23, -0.23), detfn = 'hn', model = 'inh', xlim = c(-900, 1200), ylim = c(-900, 1200))

df1 

captu <- df1$capthist_df

secr.traps <- lapply(traps_list, function(df) {
  colnames(df) <- c("x", "y")
  read.traps(data = df, detector = "proximity")
})



capt_hist <- make.capthist(captu, secr.traps, 
                           fmt = "trapID")
sessioncov <- df1$session.cov

buffer = 300
mask <- make.mask(secr.traps, buffer = buffer)

scr <- secr.fit(capt_hist, model = list(sigma ~ x.s), CL = TRUE, sessioncov = sessioncov)

vcov(scr)

