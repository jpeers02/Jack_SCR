set.seed(839151506)

sim_scr_data_multi_covsession(3, beta1 = 2, cov_session = rbinom(9, 1, 0.5),
                              g0 = 0.8, sigma = 80, xlim = c(-900, 1200), ylim = c(-900, 1200), traps_list = traps_list)

cov_ses <- rbinom(9, 1, 0.5)

df <- sim_scr_data_multi_covsession(3, beta1 = 2, cov_session = cov_ses,
                              g0 = 0.8, sigma = 80, xlim = c(-900, 1200), ylim = c(-900, 1200), traps_list = traps_list)

df1 <- sim_adapt_data_multi(beta0 = log(3), beta1 = 1.5, cov_session = cov_ses,
                                    g0 = 0.8, sigma = 80, xlim = c(-900, 1200), ylim = c(-900, 1200), detfn = 'hn', traps_list = traps_list, model = 'inh')

captu <- df1$capthist_df
sessioncov <- df1$session.cov
sessioncov

data_acre <- read.acre(captu, traps_list, control.mask = list(buffer = 400), session.cov = sessioncov)
fit_acre <- fit.acre(data_acre, model = list(D =~ cov))

fit_acre$coefficients
trial1 <- model.sim_inh(5, D0 = 1, beta1 = 3, cov_session = cov_ses, g0 = 0.8, sigma = 80,
              xlim = c(-900, 1200), ylim = c(-900, 1200), detfn = 'hn', traps = traps_list, buffer = 400, method = 'acre')
trial1

fit2 <- fit.acre(data_acre)
esas = fit2$coefficients[substr(names(fit2$coefficients), 1, 3) == 'esa']
df1$capthist_df
data_acre$capt
n = sapply(data_acre$capt, function(x) nrow(x$bincapt))
covs = sessioncov$cov
trialglm = glm(n ~ covs, offset = log(esas), family = poisson)
summary(trialglm)
covs

# Changing over space 
# Simulate a gaussian field 


