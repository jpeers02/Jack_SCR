cov_x <- sapply(traps_list, function(tr) tr[5, "x"])
cov_y <- sapply(traps_list, function(tr) tr[5, "y"])

cov_x
cov_y


covs_session <- cbind(
  x = (cov_x - mean(cov_x)) / sd(cov_x),
  y = (cov_y - mean(cov_y)) / sd(cov_y)
)




sim_scr_data_multi_covsession(beta0 = log(3), betas = c(log(2), log(3)), covs_session = covs_session, g0 = 0.8, sigma = 75, 
                              xlim = c(-900, 1200), ylim = c(-900, 1200), traps_list)


trial1 <- sim_adapt_data_multi(g0 = 0.8, sigma = 75, xlim = c(-900, 1200), ylim = c(-900, 1200),beta0 = log(3), 
                     betas = c(log(2), log(3)), covs_session = covs_session, traps_list = traps_list, model = 'inh')


captu <- trial1$capthist_df
captu
sessioncov <- trial1$session.cov

sessioncov
colnames(sessioncov) <- c('x.s', 'y.s', 'session')

sessioncov

data_acre <- read.acre(captu, traps_list, control.mask = list(buffer = 400), session.cov = sessioncov)
# fit_acre <- fit.acre(data_acre)
# summary(fit_acre)
fi2 <- fit.acre(data_acre, model=list(D=~x.s + y.s)) # 

esas = fi2$coefficients[substr(names(fi2$coefficients), 1, 3) == 'esa']
n = sapply(data_acre$capt, function(x) nrow(x$bincapt))
x = sessioncov$x.s
y = sessioncov$y.s
trialglm = glm(n ~ x + y, offset = log(esas), family = poisson)
summary(trialglm)

# Next: Modelling & simulating sigma as a function of x and y
# Control alpha to be 30 to 100 (distance between traps)
# 9 different models to detection data. 


fi2 <- fit.acre(data_acre, fix=list(D = 1)) # 

trial2 <- sim_adapt_data_multi(g0 = 0.8, xlim = c(-900, 1200), ylim = c(-900, 1200),beta0 = log(5), 
                               betas = c(0.5, 0.5), covs_session = covs_session, traps_list = traps_list, 
                               model = 'inh', alpha0 = 4.05, alphas =c(-0.23, -0.23))



captu1 <- trial2$capthist_df

sessioncov1 <- trial2$session.cov

colnames(sessioncov1) <- c('x.s', 'y.s', 'session')
colnames(covs_session) <- c('x.s', 'y.s')

data.acre2 <- read.acre(captu1, traps_list, control.mask = list(buffer = 100 * 5), session.cov = sessioncov1)

fit3 <-  fit.acre(data.acre2, model=list(sigma =~ x.s + y.s))


esas2 = fit3$coefficients[substr(names(fit3$coefficients), 1, 3) == 'esa']
n2 = sapply(data.acre2$capt, function(x) nrow(x$bincapt))
x <- sessioncov1$x.s
y <- sessioncov1$y.s
trialglm2 = glm(n2 ~ x + y, offset = log(esas2), family = poisson)
summary(trialglm2)
summary(fit3)

  
fit4 <-  fit.acre(data.acre2, model=list(D=~ x.s + y.s, sigma =~ x.s + y.s))

summary(fit4)

# Secr - make same model in secr - using a conditional likelihood to single out sigma. 
  
  

  
  
  
