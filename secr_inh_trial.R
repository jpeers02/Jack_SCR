

library(secr)


secr.traps <-  lapply(traps_list, function(df) {
  colnames(df) <- c("x", "y")
  read.traps(data = df, detector = "proximity")
})


sim_data <- trial2
capthist <- sim_data$capthist_df
capt_hist <- make.capthist(capthist, secr.traps, fmt = "trapID")

mask <- make.mask(secr.traps, buffer = 500)
fit_secr <- secr.fit(capt_hist, mask = mask, model = list(sigma ~ x.s + y.s),
                     sessioncov = sessioncov1, CL = TRUE)

# Simulate capture proportion of time CI captures the true vlaue - do 100 - 1000 simulates of all in one vs two stage. 

summary(fit_secr)$coef
esa(fit_secr)
esa_s <- (numeric(9))

for(i in 1:9){
  esa_s[i] <- esa(fit_secr, sessnum = i)[1] 
  
}
esa_s

n3 <- sapply(fit_secr$capthist, function(x) dim(x)[1])

n2 = sapply(fit_secr$capt, function(x) nrow(x$bincapt))
capt_hist

covs_session
x
trialglm3 <- glm(n3 ~ x + y, offset = log(esa_s), family = poisson)
summary(trialglm3)
summary(fit4)
