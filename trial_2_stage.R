
source('SCR_func.R')

offsets <- expand.grid(dx = c(-700, 0, 700),
                       dy = c(-700, 0, 700))

traps_list <- lapply(1:nrow(offsets), function(i) {
  traps <- as.data.frame(traps.1)  
  traps$x <- traps[,1] + offsets$dx[i]
  traps$y <- traps[,2] + offsets$dy[i]
  traps <- traps[, c("x", "y")]     
  return(traps)
})

cov_x <- sapply(traps_list, function(tr) tr[5, "x"])
cov_y <- sapply(traps_list, function(tr) tr[5, "y"])

cov_x
cov_y


covs_session <- cbind(
  x = (cov_x - mean(cov_x)) / sd(cov_x),
  y = (cov_y - mean(cov_y)) / sd(cov_y)
)

colnames(covs_session) <- c('x.s', 'y.s')

covs_session

cov_x <- sapply(traps_list, function(tr) tr[5, "x"])
cov_x_std <- (cov_x - mean(cov_x)) / sd(cov_x)

covs_session <- matrix(cov_x_std, ncol = 1, dimnames = list(NULL, "x.s"))

covs_session

trial.1 <- model.sim_inh(100, beta0 = log(3), betas = c(0.5), 
              covs_session = covs_session, traps_list = traps_list, 
              alpha0 = 4.05, alphas = c(-0.23), 
              stage = 3, detfn = 'hn', buffer = 500, xlim = c(-900, 1200), ylim = c(-900, 1200), 
              method = 'secr', g0 = 0.8, D_formula = 'x.s', sigma_formula = 'x.s', spacing = 25
              )
write_csv(trial.1$tidy, 'trial_2.csv')

library(tidyverse)

df <- read_csv("~/2.stage.model_trial.csv")
df <- df %>%
  mutate(estimate = ifelse(parameter == "g0", plogis(estimate), estimate)) %>%
  filter(estimate != 0)

df %>%
  filter(parameter == 'D.x.s') %>%
  ggplot(aes(x = method, y = estimate, fill = method)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.5, linetype = 'dashed', colour = 'red') +
  labs(title = "D.x.s estimates by method", y = "Estimate", x = "Method") + 
  theme_bw()

df %>%
  filter(parameter == 'D.y.s') %>%
  ggplot(aes(x = method, y = estimate, fill = method)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.5, linetype = 'dashed', colour = 'red') +
  labs(title = "D.y.s estimates by method", y = "Estimate", x = "Method") + 
  theme_bw()

df %>%
  ggplot(aes(x = method, y = estimate, fill = method)) +
  geom_boxplot() +
  geom_hline(data = df %>% 
               filter(estimate != 0) %>%
               distinct(parameter, true), 
             aes(yintercept = true), linetype = 'dashed', colour = 'red') +
  facet_wrap(~ parameter, scales = 'free_y') +
  labs(title = "Estimates by method and parameter", y = "Estimate", x = "Method") +
  theme_bw() +
  theme(legend.position = 'none')


## Parallel Processing 