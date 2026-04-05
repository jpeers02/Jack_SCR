library(tidyverse)
library(acre)
library(secr)
sim_results <- sim_adapt_data(D=0.75, g0=0.9, sigma=75, xlim=c(-500,900), ylim=c(-500,900), traps=test.data$traps, session_id=1)
captu <- sim_results$capthist_df

acredf <- read.acre(captu, test.data$traps, control.mask = list(buffer = 300))
acredf
n_sims <- 250
results <- vector("list", n_sims)


  
  sim_data <- sim_adapt_data(D=5, g0=0.75, sigma=75, xlim=c(-500,900), ylim=c(-500,900), traps=test.data$traps, session_id=1, detfn = 'hn') 
  captu <- sim_data$capthist_df
  df <- read.acre(captu, test.data$traps, control.mask = list(buffer = 300))
  fit <- fit.acre(df)

coef_acre <- summary(fit)$coefs

se_acre <- summary(fit)$coefs_se |> as.vector()
coef_acre
df_acre <- data.frame(
  params = names(coef_acre), 
  estimate  = as.vector(coef_acre),
  se        = se_acre
)
names(coef_acre)

  
results
all_results <- bind_rows(results, .id = "sim") %>% mutate(D = as.numeric(D), g0 = as.numeric(g0), sigma = as.numeric(sigma))
summary(all_results$sigma)
all_results_1 <- all_results %>%
  pivot_longer(cols = c("g0", "D", "sigma"), 
               names_to = "parameter", 
               values_to = "estimate")

ggplot(all_results_1, aes(x = parameter, y = estimate)) +
  geom_jitter(width = 0.2, alpha = 0.6, color = "steelblue") +
  facet_wrap(~ parameter, scales = "free_y") +
  labs(title = "Jitter plots of parameter estimates",
       x = "",
       y = "Estimate") +
  theme_minimal() +
  theme(strip.text = element_text(size = 12))

secr_data <- sim_adapt_data(5, 0.9, 75, c(-500, 900), c(-500, 900), test.data$traps, session_id = 1)
secr_capt <- secr_data$capthist_df
trapdf <- as.data.frame(test.data$traps)
colnames(trapdf) <- c("x", "y")
traps_secr <- read.traps(data = trapdf, detector = "proximity")
capt_hist <- make.capthist(secr_capt, traps = traps_secr, fmt = "trapID")
capt_hist 
secr_fit <- secr.fit(capt_hist, buffer = 300, start = list(D = 5, g0 = 0.9, sigma = 75))
