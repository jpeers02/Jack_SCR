library(tidyverse)

trialdf <- read_csv("~/secrvsacre.csv")
true_vals <- c(D = 3, g0 = 0.75, sigma = 80)

trialdf <- trialdf %>%
  mutate(
    true = case_when(
      parameter %in% c("D") ~ true_vals["D"],
      parameter %in% c("g0") ~ true_vals["g0"],
      parameter %in% c("sigma") ~ true_vals["sigma"]
    ),
    pct_diff = 100 * ((estimate - true) / true)
  )
trialdf
trialdf %>%
  ggplot(aes(x = parameter, y = pct_diff, fill = parameter)) + 
  geom_violin(alpha = 0.4) + 
  facet_wrap(~method) + 
  theme_minimal() + 
  geom_hline(yintercept = 0, col = "purple") +
  geom_hline(yintercept = c(-20, 20), linetype = "dashed", color = "gray")

