# randomly sample true parameter values

library(tidyverse)
n_scenarios <- 100

beta0_vals <- runif(n_scenarios, log(1), log(10))  # D range 1-10
betas_vals <- runif(n_scenarios, 0.1, 1.0)         # gradient strength

# run one sim per scenario
scenario_results <- mapply(function(b0, b1){
  model.sim_inh(1, beta0 = b0, betas = c(b1), 
                covs_session = covs_session, traps_list = traps_list, 
                alpha0 = 4.05, alphas = c(-0.23), 
                stage = 3, detfn = 'hn', buffer = 500, xlim = c(-900, 1200), ylim = c(-900, 1200), 
                method = 'secr', g0 = 0.8, D_formula = 'x.s', sigma_formula = 'x.s', spacing = 25
  )$tidy
}, beta0_vals, betas_vals, SIMPLIFY = FALSE)

scenario_df <- bind_rows(scenario_results) %>%
  mutate(beta0_true = rep(beta0_vals, each = nrow(scenario_results[[1]])),
         beta1_true = rep(betas_vals, each = nrow(scenario_results[[1]])))


scenario_df <- read_csv("~/monte_2_stage_sim.csv")

# coverage by method
scenario_df %>%
  filter(grepl('^D', parameter)) %>%
  mutate(
    lower = estimate - 1.96 * se,
    upper = estimate + 1.96 * se,
    covered = true >= lower & true <= upper
  ) %>%
  group_by(method, parameter) %>%
  summarise(coverage = mean(covered), .groups = 'drop')


scenario_df %>%
  filter(grepl('^D', parameter)) %>%
  ggplot(aes(x = true, y = estimate, colour = method)) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', colour = 'red') +
  facet_wrap(~ parameter, scales = 'free') +
  theme_classic() + 
  labs(x = "True", y = 'Estimate') + 
  geom_point(alpha = 0.5, position = position_jitter(width = 0.01)) +
  scale_colour_manual(name = 'Model',
                      values = c("secr_glm" = "purple", "secr_joint" = "green"))


scenario_df %>%
  distinct(method, time, beta0_true, beta1_true) %>%
  mutate(mean_D = exp(beta0_true)) %>%
  ggplot(aes(x = mean_D, y = time, colour = method)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = 'lm') +
  scale_colour_manual(values = c("secr_glm" = "purple", "secr_joint" = "green"),
                      labels = c("secr_glm" = "Two-stage Model", "secr_joint" = "All-In-One Model")) +
  labs(x = "Baseline Density (animals/ha)", y = "Time (s)", 
       title = "Computation Time vs Baseline Density",
       colour = "Model") +
  theme_classic()

# Next steps 



scenario_df %>%
  distinct(method, time, beta1_true) %>%
  ggplot(aes(x = beta1_true, y = time, colour = method)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = 'lm') +
  scale_colour_manual(values = c("secr_glm" = "purple", "secr_joint" = "green"),
                      labels = c("secr_glm" = "Two-stage Model", "secr_joint" = "All-In-One Model")) +
  labs(x = "Beta1 (gradient strength)", y = "Time (s)", colour = "Model") +
  theme_classic()

scenario_df %>%
  filter(grepl('^D', parameter)) %>%
  ggplot(aes(x = method, y = se, fill = method)) +
  geom_violin(alpha = 0.6) +
  geom_boxplot(width = 0.1, fill = 'white', alpha = 0.8) +
  facet_wrap(~ parameter, scales = 'free_y') +
  scale_fill_manual(values = c("secr_glm" = "purple", "secr_joint" = "green"),
                    labels = c("secr_glm" = "Two-stage Model", "secr_joint" = "All-In-One Model")) +
  scale_x_discrete(labels = c("secr_glm" = "Two-stage", "secr_joint" = "All-In-One")) +
  labs(y = "Standard Error", x = "Model", fill = "Model") +
  theme_classic()



