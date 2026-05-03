
library(tidyverse)
df <- read_csv("C:/Users/Jack/Downloads/varprop.csv")

df |> glimpse()
df %>%
  filter(grepl('^D', parameter)) %>%
  ggplot(aes(x = method, y = se, fill = method)) +
  geom_violin(alpha = 0.6) +
  geom_boxplot(width = 0.1, fill = 'white', alpha = 0.8) +
  facet_wrap(~ parameter, scales = 'free_y') +
  scale_fill_manual(
    name = "Model",
    values = c("secr_glm" = "#2ecc71", 
               "secr_joint" = "#e74c3c",
               "secr_varprop" = "#9b59b6"),
    labels = c("secr_glm" = "Two-stage GLM",
               "secr_joint" = "All-In-One",
               "secr_varprop" = "Two-stage LMM")
  ) +
  scale_x_discrete(labels = c("secr_glm" = "Two-stage GLM",
                              "secr_joint" = "All-In-One",
                              "secr_varprop" = "Two-stage LMM")) + 
  labs(x = 'Method', y =  'Standard Error', 
       title = 'Standard Errors for Coefficients by Model') +
  theme_bw()

df %>%
  filter(grepl('^D', parameter)) |> 
  ggplot(aes(x = time, fill = method)) + 
  geom_density(alpha = 0.5, color = NA) + 
  theme_classic() + 
  scale_fill_manual(name = 'Model', values = c('green', '#e74c3c', 'purple') ,labels = c('secr_glm' = 'Two Stage GLM', 
                               'secr_joint' = 'All-In-One Model', 
                               'secr_varprop' = 'Two Stage Linear Mixed Model' )) + 
  labs(x = 'Time', y = 'Density', 
       title = 'Time Distribution For Each Model')


coverage_df <- df %>%
  filter(grepl('^D', parameter)) %>%
  mutate(
    lower = estimate - 1.96 * se,
    upper = estimate + 1.96 * se,
    covered = true >= lower & true <= upper
  ) %>%
  group_by(method, parameter) %>%
  summarise(coverage = mean(covered), mean_se = mean(se), .groups = 'drop')

coverage_df %>%
  ggplot(aes(x = method, y = coverage, fill = method)) +
  geom_col() +
  geom_hline(yintercept = 0.95, linetype = 'dashed', colour = 'black') +
  facet_wrap(~ parameter) +
  scale_fill_manual(
    values = c("secr_glm" = "#2ecc71", 
               "secr_joint" = "#e74c3c",
               "secr_varprop" = "#9b59b6"),
    labels = c("secr_glm" = "Two-stage GLM",
               "secr_joint" = "All-In-One Model",
               "secr_varprop" = "Two-stage LMM")
  ) +
  scale_x_discrete(labels = c("secr_glm" = "Two-stage GLM",
                              "secr_joint" = "All-In-One",
                              "secr_varprop" = "Two-stage LMM")) +
  labs(y = "Coverage", x = "Model", title = "95% CI Coverage by Method") +
  theme_bw() +
  theme(legend.position = 'none')

# Test harmonic regression -> capturing seasonal changes
# Writing 
# GAMs - 
# Email Dave - Show all in one model - translated delta method to glmm for SCR - Show off plots - moving onto smoothing, like to replicate in mgcv - ask for some info on steps moving forward - give some code snippets

