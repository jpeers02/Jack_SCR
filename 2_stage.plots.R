

df <- trial.1$tidy

df %>%
  ggplot(aes(x = method, y = estimate, fill = method)) +
  geom_violin() +
  geom_hline(data = df %>% 
               filter(estimate != 0) %>%
               distinct(parameter, true), 
             aes(yintercept = true), linetype = 'dashed', colour = 'red') +
  facet_wrap(~ parameter, scales = 'free_y') +
  labs(title = "Estimates by method and parameter", y = "Estimate", x = "Method") +
  theme_bw() +
  theme(legend.position = 'none')


df %>% 
  distinct(sim, method, time) %>%
  ggplot(aes(x = time, fill = method)) + 
  geom_density(alpha = 0.6, position = 'identity') +
  scale_fill_manual(values = c("secr_glm" = "purple", "secr_joint" = "green"),
                    labels = c("secr_glm" = "Two-stage GLM", "secr_joint" = "Joint Model")
  ) +
  labs(title = "Computation time by method", x = "Time (s)", y = "Count") +
  theme_minimal()


df %>%
  filter(unique(time))
