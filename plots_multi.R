
multidata <- read_csv("~/multidatatest.csv")

multidata %>% 
  glimpse()
multidata %>% ggplot(aes(x = parameter, y = pct_diff, fill = parameter)) +
  geom_violin(alpha = 0.4) +
  ylab("% Difference") + 
  facet_wrap(~method) +
  xlab("Parameter") +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  geom_hline(yintercept = 0, col = "purple") +
  geom_hline(yintercept = c(-20, 20), linetype = "dashed", color = "gray")



