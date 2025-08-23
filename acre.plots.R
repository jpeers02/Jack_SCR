library(tidyverse)

# acre.sim <- acre_model.sim(500, 3, 0.75, 60, c(-500, 900), c(-500, 900), traps.1, buffer = 300)


# Will read the data from csv file if using again rather than running the simulation again.

# wide <- as.data.frame(acre.sim$wide)
# acre.sim$results
# write_csv(wide, "acre.wide.csv")

# tidy <- as.data.frame(acre.sim$tidy)

# write_csv(tidy, "acre.tidy.csv")

# tidy
pctmulti <- read_csv("~/pctdiff_multi.csv")
tidymulti <- read_csv("~/multitidy_scr.csv")
tidy <- read_csv("~/acre.tidy.csv")
wide <- read_csv("~/acre.wide.csv")
# Computing % difference from true values

wide <- wide %>% 
  mutate(p.g0 = ((g0 - 0.75)/0.75) * 100, 
         p.D = ((D - 3)/3) * 100, 
         p.sig = ((sigma - 60)/60) * 100)

wide %>% glimpse()
p.tidy <- wide %>% 
  pivot_longer(
    cols = c("p.g0", "p.D", "p.sig"),   
    names_to = "parameter",             
    values_to = "estimate"              
  )

ggplot(p.tidy, aes(x = parameter, y = estimate, col = parameter)) +
  geom_boxplot(alpha = 0.3) +
  geom_jitter(width = 0.2, size = 2) +
  ylab("% Difference") +
  xlab("Parameter") +
  scale_x_discrete(labels = c("p.g0" = "g0", "p.D" = "D", "p.sig" = "Sigma")) +
  scale_color_manual(
    values = c("p.g0" = "green", "p.D" = "red", "p.sig" = "blue"),
    labels = c("p.g0" = "g0", "p.D" = "D", "p.sig" = "Sigma")
  ) +
  theme_minimal() +
  theme(legend.title = element_blank()) +
    geom_hline(yintercept = 0, col = "purple", size = 1) 

ggplot(p.tidy, aes(x = estimate)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  facet_wrap(~ parameter, scales = "free") +
  xlab("% Difference") +
  theme_minimal()

ggplot(p.tidy, aes(x = parameter, y = estimate, fill = parameter)) +
  geom_violin(alpha = 0.4) +
  scale_fill_manual(values = c("green", "red", "blue")) +
  ylab("% Difference") +
  xlab("Parameter") +
  ggtitle("Parameters % Difference from True Value") +
  scale_x_discrete(labels = c("p.g0" = "g0", "p.D" = "D", "p.sig" = "Sigma")) +
  scale_fill_manual(
    values = c("p.g0" = "green", "p.D" = "red", "p.sig" = "blue"),
    labels = c("p.g0" = "g0", "p.D" = "D", "p.sig" = "Sigma")
  ) +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  geom_hline(yintercept = 0, col = "purple", size = 1) +
  geom_hline(yintercept = c(-20, 20), linetype = "dashed", color = "gray")

pctmulti %>% ggplot(aes(x = parameter, y = estimate, fill = parameter)) +
  geom_violin(alpha = 0.4) +
  ylab("% Difference") + 
  facet_wrap(~sim_id, 
             labeller = labeller(sim_id = function(x) paste0("D = ", x))) +
  xlab("Parameter") +
  ggtitle("Parameters % Difference from True Value - Different D values") +
  scale_x_discrete(labels = c("p.g0" = "g0", "p.D" = "D", "p.sig" = "Sigma")) +
  scale_fill_manual(
    values = c("p.g0" = "green", "p.D" = "red", "p.sig" = "blue"),
    labels = c("p.g0" = "g0", "p.D" = "D", "p.sig" = "Sigma")
  ) +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  geom_hline(yintercept = 0, col = "purple") +
  geom_hline(yintercept = c(-20, 20), linetype = "dashed", color = "gray")

