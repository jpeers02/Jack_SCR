
covs_session_x <- covs_session[, -2, drop = FALSE]
colnames(covs_session_x) <- 'x.s'

trial_varprop <- model.sim_inh(100, beta0 = log(3), betas = c(0.5), 
                               covs_session = covs_session_x, traps_list = traps_list,
                               alpha0 = 4.05, alphas = c(-0.23),
                               stage = 3, detfn = 'hn', buffer = 500,
                               xlim = c(-900, 1200), ylim = c(-900, 1200),
                               method = 'secr', g0 = 0.8,
                               D_formula = 'x.s', sigma_formula = 'x.s',
                               varprop = TRUE)

trial_varprop$tidy


library(patchwork)

# ---- Data Setup ----
set.seed(123)
forest <- rbinom(9, 1, 0.3)

covs_session <- cbind(covs_session, forest)


# simulation data for visualisation plots
df1 <- sim_adapt_data_multi(beta0 = log(0.5), betas = c(0.5, 0.5, 0.5), 
                            covs_session = covs_session_x, traps_list = traps_list, 
                            g0 = 0.8, alpha0 = 4.05, alphas = c(-0.23, -0.23, -0.23), 
                            detfn = 'hn', model = 'inh', 
                            xlim = c(-900, 1200), ylim = c(-900, 1200))

# extract activity centers and detections
all_centers <- do.call(rbind, lapply(df1$sim_data$sessions, function(s) s$activity_centers))
detected_centers <- do.call(rbind, lapply(df1$sim_data$sessions, function(s){
  s$activity_centers[rowSums(s$capt) > 0, ]
}))
all_traps <- do.call(rbind, traps_list)

# data frames
df_all <- as.data.frame(all_centers)
colnames(df_all) <- c("x", "y")

df_detected <- as.data.frame(detected_centers)
colnames(df_detected) <- c("x", "y")

df_traps <- as.data.frame(all_traps)
colnames(df_traps) <- c("x", "y")

# density surface grid
grid <- expand.grid(
  x = seq(-900, 1200, length.out = 200),
  y = seq(-900, 1200, length.out = 200)
)
grid$D <- exp(log(0.5) + 0.5 * (grid$x - mean(cov_x)) / sd(cov_x))

# ---- Colour scheme ----
cols <- c("secr_glm"     = "#2ecc71",
          "secr_joint"   = "#e74c3c",
          "secr_varprop" = "#9b59b6")

labs <- c("secr_glm"     = "Two-stage GLM",
          "secr_joint"   = "All-In-One",
          "secr_varprop" = "Two-stage LMM")

# ---- Plot 1: 3 panel density/detection ----
p1 <- ggplot() +
  geom_tile(data = grid, aes(x = x, y = y, fill = D)) +
  geom_point(data = df_traps, aes(x = x, y = y),
             colour = "red", shape = 4, size = 2, stroke = 1.2) +
  scale_fill_viridis_c(name = "Density\n(animals/ha)") +
  labs(title = "True Density Surface", x = "X coordinate", y = "Y coordinate") +
  theme_classic()

p2 <- ggplot() +
  geom_point(data = df_all, aes(x = x, y = y), 
             colour = "steelblue", alpha = 0.3, size = 0.8) +
  geom_point(data = df_traps, aes(x = x, y = y),
             colour = "red", shape = 4, size = 2, stroke = 1.2) +
  labs(title = paste0("True Population (N = ", nrow(df_all), ")"), 
       x = "X coordinate", y = "Y coordinate") +
  theme_classic()

p3 <- ggplot() +
  geom_point(data = df_detected, aes(x = x, y = y), 
             colour = "steelblue", alpha = 0.5, size = 0.8) +
  geom_point(data = df_traps, aes(x = x, y = y),
             colour = "red", shape = 4, size = 2, stroke = 1.2) +
  labs(title = paste0("Detected Animals (n = ", nrow(df_detected), ")"), 
       x = "X coordinate", y = "Y coordinate") +
  theme_classic()

plot_density_detection <- p1 + p2 + p3 + 
  plot_layout(guides = 'collect') +
  plot_annotation(title = "Spatially Varying Density and Detection")

