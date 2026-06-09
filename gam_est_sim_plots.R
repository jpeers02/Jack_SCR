# ============================================================================
# DENSITY SURFACE PLOTS — 1D and 2D
# Works with sim_grf_500_gam_FULL2.rds (predictions on log-density scale,
# l_esas = 0, so exp(pred) = density)
# ============================================================================
library(dplyr)
library(ggplot2)
library(viridis)

sim_data <- readRDS("sim_grf_500_gam_FULL2.rds")

# ----------------------------------------------------------------------------
# 1D: estimated density curve vs true, with CI ribbon
# ----------------------------------------------------------------------------

grf_seq <- sim_data$grf_seq

# Mean estimated curve (back-transform log-density -> density)
mean_dens_gam <- colMeans(exp(sim_data$pred_1d_gam), na.rm = TRUE)
mean_dens_vp  <- colMeans(exp(sim_data$pred_1d_vp),  na.rm = TRUE)
true_dens     <- exp(sim_data$true_log_1d)

# Mean per-sim CI bounds (back-transformed), then averaged -> "typical CI"
ci_lo_gam <- colMeans(exp(sim_data$pred_1d_gam - 1.96 * sim_data$se_1d_gam), na.rm = TRUE)
ci_hi_gam <- colMeans(exp(sim_data$pred_1d_gam + 1.96 * sim_data$se_1d_gam), na.rm = TRUE)
ci_lo_vp  <- colMeans(exp(sim_data$pred_1d_vp  - 1.96 * sim_data$se_1d_vp),  na.rm = TRUE)
ci_hi_vp  <- colMeans(exp(sim_data$pred_1d_vp  + 1.96 * sim_data$se_1d_vp),  na.rm = TRUE)

df_1d <- data.frame(
  grf    = rep(grf_seq, 2),
  mean   = c(mean_dens_gam, mean_dens_vp),
  lo     = c(ci_lo_gam, ci_lo_vp),
  hi     = c(ci_hi_gam, ci_hi_vp),
  method = rep(c("GAM", "GAM + VarProp"), each = length(grf_seq))
)
df_true <- data.frame(grf = grf_seq, true = true_dens)

p_1d <- ggplot(df_1d, aes(x = grf)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = method), alpha = 0.2) +
  geom_line(aes(y = mean, colour = method), linewidth = 1) +
  geom_line(data = df_true, aes(y = true), colour = "black",
            linetype = "dashed", linewidth = 1) +
  scale_colour_manual(values = c("GAM" = "#9b59b6", "GAM + VarProp" = "#2ecc71")) +
  scale_fill_manual(values   = c("GAM" = "#9b59b6", "GAM + VarProp" = "#2ecc71")) +
  labs(x = "GRF covariate (scaled)", y = "Density (animals/ha)",
       title = "Estimated Density Along the GRF Gradient",
       subtitle = "Dashed = true; ribbon = mean 95% CI") +
  theme_classic(base_family = "Times New Roman")

ggsave("density_1d.png", p_1d, width = 8, height = 5, dpi = 150)

# ----------------------------------------------------------------------------
# 1D spaghetti: 500 individual estimated curves + mean + true
# ----------------------------------------------------------------------------

spag_df <- do.call(rbind, lapply(1:nrow(sim_data$pred_1d_vp), function(i){
  data.frame(grf = grf_seq, dens = exp(sim_data$pred_1d_vp[i,]), sim = i)
}))

p_spag <- ggplot() +
  geom_line(data = spag_df, aes(x = grf, y = dens, group = sim),
            colour = "grey70", alpha = 0.15, linewidth = 0.3) +
  geom_line(data = data.frame(grf = grf_seq, dens = mean_dens_vp),
            aes(x = grf, y = dens), colour = "#9b59b6", linewidth = 1.2) +
  geom_line(data = df_true, aes(x = grf, y = true),
            colour = "black", linetype = "dashed", linewidth = 1) +
  labs(x = "GRF covariate (scaled)", y = "Density (animals/ha)",
       title = "500 Estimated Density Curves (GAM + VarProp)",
       subtitle = "Purple = mean; dashed = true") +
  theme_classic(base_family = "Times New Roman")

ggsave("density_1d_spaghetti.png", p_spag, width = 8, height = 5, dpi = 150)

# ----------------------------------------------------------------------------
# 2D HEATMAP: True / GAM / GAM+VarProp side by side (density scale)
# ----------------------------------------------------------------------------

nd_grid <- sim_data$newdata_grid
sess_df <- data.frame(x = sim_data$sess.locs[,1], y = sim_data$sess.locs[,2])

mean_grid_gam <- colMeans(exp(sim_data$pred_grid_gam), na.rm = TRUE)
mean_grid_vp  <- colMeans(exp(sim_data$pred_grid_vp),  na.rm = TRUE)
true_grid     <- exp(sim_data$true_log_grid)

heat_df <- data.frame(
  x      = rep(nd_grid$x, 3),
  y      = rep(nd_grid$y, 3),
  dens   = c(true_grid, mean_grid_gam, mean_grid_vp),
  method = rep(c("True", "GAM", "GAM + VarProp"), each = nrow(nd_grid))
) %>%
  mutate(method = factor(method, levels = c("True", "GAM", "GAM + VarProp")))

p_heat <- ggplot(heat_df, aes(x = x, y = y, fill = dens)) +
  geom_tile() +
  geom_point(data = sess_df, aes(x = x, y = y), inherit.aes = FALSE,
             colour = "red", shape = 4, size = 1.5, stroke = 0.8) +
  facet_wrap(~ method, nrow = 1) +
  scale_fill_viridis_c(name = "Density\n(animals/ha)") +
  coord_equal() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(title = "Estimated Density Surface", x = "Easting", y = "Northing") +
  theme_minimal(base_family = "Times New Roman") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank())

ggsave("density_2d_heatmap.png", p_heat, width = 12, height = 4.5, dpi = 150)

# ----------------------------------------------------------------------------
# 2D UNCERTAINTY HEATMAP: mean SE across sims (VarProp), spatial
# ----------------------------------------------------------------------------

mean_se_grid_vp <- colMeans(sim_data$se_grid_vp, na.rm = TRUE)

unc_df <- data.frame(x = nd_grid$x, y = nd_grid$y, se = mean_se_grid_vp)

p_unc <- ggplot(unc_df, aes(x = x, y = y, fill = se)) +
  geom_tile() +
  geom_point(data = sess_df, aes(x = x, y = y), inherit.aes = FALSE,
             colour = "white", shape = 4, size = 1.5, stroke = 0.8) +
  scale_fill_viridis_c(name = "SE\n(log density)", option = "magma") +
  coord_equal() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(title = "Spatial Uncertainty (GAM + VarProp)", x = "Easting", y = "Northing") +
  theme_minimal(base_family = "Times New Roman") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank())

ggsave("density_2d_uncertainty.png", p_unc, width = 6, height = 5, dpi = 150)

cat("Saved: density_1d.png, density_1d_spaghetti.png,",
    "density_2d_heatmap.png, density_2d_uncertainty.png\n")