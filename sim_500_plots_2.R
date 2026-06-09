# ============================================================================
# COVERAGE: Spatial (2D heatmap)
# ============================================================================

# Compute 95% CI for each grid point
ci_lower_gam <- apply(sim_grf_500$pred_grid_gam, 2, quantile, 0.025)
ci_upper_gam <- apply(sim_grf_500$pred_grid_gam, 2, quantile, 0.975)

ci_lower_varprop <- apply(sim_grf_500$pred_grid_varprop, 2, quantile, 0.025)
ci_upper_varprop <- apply(sim_grf_500$pred_grid_varprop, 2, quantile, 0.975)

# Check if true value falls in CI
coverage_gam <- (true_D_centered >= ci_lower_gam) & (true_D_centered <= ci_upper_gam)
coverage_varprop <- (true_D_centered >= ci_lower_varprop) & (true_D_centered <= ci_upper_varprop)

# Prepare data
cov_plot <- data.frame(
  x = rep(newdata$x, 2),
  y = rep(newdata$y, 2),
  coverage = c(as.numeric(coverage_gam), as.numeric(coverage_varprop)),
  method = rep(c("GAM", "GAM + VarProp"), each = nrow(newdata))
)

sess_df <- data.frame(x = sim_grf_500$sess.locs[, 1], y = sim_grf_500$sess.locs[, 2])

# Plot
ggplot(cov_plot, aes(x = x, y = y, fill = coverage)) +
  geom_tile() +
  geom_point(data = sess_df, aes(x = x, y = y), 
             inherit.aes = FALSE, colour = "red", shape = 4, size = 2) +
  facet_wrap(~ method, nrow = 1) +
  scale_fill_gradient(name = "Coverage", low = "white", high = "darkgreen",
                      limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  coord_equal() +
  scale_x_continuous(expand = FALSE) +
  scale_y_continuous(expand = FALSE) +
  labs(title = "95% CI Coverage: True Value in Interval?") +
  theme_minimal(base_family = "Times New Roman") +
  theme(strip.text = element_text(size = 12, face = "bold"))

# ============================================================================
# COVERAGE: 1D curve with CI ribbon
# ============================================================================

ci_lower_1d_gam <- apply(pred_gam_1d, 2, quantile, 0.025)
ci_upper_1d_gam <- apply(pred_gam_1d, 2, quantile, 0.975)

ci_lower_1d_vp <- apply(pred_varprop_1d, 2, quantile, 0.025)
ci_upper_1d_vp <- apply(pred_varprop_1d, 2, quantile, 0.975)

# Check coverage at each GRF point
cov_1d_gam <- (true_1d >= ci_lower_1d_gam) & (true_1d <= ci_upper_1d_gam)
cov_1d_vp <- (true_1d >= ci_lower_1d_vp) & (true_1d <= ci_upper_1d_vp)

cat("GAM 1D coverage:", mean(cov_1d_gam), "\n")
cat("VarProp 1D coverage:", mean(cov_1d_vp), "\n")

# Plot with CI ribbon
par(mfrow = c(1, 2))

# GAM
plot(grf_seq, true_1d, type = "n", 
     main = "GAM: 95% CI Coverage",
     xlab = "GRF Covariate", ylab = "log Density",
     ylim = c(min(ci_lower_1d_gam), max(ci_upper_1d_gam)))
polygon(c(grf_seq, rev(grf_seq)), 
        c(ci_upper_1d_gam, rev(ci_lower_1d_gam)),
        col = rgb(0.8, 0.8, 0.8, 0.5), border = NA)
lines(grf_seq, colMeans(pred_gam_1d), col = "purple", lwd = 2)
lines(grf_seq, true_1d, col = "black", lwd = 2, lty = 2)
legend("topleft", c("True", "Mean", "95% CI"), 
       col = c("black", "purple", "gray"), lty = c(2, 1, NA), 
       fill = c(NA, NA, "gray"), lwd = 2)

# VarProp
plot(grf_seq, true_1d, type = "n", 
     main = "GAM + VarProp: 95% CI Coverage",
     xlab = "GRF Covariate", ylab = "log Density",
     ylim = c(min(ci_lower_1d_vp), max(ci_upper_1d_vp)))
polygon(c(grf_seq, rev(grf_seq)), 
        c(ci_upper_1d_vp, rev(ci_lower_1d_vp)),
        col = rgb(0.8, 0.8, 0.8, 0.5), border = NA)
lines(grf_seq, colMeans(pred_varprop_1d), col = "purple", lwd = 2)
lines(grf_seq, true_1d, col = "black", lwd = 2, lty = 2)
legend("topleft", c("True", "Mean", "95% CI"), 
       col = c("black", "purple", "gray"), lty = c(2, 1, NA), 
       fill = c(NA, NA, "gray"), lwd = 2)

par(mfrow = c(1, 1))