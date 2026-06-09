library(ggplot2)
library(geoR)
library(dplyr)

# ============================================================================
# Load simulation results
# ============================================================================

sim_grf_500 <- readRDS("sim_grf_500_gam_FULL.rds")

# ============================================================================
# Define D.calc function
# ============================================================================

D.calc <- function(grf_val){
  b0 <- -log(44426388/200000)
  alpha <- 1
  gamma <- 3
  tau <- 0
  exp(b0 + alpha * (1 - exp(-gamma * (grf_val - tau)))) + 0.005
}

# ============================================================================
# Regenerate GRF and spatial grid
# ============================================================================

set.seed(506)
x.range <- c(0, 100000)
y.range <- c(0, 100000)

# Create grid
xs <- seq(x.range[1], x.range[2], length.out = 100)
ys <- seq(y.range[1], y.range[2], length.out = 100)
xx <- rep(xs, length(ys))
yy <- rep(ys, each = length(xs))
mask.full <- cbind(xx, yy)

# Regenerate GRF
grf_sim <- grf(nrow(mask.full),
               grid = mask.full,
               cov.model = "matern",
               cov.pars = c(2, 100000 * 1.2),
               kappa = 1.5)
grf_vals <- grf_sim$data

# ============================================================================
# Compute nearest grid indices and true density
# ============================================================================

newdata <- expand.grid(
  x = xs,
  y = ys
)

nearest_grid <- apply(newdata, 1, function(pt){
  dists <- sqrt((mask.full[,1] - pt[1])^2 + (mask.full[,2] - pt[2])^2)
  which.min(dists)
})

nearest_sess <- apply(sim_grf_500$sess.locs, 1, function(pt){
  dists <- sqrt((mask.full[,1] - pt[1])^2 + (mask.full[,2] - pt[2])^2)
  which.min(dists)
})
grf_center <- mean(grf_vals[nearest_sess])
grf_scale <- sd(grf_vals[nearest_sess])

# Use SCALED GRF for true_D
grf_scaled <- (grf_vals[nearest_grid] - grf_center) / grf_scale
log_D_true <- log(D.calc(grf_scaled))
true_D_centered <- log_D_true - mean(log_D_true)

# ============================================================================
# Extract and aggregate GAM predictions
# ============================================================================

mean_grid_gam <- colMeans(sim_grf_500$pred_grid_gam)
mean_grid_varprop <- colMeans(sim_grf_500$pred_grid_varprop)

# ============================================================================
# Prepare plotting data
# ============================================================================

plot_df <- data.frame(
  x = rep(newdata$x, 3),
  y = rep(newdata$y, 3),
  fitted = c(mean_grid_gam, mean_grid_varprop, true_D_centered),
  method = rep(c("GAM", "GAM + VarProp", "True"), each = nrow(newdata))
) %>%
  mutate(method = factor(method, levels = c("True", "GAM", "GAM + VarProp")))

# Session centroids
sess.locs <- sim_grf_500$sess.locs
sess_df <- data.frame(x = sess.locs[, 1], y = sess.locs[, 2])

# ============================================================================
# Plot
# ============================================================================

dens_plot <- ggplot(plot_df, aes(x = x, y = y, fill = fitted)) +
  geom_tile() +
  geom_point(data = sess_df, aes(x = x, y = y), 
             inherit.aes = FALSE, colour = "red", shape = 4, size = 2) +
  facet_wrap(~ method, nrow = 1) +
  scale_fill_viridis_c(name = "log Density (centered)") +
  coord_equal() +
  scale_x_continuous(expand = FALSE) +
  scale_y_continuous(expand = FALSE) +
  labs(title = "Estimated Density Surface: GAM Results (500 sims)") +
  theme_minimal(base_family = "Times New Roman") +
  theme(strip.text = element_text(size = 12, face = "bold"),
        legend.position = "right")

# Extract data
pred_gam_1d <- sim_grf_500$pred_gam_1d_density
pred_varprop_1d <- sim_grf_500$pred_varprop_1d_density
grf_seq <- sim_grf_500$grf_seq

# Compute true curve (centered on log scale)
grf_scaled <- (grf_seq - grf_center) / grf_scale
log_D_true <- log(D.calc(grf_scaled))
true_1d <- log_D_true - mean(log_D_true)

# ============================================================================
# PLOT 1: GAM
# ============================================================================

plot(grf_seq, true_1d, type = "n", 
     main = "GAM: 500 Simulations",
     xlab = "GRF Covariate (scaled)", ylab = "log Density (centered)",
     ylim = c(min(pred_gam_1d), max(pred_gam_1d)))

# Spaghetti lines
for(i in 1:500) {
  lines(grf_seq, pred_gam_1d[i, ], col = rgb(0.5, 0.5, 0.5, 0.1))
}

# Mean
lines(grf_seq, colMeans(pred_gam_1d), col = "purple", lwd = 2)

# True
lines(grf_seq, true_1d, col = "black", lwd = 2, lty = 2)

legend("topleft", c("True", "Mean fit"), col = c("black", "purple"), lty = c(2, 1), lwd = 2)

# ============================================================================
# PLOT 2: GAM + VarProp
# ============================================================================

plot(grf_seq, true_1d, type = "n", 
     main = "GAM + VarProp: 500 Simulations",
     xlab = "GRF Covariate (scaled)", ylab = "log Density (centered)",
     ylim = c(min(pred_varprop_1d), max(pred_varprop_1d)))

# Spaghetti lines
for(i in 1:500) {
  lines(grf_seq, pred_varprop_1d[i, ], col = rgb(0.5, 0.5, 0.5, 0.1))
}

# Mean
lines(grf_seq, colMeans(pred_varprop_1d), col = "purple", lwd = 2)

legend("topleft", c("Mean fit"), col = c ("purple"), lty = c(1), lwd = 2)