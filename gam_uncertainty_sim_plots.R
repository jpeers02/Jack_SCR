# ============================================================================
# GAM SIM — COVERAGE & UNCERTAINTY CALIBRATION ANALYSIS
# ============================================================================
library(dplyr)
library(ggplot2)

sim_data <- readRDS("sim_grf_500_gam_FULL2.rds")

# ----------------------------------------------------------------------------
# helper: coverage = fraction of per-sim CIs containing the truth
# all on log-density scale; CIs = est +/- 1.96 * se
# ----------------------------------------------------------------------------
cov_fn <- function(pred, se, truth){
  lo <- pred - 1.96 * se
  hi <- pred + 1.96 * se
  # truth recycled across rows (sims); compare elementwise
  inside <- sweep(lo, 2, truth, "<=") & sweep(hi, 2, truth, ">=")
  inside
}

# ----------------------------------------------------------------------------
# SESSION-LEVEL (headline: 36 sessions x 500 sims)
# ----------------------------------------------------------------------------
in_gam_sess <- cov_fn(sim_data$pred_sess_gam, sim_data$se_sess_gam, sim_data$true_log_sess)
in_vp_sess  <- cov_fn(sim_data$pred_sess_vp,  sim_data$se_sess_vp,  sim_data$true_log_sess)

cat("=== SESSION-LEVEL COVERAGE ===\n")
cat(sprintf("  GAM:          %.1f%% (n = %d CIs)\n",
            100*mean(in_gam_sess, na.rm=TRUE), sum(!is.na(in_gam_sess))))
cat(sprintf("  GAM+VarProp:  %.1f%%\n", 100*mean(in_vp_sess, na.rm=TRUE)))

# per-session coverage (averaged across sims) — diagnostic
cov_by_sess_gam <- colMeans(in_gam_sess, na.rm=TRUE)
cov_by_sess_vp  <- colMeans(in_vp_sess,  na.rm=TRUE)
cat(sprintf("  GAM per-session range: %.2f - %.2f\n",
            min(cov_by_sess_gam), max(cov_by_sess_gam)))

# ----------------------------------------------------------------------------
# 1D COVERAGE (along grf gradient)
# ----------------------------------------------------------------------------
in_gam_1d <- cov_fn(sim_data$pred_1d_gam, sim_data$se_1d_gam, sim_data$true_log_1d)
in_vp_1d  <- cov_fn(sim_data$pred_1d_vp,  sim_data$se_1d_vp,  sim_data$true_log_1d)

cat("\n=== 1D COVERAGE (along GRF gradient) ===\n")
cat(sprintf("  GAM:          %.1f%%\n", 100*mean(in_gam_1d, na.rm=TRUE)))
cat(sprintf("  GAM+VarProp:  %.1f%%\n", 100*mean(in_vp_1d, na.rm=TRUE)))

# ----------------------------------------------------------------------------
# GRID COVERAGE (full spatial surface)
# ----------------------------------------------------------------------------
in_gam_grid <- cov_fn(sim_data$pred_grid_gam, sim_data$se_grid_gam, sim_data$true_log_grid)
in_vp_grid  <- cov_fn(sim_data$pred_grid_vp,  sim_data$se_grid_vp,  sim_data$true_log_grid)

cat("\n=== GRID COVERAGE (full spatial surface) ===\n")
cat(sprintf("  GAM:          %.1f%%\n", 100*mean(in_gam_grid, na.rm=TRUE)))
cat(sprintf("  GAM+VarProp:  %.1f%%\n", 100*mean(in_vp_grid, na.rm=TRUE)))

# ----------------------------------------------------------------------------
# SE COMPARISON: does varprop inflate SEs?
# ----------------------------------------------------------------------------
cat("\n=== SE COMPARISON (mean SE, session-level) ===\n")
cat(sprintf("  GAM mean SE:         %.4f\n", mean(sim_data$se_sess_gam, na.rm=TRUE)))
cat(sprintf("  GAM+VarProp mean SE: %.4f\n", mean(sim_data$se_sess_vp,  na.rm=TRUE)))
cat(sprintf("  Ratio (VP/GAM):      %.3f\n",
            mean(sim_data$se_sess_vp, na.rm=TRUE)/mean(sim_data$se_sess_gam, na.rm=TRUE)))

# ----------------------------------------------------------------------------
# ESA SANITY CHECK (should be constant within sim, vary across sims)
# ----------------------------------------------------------------------------
cat("\n=== ESA CHECK ===\n")
within_sim_sd <- apply(sim_data$esa_mat, 1, sd)  # should be ~0 (sigma~1)
across_sim_sd <- apply(sim_data$esa_mat, 2, sd)  # should be > 0
cat(sprintf("  Within-sim SD (should be ~0):   mean %.4f\n", mean(within_sim_sd, na.rm=TRUE)))
cat(sprintf("  Across-sim SD (should be > 0):  mean %.2f\n", mean(across_sim_sd, na.rm=TRUE)))

# ----------------------------------------------------------------------------
# TIMING
# ----------------------------------------------------------------------------
cat("\n=== TIMING ===\n")
timing_summary <- sim_data$timing_df %>%
  group_by(method) %>%
  summarise(mean_sec = mean(time_sec), median_sec = median(time_sec))
print(timing_summary)

# ----------------------------------------------------------------------------
# PLOT: pointwise coverage along the GRF gradient
# ----------------------------------------------------------------------------
cov_1d_df <- data.frame(
  grf = rep(sim_data$grf_seq, 2),
  coverage = c(colMeans(in_gam_1d, na.rm=TRUE), colMeans(in_vp_1d, na.rm=TRUE)),
  method = rep(c("GAM", "GAM + VarProp"), each = length(sim_data$grf_seq))
)

p_cov <- ggplot(cov_1d_df, aes(x = grf, y = coverage, colour = method)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0.95, linetype = "dashed", colour = "black") +
  scale_colour_manual(values = c("GAM" = "#9b59b6", "GAM + VarProp" = "#2ecc71")) +
  coord_cartesian(ylim = c(0.7, 1)) +
  labs(x = "GRF covariate (scaled)", y = "Pointwise coverage",
       title = "Pointwise 95% CI Coverage Along the Density Gradient") +
  theme_classic(base_family = "Times New Roman")

ggsave("gam_coverage_1d.png", p_cov, width = 8, height = 5, dpi = 150)
cat("\nSaved gam_coverage_1d.png\n")