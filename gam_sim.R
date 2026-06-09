grf_seq <- seq(min(covs_session_grf$grf), max(covs_session_grf$grf), length.out = 100)

# 2D: Spatial grid for heatmaps
newdata_grid <- expand.grid(
  x = seq(x.range[1], x.range[2], length.out = 100),
  y = seq(y.range[1], y.range[2], length.out = 100)
)

# Assign GRF values to grid points via nearest mask point
nearest_grid <- apply(newdata_grid, 1, function(pt){
  dists <- sqrt((mask.full[,1] - pt[1])^2 + (mask.full[,2] - pt[2])^2)
  which.min(dists)
})

# GRF scaling parameters (from original sessions)
nearest_sess <- apply(sess.locs, 1, function(pt){
  dists <- sqrt((mask.full[,1] - pt[1])^2 + (mask.full[,2] - pt[2])^2)
  which.min(dists)
})
grf_center <- mean(grf_vals[nearest_sess])
grf_scale <- sd(grf_vals[nearest_sess])

# Apply scaling to grid
newdata_grid$grf <- (grf_vals[nearest_grid] - grf_center) / grf_scale

# ============================================================================
# STORAGE: Pre-allocate matrices
# ============================================================================

n_sims <- 500

# 1D predictions (for spaghetti plot)
pred_gam_1d <- matrix(NA, n_sims, length(grf_seq))
pred_varprop_1d <- matrix(NA, n_sims, length(grf_seq))

# 2D fitted surfaces (at session locations, 36 per sim)
fitted_gam_2d <- list()
fitted_varprop_2d <- list()

# 2D grid predictions (for heatmaps, 10000 points per sim)
pred_grid_gam <- matrix(NA, n_sims, nrow(newdata_grid))
pred_grid_varprop <- matrix(NA, n_sims, nrow(newdata_grid))

# Timing
timing_df <- data.frame()

# ============================================================================
# MAIN LOOP: Simulate, fit, predict
# ============================================================================

for(i in 1:n_sims) {
  
  # ---- STEP 1: Simulate capture data (uses next RNG state) ----
  sim_data_i <- sim_adapt_data_multi_2(
    beta0 = NULL, 
    betas = c(),
    covs_session = as.matrix(covs_session_grf[, "grf", drop = FALSE]),
    traps_list = traps_list_grf,
    g0 = 0.95, 
    alpha0 = log(1500), 
    alphas = c(0),
    detfn = 'hn', 
    model = 'inh',
    xlim = x.range, 
    ylim = y.range,
    D_session = D_session_grf
  )
  
  # ---- STEP 2: Prepare secr objects ----
  secr.traps <- lapply(traps_list_grf, function(df) {
    colnames(df) <- c("x", "y")
    read.traps(data = df, detector = "proximity")
  })
  
  capt_i <- sim_data_i$capthist_df
  secr_c <- make.capthist(capt_i, secr.traps, fmt = 'trapID')
  mask <- make.mask(secr.traps, buffer = 1500 * 5)
  
  # ---- STEP 3: Fit secr stage 1 (CL=TRUE) ----
  s_fit_i <- secr.fit(secr_c, mask = mask, CL = TRUE)
  
  # ---- STEP 4: Extract ESA, counts, covariate matrix ----
  glmmfit_i <- density_fit(s_fit_i, 'grf', sim_data_i$session.cov, 
                           traps_list = traps_list_grf)
  n_i <- glmmfit_i$n
  esas_i <- glmmfit_i$esa_s
  l_esas_i <- log(esas_i)
  Z_i <- glmmfit_i$Z
  
  # ---- STEP 5: Prepare data for GAM fitting ----
  df_gam_i <- data.frame(
    esas = esas_i, 
    l_esas = l_esas_i, 
    covs_session_grf, 
    n = n_i
  )
  df_gam_i$XX <- Z_i
  
  # ---- STEP 6: Fit GAM WITHOUT variance propagation ----
  time_gam <- system.time({
    fit_gam_i <- gam(
      n ~ s(grf, k = 8) + offset(l_esas),
      family = poisson, 
      data = df_gam_i,
      method = "REML"
    )
  })["elapsed"]
  
  # ---- STEP 7: Fit GAM WITH variance propagation ----
  ESA_hessian_i <- solve(vcov(s_fit_i, type = "linked"))
  time_varprop <- system.time({
    fit_varprop_i <- gam(
      n ~ s(grf, k = 8) + offset(l_esas) + XX,
      paraPen = list(XX = list(ESA_hessian_i)),
      family = poisson,
      data = df_gam_i,
      method = "REML"
    )
  })["elapsed"]
  
  # ---- STEP 8: 1D predictions (spaghetti plot) ----
  newdata_1d <- data.frame(
    grf = grf_seq, 
    l_esas = mean(l_esas_i)
  )
  newdata_1d$XX <- matrix(0, nrow(newdata_1d), ncol(Z_i))
  pred_gam_1d[i, ] <- predict(fit_gam_i, 
                              newdata = newdata_1d,
                              type = "link")
  pred_varprop_1d[i, ] <- predict(fit_varprop_i,
                                  newdata = newdata_1d,
                                  type = "link",
                                  exclude = "XX")
  
  # ---- STEP 9: 2D fitted surfaces (at session locations) ----
  fitted_gam_2d[[i]] <- predict(fit_gam_i, type = "response")
  fitted_varprop_2d[[i]] <- predict(fit_varprop_i, 
                                    type = "response", 
                                    exclude = "XX")
  
  # ---- STEP 10: 2D grid predictions (for heatmaps) ----
  newdata_grid_i <- newdata_grid
  newdata_grid_i$l_esas <- mean(l_esas_i)
  newdata_grid_i$XX <- matrix(0, nrow(newdata_grid_i), ncol(Z_i))
  
  # Predict smooth terms (not response, to get log-density)
  pred_smooth_gam <- predict(fit_gam_i, 
                             newdata = newdata_grid_i, 
                             type = "terms")
  pred_smooth_varprop <- predict(fit_varprop_i,
                                 newdata = newdata_grid_i, 
                                 type = "terms",
                                 exclude = "XX")
  
  pred_grid_gam[i, ] <- pred_smooth_gam[, 1]
  pred_grid_varprop[i, ] <- pred_smooth_varprop[, 1]
  
  # ---- STEP 11: Store timing ----
  timing_df <- rbind(timing_df, 
                     data.frame(
                       sim = i, 
                       method = c("GAM", "GAM_VarProp"),
                       time_sec = c(time_gam, time_varprop)
                     ))
  
  cat('Completed', i, 'of', n_sims, '\n')
}

# ============================================================================
# AGGREGATE & SAVE
# ============================================================================

# Back-transform 1D predictions to density scale
mean_l_esas <- mean(log(esas_i))  # use last sim's mean for scaling
pred_gam_1d_density <- exp(pred_gam_1d - mean_l_esas)
pred_varprop_1d_density <- exp(pred_varprop_1d - mean_l_esas)

# Aggregate 2D fitted surfaces (at session locations)
fitted_gam_2d_mat <- do.call(rbind, fitted_gam_2d)  # 500 x 36
fitted_varprop_2d_mat <- do.call(rbind, fitted_varprop_2d)

mean_gam_2d <- colMeans(fitted_gam_2d_mat)
ci_lower_gam_2d <- apply(fitted_gam_2d_mat, 2, quantile, 0.025)
ci_upper_gam_2d <- apply(fitted_gam_2d_mat, 2, quantile, 0.975)

mean_varprop_2d <- colMeans(fitted_varprop_2d_mat)
ci_lower_varprop_2d <- apply(fitted_varprop_2d_mat, 2, quantile, 0.025)
ci_upper_varprop_2d <- apply(fitted_varprop_2d_mat, 2, quantile, 0.975)

# Aggregate 2D grid predictions (for heatmaps)
mean_grid_gam <- colMeans(pred_grid_gam)
ci_lower_grid_gam <- apply(pred_grid_gam, 2, quantile, 0.025)
ci_upper_grid_gam <- apply(pred_grid_gam, 2, quantile, 0.975)

mean_grid_varprop <- colMeans(pred_grid_varprop)
ci_lower_grid_varprop <- apply(pred_grid_varprop, 2, quantile, 0.025)
ci_upper_grid_varprop <- apply(pred_grid_varprop, 2, quantile, 0.975)

# Save everything
saveRDS(list(
  # 1D predictions
  pred_gam_1d_density = pred_gam_1d_density,
  pred_varprop_1d_density = pred_varprop_1d_density,
  grf_seq = grf_seq,
  
  # 2D session-level aggregates
  fitted_gam_2d_mat = fitted_gam_2d_mat,
  fitted_varprop_2d_mat = fitted_varprop_2d_mat,
  mean_gam_2d = mean_gam_2d,
  ci_lower_gam_2d = ci_lower_gam_2d,
  ci_upper_gam_2d = ci_upper_gam_2d,
  mean_varprop_2d = mean_varprop_2d,
  ci_lower_varprop_2d = ci_lower_varprop_2d,
  ci_upper_varprop_2d = ci_upper_varprop_2d,
  sess.locs = sess.locs,
  
  # 2D grid aggregates (for heatmaps)
  pred_grid_gam = pred_grid_gam,
  pred_grid_varprop = pred_grid_varprop,
  mean_grid_gam = mean_grid_gam,
  ci_lower_grid_gam = ci_lower_grid_gam,
  ci_upper_grid_gam = ci_upper_grid_gam,
  mean_grid_varprop = mean_grid_varprop,
  ci_lower_grid_varprop = ci_lower_grid_varprop,
  ci_upper_grid_varprop = ci_upper_grid_varprop,
  newdata_grid = newdata_grid,
  
  # Metadata
  true_D_vals = D_session_grf,
  covs_session_grf = covs_session_grf,
  timing_df = timing_df
), "sim_grf_500_gam_FULL.rds")

cat("Done! Saved to sim_grf_500_gam_FULL.rds\n")
cat("Mean time per sim (GAM):", mean(timing_df$time_sec[timing_df$method == "GAM"]), "sec\n")
cat("Mean time per sim (GAM+VarProp):", mean(timing_df$time_sec[timing_df$method == "GAM_VarProp"]), "sec\n")
