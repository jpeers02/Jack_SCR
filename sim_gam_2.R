# ============================================================================
# GAM SIMULATION — FULL RE-RUN
# Captures EVERYTHING needed for proper coverage analysis:
#   - point estimates (log-density, l_esas = 0) at 1D / grid / session levels
#   - SEs (se.fit) at every level for BOTH gam and gam+varprop
#   - per-sim ESAs
#   - timing
# Coverage logic mirrors the basic/ME sims: per-sim CI = est +/- 1.96*se,
# checked against the TRUE log-density, then averaged across sims.
# ============================================================================

library(secr)
library(mgcv)
library(geoR)
library(dplyr)
source("SCR_func.R")   # sim_adapt_data_multi_2, density_fit, etc.

# ----------------------------------------------------------------------------
# SETUP: grid, GRF, sessions, traps  (unchanged from sim_GRF.R)
# ----------------------------------------------------------------------------

x.range <- c(0, 100000)
y.range <- c(0, 100000)
xs <- seq(x.range[1], x.range[2], length.out = 100)
ys <- seq(y.range[1], y.range[2], length.out = 100)
xx <- rep(xs, length(ys))
yy <- rep(ys, each = length(xs))
mask.full <- cbind(xx, yy)

pars <- c(0.95, 1500, 50)
sess.xs <- seq(x.range[1] + 3*pars[2], x.range[2] - 3*pars[2], length.out = 6)
sess.ys <- seq(y.range[1] + 3*pars[2], y.range[2] - 3*pars[2], length.out = 6)
sess.locs <- expand.grid(sess.xs, sess.ys)
n.sessions <- nrow(sess.locs)

traps_list_grf <- vector(mode = "list", length = n.sessions)
for(i in 1:n.sessions){
  trap_grid <- expand.grid(
    x = sess.locs[i, 1] + c(-1000, 0, 1000),
    y = sess.locs[i, 2] + c(-1000, 0, 1000)
  )
  traps_list_grf[[i]] <- as.data.frame(trap_grid)
}

# Fixed GRF (same surface every run)
set.seed(506)
grf_sim <- grf(nrow(mask.full), grid = mask.full,
               cov.model = "matern", cov.pars = c(2, 100000 * 1.2), kappa = 1.5)
grf_vals <- grf_sim$data

nearest <- apply(sess.locs, 1, function(pt){
  dists <- sqrt((mask.full[,1] - pt[1])^2 + (mask.full[,2] - pt[2])^2)
  which.min(dists)
})

covs_session_grf <- data.frame(
  grf = as.numeric(scale(grf_vals[nearest])),
  x.s = as.numeric(scale(sess.locs[,1])),
  y.s = as.numeric(scale(sess.locs[,2]))
)

D.calc <- function(grf_val){
  b0 <- -log(44426388/200000)
  alpha <- 1; gamma <- 3; tau <- 0
  exp(b0 + alpha * (1 - exp(-gamma * (grf_val - tau)))) + 0.005
}
D_session_grf <- D.calc(covs_session_grf$grf)

# GRF scaling params (raw -> scaled), for mapping the grid
grf_center <- mean(grf_vals[nearest])
grf_scale  <- sd(grf_vals[nearest])

# ----------------------------------------------------------------------------
# PREDICTION GRIDS (computed once)
# ----------------------------------------------------------------------------

# 1D: grf_seq is on the SCALED grf range (same units the GAM was fit on)
grf_seq <- seq(min(covs_session_grf$grf), max(covs_session_grf$grf), length.out = 100)

# 2D spatial grid
newdata_grid <- expand.grid(x = xs, y = ys)
nearest_grid <- apply(newdata_grid, 1, function(pt){
  dists <- sqrt((mask.full[,1] - pt[1])^2 + (mask.full[,2] - pt[2])^2)
  which.min(dists)
})
# map raw GRF -> scaled, matching how the GAM covariate was built
newdata_grid$grf <- (grf_vals[nearest_grid] - grf_center) / grf_scale

# TRUE log-densities (the targets for coverage)
true_log_1d   <- log(D.calc(grf_seq))                 # grf_seq already scaled
true_log_grid <- log(D.calc(newdata_grid$grf))        # newdata_grid$grf scaled above
true_log_sess <- log(D.calc(covs_session_grf$grf))    # per-session

# ----------------------------------------------------------------------------
# STORAGE
# ----------------------------------------------------------------------------

n_sims  <- 500
n_grid  <- nrow(newdata_grid)
n_sess  <- n.sessions
n_1d    <- length(grf_seq)

# 1D (log-density scale): estimate + SE
pred_1d_gam <- matrix(NA, n_sims, n_1d);  se_1d_gam <- matrix(NA, n_sims, n_1d)
pred_1d_vp  <- matrix(NA, n_sims, n_1d);  se_1d_vp  <- matrix(NA, n_sims, n_1d)

# 2D grid (log-density scale): estimate + SE
pred_grid_gam <- matrix(NA, n_sims, n_grid);  se_grid_gam <- matrix(NA, n_sims, n_grid)
pred_grid_vp  <- matrix(NA, n_sims, n_grid);  se_grid_vp  <- matrix(NA, n_sims, n_grid)

# Session level (log-density scale): estimate + SE
pred_sess_gam <- matrix(NA, n_sims, n_sess);  se_sess_gam <- matrix(NA, n_sims, n_sess)
pred_sess_vp  <- matrix(NA, n_sims, n_sess);  se_sess_vp  <- matrix(NA, n_sims, n_sess)

# ESAs (per session, per sim)
esa_mat <- matrix(NA, n_sims, n_sess)

# Counts (per session, per sim) — handy for diagnostics
n_mat <- matrix(NA, n_sims, n_sess)

timing_df <- data.frame()

# ----------------------------------------------------------------------------
# MAIN LOOP
# ----------------------------------------------------------------------------

secr.traps <- lapply(traps_list_grf, function(df){
  colnames(df) <- c("x", "y"); read.traps(data = df, detector = "proximity")
})

for(i in 1:n_sims){
  
  ok <- tryCatch({
    
    # ---- Simulate data (fixed GRF, new captures each sim) ----
    sim_data_i <- sim_adapt_data_multi_2(
      beta0 = NULL, betas = c(),
      covs_session = as.matrix(covs_session_grf[, "grf", drop = FALSE]),
      traps_list = traps_list_grf,
      g0 = 0.95, alpha0 = log(1500), alphas = c(0),
      detfn = 'hn', model = 'inh',
      xlim = x.range, ylim = y.range,
      D_session = D_session_grf
    )
    
    capt_i  <- sim_data_i$capthist_df
    secr_c  <- make.capthist(capt_i, secr.traps, fmt = 'trapID')
    mask    <- make.mask(secr.traps, buffer = 1500 * 5)
    
    # ---- Stage 1: secr (conditional likelihood) ----
    s_fit_i <- secr.fit(secr_c, mask = mask, CL = TRUE, trace = FALSE)
    
    # ---- Extract ESA / n / Z ----
    glmmfit_i <- density_fit(s_fit_i, 'grf', sim_data_i$session.cov,
                             traps_list = traps_list_grf)
    n_i     <- glmmfit_i$n
    esas_i  <- glmmfit_i$esa_s
    l_esas_i<- log(esas_i)
    Z_i     <- glmmfit_i$Z
    
    df_gam_i <- data.frame(esas = esas_i, l_esas = l_esas_i,
                           covs_session_grf, n = n_i)
    df_gam_i$XX <- Z_i
    
    # ---- Stage 2: GAM (no varprop) ----
    t_gam <- system.time({
      fit_gam_i <- gam(n ~ s(grf, k = 8) + offset(l_esas),
                       family = poisson, data = df_gam_i, method = "REML")
    })["elapsed"]
    
    # ---- Stage 2: GAM + variance propagation ----
    ESA_hess_i <- solve(vcov(s_fit_i, type = "linked"))
    t_vp <- system.time({
      fit_vp_i <- gam(n ~ s(grf, k = 8) + offset(l_esas) + XX,
                      paraPen = list(XX = list(ESA_hess_i)),
                      family = poisson, data = df_gam_i, method = "REML")
    })["elapsed"]
    
    # ======================================================================
    # PREDICTIONS — link scale, l_esas = 0  =>  linear predictor = log density
    # se.fit = TRUE everywhere for proper per-sim coverage
    # ======================================================================
    
    XX_zero_1d   <- matrix(0, n_1d,   ncol(Z_i))
    XX_zero_grid <- matrix(0, n_grid, ncol(Z_i))
    XX_zero_sess <- matrix(0, n_sess, ncol(Z_i))
    
    # --- 1D ---
    nd_1d <- data.frame(grf = grf_seq, l_esas = 0); nd_1d$XX <- XX_zero_1d
    p_gam <- predict(fit_gam_i, newdata = nd_1d, type = "link", se.fit = TRUE)
    p_vp  <- predict(fit_vp_i,  newdata = nd_1d, type = "link", se.fit = TRUE, exclude = "XX")
    pred_1d_gam[i,] <- p_gam$fit; se_1d_gam[i,] <- p_gam$se.fit
    pred_1d_vp[i,]  <- p_vp$fit;  se_1d_vp[i,]  <- p_vp$se.fit
    
    # --- 2D grid ---
    nd_grid <- newdata_grid; nd_grid$l_esas <- 0; nd_grid$XX <- XX_zero_grid
    g_gam <- predict(fit_gam_i, newdata = nd_grid, type = "link", se.fit = TRUE)
    g_vp  <- predict(fit_vp_i,  newdata = nd_grid, type = "link", se.fit = TRUE, exclude = "XX")
    pred_grid_gam[i,] <- g_gam$fit; se_grid_gam[i,] <- g_gam$se.fit
    pred_grid_vp[i,]  <- g_vp$fit;  se_grid_vp[i,]  <- g_vp$se.fit
    
    # --- Session level ---
    nd_sess <- data.frame(grf = covs_session_grf$grf, l_esas = 0); nd_sess$XX <- XX_zero_sess
    s_gam <- predict(fit_gam_i, newdata = nd_sess, type = "link", se.fit = TRUE)
    s_vp  <- predict(fit_vp_i,  newdata = nd_sess, type = "link", se.fit = TRUE, exclude = "XX")
    pred_sess_gam[i,] <- s_gam$fit; se_sess_gam[i,] <- s_gam$se.fit
    pred_sess_vp[i,]  <- s_vp$fit;  se_sess_vp[i,]  <- s_vp$se.fit
    
    # --- ESAs, counts, timing ---
    esa_mat[i,] <- esas_i
    n_mat[i,]   <- n_i
    timing_df <- rbind(timing_df, data.frame(
      sim = i, method = c("GAM", "GAM_VarProp"),
      time_sec = c(as.numeric(t_gam), as.numeric(t_vp))
    ))
    
    TRUE
  }, error = function(e){
    cat("  sim", i, "failed:", conditionMessage(e), "\n"); FALSE
  })
  
  cat('Completed', i, 'of', n_sims, if(!isTRUE(ok)) '(FAILED)' else '', '\n')
}

# ----------------------------------------------------------------------------
# SAVE EVERYTHING
# ----------------------------------------------------------------------------

saveRDS(list(
  # grids / metadata
  grf_seq        = grf_seq,
  newdata_grid   = newdata_grid,
  sess.locs      = sess.locs,
  covs_session_grf = covs_session_grf,
  grf_center     = grf_center,
  grf_scale      = grf_scale,
  
  # truth (log-density scale)
  true_log_1d    = true_log_1d,
  true_log_grid  = true_log_grid,
  true_log_sess  = true_log_sess,
  D_session_grf  = D_session_grf,
  
  # 1D estimates + SEs
  pred_1d_gam = pred_1d_gam, se_1d_gam = se_1d_gam,
  pred_1d_vp  = pred_1d_vp,  se_1d_vp  = se_1d_vp,
  
  # grid estimates + SEs
  pred_grid_gam = pred_grid_gam, se_grid_gam = se_grid_gam,
  pred_grid_vp  = pred_grid_vp,  se_grid_vp  = se_grid_vp,
  
  # session estimates + SEs
  pred_sess_gam = pred_sess_gam, se_sess_gam = se_sess_gam,
  pred_sess_vp  = pred_sess_vp,  se_sess_vp  = se_sess_vp,
  
  # ESAs, counts, timing
  esa_mat   = esa_mat,
  n_mat     = n_mat,
  timing_df = timing_df
), "sim_grf_500_gam_FULL2.rds")

cat("\nDone. Saved to sim_grf_500_gam_FULL2.rds\n")