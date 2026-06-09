# 6x6 = 36 sessions
library(geoR)
library(secr)
library(CircStats)
# sourceCpp("fitting-functions/fitting-functions.cpp")
# source("fitting-functions/fitting-functions.R")
source('SCR_func.R')
library(mgcv)

## Creating a big grid of mask points.
x.range <- c(0, 100000)
y.range <- c(0, 100000)
xs <- seq(x.range[1], x.range[2], length.out = 100)
ys <- seq(y.range[1], y.range[2], length.out = 100)
xx <- rep(xs, length(ys))
yy <- rep(ys, each = length(xs))
mask.full <- cbind(xx, yy)
mask.area <- sqrt((max(x.range) - min(x.range))^2) *
  sqrt((max(y.range) - min(y.range))^2)

## Detection function parameter values. Order is g0, sigma, kappa.
pars <- c(0.95, 1500, 50)

sess.xs <- seq(x.range[1] + 3*pars[2], x.range[2] - 3*pars[2], length.out = 6)
sess.ys <- seq(y.range[1] + 3*pars[2], y.range[2] - 3*pars[2], length.out = 6)
sess.locs <- expand.grid(sess.xs, sess.ys)
n.sessions <- nrow(sess.locs)  # 36

# generating trap locations
traps_list_grf <- vector(mode = "list", length = n.sessions)
for(i in 1:n.sessions){
  traps_list_grf[[i]] <- data.frame(
    x = rep(sess.locs[i, 1], 3) + rep(c(-1000, 0, 1000), each = 3),
    y = rep(sess.locs[i, 2], 3) + rep(c(-1000, 0, 1000), each = 3)
  )
}

traps_list_grf <- vector(mode = "list", length = n.sessions)
for(i in 1:n.sessions){
  trap_grid <- expand.grid(
    x = sess.locs[i, 1] + c(-1000, 0, 1000),
    y = sess.locs[i, 2] + c(-1000, 0, 1000)
  )
  traps_list_grf[[i]] <- as.data.frame(trap_grid)
}

# simulate GRF over the full grid
set.seed(506)
grf_sim <- grf(nrow(mask.full),
               grid = mask.full,
               cov.model = "matern",
               cov.pars = c(2, 100000 * 1.2),  # from Ben's code
               kappa = 1.5)

grf_vals <- grf_sim$data

# assign GRF value to each session centroid via nearest grid point
nearest <- apply(sess.locs, 1, function(pt){
  dists <- sqrt((mask.full[,1] - pt[1])^2 + (mask.full[,2] - pt[2])^2)
  which.min(dists)
})

# session covariates
covs_session_grf <- data.frame(
  grf = as.numeric(scale(grf_vals[nearest])),
  x.s = as.numeric(scale(sess.locs[,1])),
  y.s = as.numeric(scale(sess.locs[,2]))
)

# true density as nonlinear function of GRF (from Ben's code)
D.calc <- function(grf_val){
  b0 <- log(1)
  alpha <- 1
  gamma <- 3
  tau <- -1
  exp(b0 + alpha * (1 - exp(-gamma * (grf_val - tau))))
}

D.calc <- function(grf_val){
  b0 <- -log(44426388/200000)
  alpha <- 1
  gamma <- 3
  tau <- -3.164784
  exp(b0 + alpha * (1 - exp(-gamma * (grf_val - tau))))
}

D.calc <- function(grf_val){
  b0 <- -log(44426388/200000)  # keep Ben's b0
  alpha <- 1
  gamma <- 3
  tau <- 0  # shift to centre of your distribution
  exp(b0 + alpha * (1 - exp(-gamma * (grf_val - tau)))) + 0.005
}



# check new curve
grf_seq <- seq(-3, 3, length.out = 100)
plot(grf_seq, D.calc(grf_seq), type = 'l',
     xlab = "GRF value", ylab = "True density")

D_session_grf <- D.calc(covs_session_grf$grf)

df_grf <- sim_adapt_data_multi_2(
  beta0 = NULL, betas = c(),
  covs_session = as.matrix(covs_session_grf[, "grf", drop = FALSE]),
  traps_list = traps_list_grf,
  g0 = 0.95, alpha0 = log(1500), alphas = c(0),
  detfn = 'hn', model = 'inh',
  xlim = x.range, ylim = y.range,
  D_session = D_session_grf
)

secr.traps <- lapply(traps_list_grf, function(df) {
  colnames(df) <- c("x", "y")
  read.traps(data = df, detector = "proximity")
})

df_grf$capthist_df -> capt
secr.c <- make.capthist(capt, secr.traps, fmt = 'trapID' )
mask <- make.mask(secr.traps, buffer = 1500 * 5)

s_fit <- secr.fit(secr.c, mask = mask, buffer = 1500 * 5, CL = TRUE)

glmmfit <- density_fit(s_fit, 'grf', df_grf$session.cov, traps_list = traps_list_grf)
n <- glmmfit$n
esas <- glmmfit$esa_s
l_esas <- log(esas)
Z <- glmmfit$Z
df_glm <- data.frame(esas, l_esas, covs_session_grf, n)
df_glm$XX <- Z
base_formula <- n ~ s(grf, k = 10) + offset(l_esas)
new_formula <- update(base_formula, . ~ . + XX)
ESA_hessian <- solve(vcov(s_fit, type = "linked"))

fit_gam <- gam(
  base_formula,
  family = poisson, 
  data=df_glm,
  method="REML"
)
## Variance propagatoin step.
fit_varprop_gam <- gam(
  new_formula,
  paraPen = list(XX = list(ESA_hessian)),
  family = poisson,
  data = df_glm,
  method="REML"
)
## Looking at summaries.
summary(fit_gam)
summary(fit_varprop_gam)


# parametric GLM - forces linear relationship
glm_grf <- glm(n ~ grf, offset = l_esas, family = poisson, data = df_glm)
summary(glm_grf)

# compare AIC
AIC(glm_grf, fit_gam)

# plot both fits against true D.calc curve
grf_seq <- seq(min(df_glm$grf), max(df_glm$grf), length.out = 100)

pred_glm <- exp(coef(glm_grf)[1] + coef(glm_grf)[2] * grf_seq)
pred_gam <- exp(predict(fit_gam, 
                        newdata = data.frame(grf = grf_seq, l_esas = mean(l_esas)),
                        type = "link") - mean(l_esas))
true_D <- D.calc(grf_seq)

plot_df <- data.frame(
  grf = rep(grf_seq, 3),
  D = c(pred_glm, pred_gam, true_D),
  model = rep(c("GLM", "GAM", "True"), each = 100)
)

ggplot(plot_df, aes(x = grf, y = D, colour = model, linetype = model)) +
  geom_line(linewidth = 1.2) +
  scale_colour_manual(values = c("GLM" = "#2ecc71", 
                                 "GAM" = "#9b59b6", 
                                 "True" = "black")) +
  scale_linetype_manual(values = c("GLM" = "solid",
                                   "GAM" = "solid", 
                                   "True" = "dashed")) +
  labs(x = "GRF covariate", y = "log(Density (animals/ha))",
       title = "True vs Fitted Density: GLM vs GAM") +
  theme_classic(base_family = "Times New Roman")


# Trial with x and y covariates 
# prediction grid over grf range
