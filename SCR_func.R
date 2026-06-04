library(spatstat)
library(secr)
library(acre)
library(numDeriv)
library(glmmTMB)
library(dplyr)

`%||%` <- function(x, y) if (!is.null(x)) x else y

f <- function(pars, sessnum){
  log(esa(test.fit.secr.cl, sessnum = sessnum, beta = pars)[1])
}
# add data set to github

# Function to compute detection probabilities.

compute_p <- function(lambda0 = NULL, g0 = NULL, dists, sigma, detfn = c('hn', 'hhn')){
  detfn <- match.arg(detfn)
  if(detfn == "hn"){
    p <- g0 * exp(-dists^2 / (2 * sigma^2))
  } else if(detfn == "hhn"){
    lamd <- lambda0 * exp(-dists^2 / (2 * sigma^2))
    p <- 1 - exp(-lamd)
  }
}

## Simulating a single session & occasion data set

sim_scr_data <- function(D, lambda0 = NULL, g0 = NULL, 
                         sigma, xlim, ylim, 
                         traps, detfn = c('hn','hhn')) {
  detfn <- match.arg(detfn)
  a_m2 <- diff(xlim) * diff(ylim)
  a_ha <- a_m2 / 10000
  N <- rpois(1, D * a_ha)
  
  activity_centers <- cbind(runif(N, xlim[1], xlim[2]),
                            runif(N, ylim[1], ylim[2]))
  
  dists <- crossdist(activity_centers[, 1], activity_centers[, 2],
                     traps[, 1], traps[, 2])
  
  p <- compute_p(dists = dists, sigma = sigma, g0 = g0, lambda0 = lambda0, detfn = detfn)
  
  capt <- matrix(rbinom(length(p), 1, p), nrow = N)
  capt_filt <- capt[rowSums(capt) > 0, ]
  list(capt = capt_filt, N = N, activity_centers = activity_centers)
}

sim_scr_data_multi <- function(D, g0 = NULL, sigma, xlim, ylim, traps_list, 
                               detfn = c('hn', 'hhn'), lambda0 = NULL) {
  
  detfn <- match.arg(detfn)
  
  n_sessions <- length(traps_list)
  
  # Allow scalar D replicated for all sessions
  if (length(D) == 1) {
    D <- rep(D, n_sessions)
  } else if (length(D) != n_sessions) {
    stop("Length of D must be 1 or equal to number of sessions (length(traps_list)).")
  }
  
  a_m2 <- diff(xlim) * diff(ylim)
  a_ha <- a_m2 / 10000
  
  capt_sessions <- vector("list", n_sessions)
  
  for (s in seq_len(n_sessions)) {

    N <- rpois(1, D[s] * a_ha)
    
    # simulate activity centers
    activity_centers <- cbind(
      runif(N, xlim[1], xlim[2]),
      runif(N, ylim[1], ylim[2])
    )
    
    traps <- traps_list[[s]]
    n_traps <- nrow(traps)
    
    dists <- crossdist(activity_centers[,1], activity_centers[,2],
                       traps[,1], traps[,2])
    
    p <- compute_p(dists = dists, sigma = sigma, g0 = g0, 
                   lambda0 = lambda0, detfn = detfn)
    
    capt_mat <- matrix(rbinom(N * n_traps, 1, p), nrow = N)
    
    detected <- rowSums(capt_mat) > 0
    capt_mat <- capt_mat[detected, , drop = FALSE]
    
    capt_sessions[[s]] <- list(
      capt = capt_mat,
      traps = traps,
      detected = sum(detected),
      N_session = N,
      activity_centers = activity_centers
    )
  }
  
  list(
    N_total = sum(sapply(capt_sessions, function(x) x$N_session)),
    sessions = capt_sessions
  )
}

sim_scr_data_multi_covsession <- function(
    beta0 = NULL,                       
    betas = c(),                
    covs_session = NULL,  
    alpha0, alphas = c(), 
    g0 = NULL,
    xlim, ylim,
    traps_list,
    detfn = c('hn', 'hhn'),
    lambda0 = NULL, D_session = NULL 
) {
  detfn <- match.arg(detfn)
  n_sessions <- length(traps_list)
  
  # coerce to matrix so indexing is always consistent
  if (is.null(covs_session)){
    covs_session <- matrix(0, nrow = n_sessions, ncol = length(betas))
  } else {
    covs_session <- as.matrix(covs_session)
  }
  
  if (nrow(covs_session) != n_sessions)
    stop("nrow(covs_session) must match number of sessions.")
  #if (ncol(covs_session) != length(betas))
    #stop("ncol(covs_session) must match length(betas).")
  
  if (length(g0) == 1 || is.null(g0)) g0 <- rep(g0, n_sessions)
  
  a_m2 <- diff(xlim) * diff(ylim)
  a_ha <- a_m2 / 10000
  
  capt_sessions <- vector("list", n_sessions)
  N_by_session <- integer(n_sessions)
  
  for (s in seq_len(n_sessions)) { 
    
    if(is.null(D_session)){
      log_D <- beta0 + sum(betas * covs_session[s, ])
      D_s <- exp(log_D)
    } else {
      D_s <- D_session[s]  # use provided session-specific density
    }
    
    N <- rpois(1, D_s * a_ha)
    N_by_session[s] <- N
    
    log_sigma <- alpha0 + sum(alphas * covs_session[s, ])
    
    sigma <- exp(log_sigma)
    
    activity_centers <- cbind(
      runif(N, xlim[1], xlim[2]),
      runif(N, ylim[1], ylim[2])
    )
    
    traps <- traps_list[[s]]
    n_traps <- nrow(traps)
    
    dists <- crossdist(activity_centers[,1], activity_centers[,2],
                       traps[,1], traps[,2])
    
    p <- compute_p(dists = dists, sigma = sigma,
                   g0 = g0[s], lambda0 = lambda0, detfn = detfn)
    
    capt_mat <- matrix(rbinom(N * n_traps, 1, p), nrow = N)
    detected <- rowSums(capt_mat) > 0
    capt_mat <- capt_mat[detected, , drop = FALSE]
    
    capt_sessions[[s]] <- list(
      capt = capt_mat,
      traps = traps,
      detected = sum(detected),
      N_session = N,
      D_session = D_s,
      covs_session = covs_session[s, ],  # store full row
      activity_centers = activity_centers
    )
    
    cat('Completed', s, 'of', n_sessions, 'sessions')
  }
  
  list(
    N_total = sum(N_by_session),
    N_by_session = N_by_session,
    sessions = capt_sessions
  )
}
# Need session covariate data frame ^

## Adapting the data set into  format required for secr and acre

# Adapt to allow for multi-session

sim_adapt_data <- function(D,lambda0 = NULL, g0 = NULL, sigma, xlim, ylim, traps, session_id = 1, detfn) {
  
  
  sim <- sim_scr_data(D, lambda0, g0, sigma, xlim, ylim, traps, detfn)
  
  capt <- sim$capt
  capthist_df <- data.frame()
  for (i in 1:nrow(capt)) {
    for (j in 1:ncol(capt)) {
      if (capt[i, j] > 0) {
        capthist_df <- rbind(capthist_df, 
                             data.frame(session = session_id,
                                        ID = paste0("ind", i),
                                        occasion = 1,
                                        stringsAsFactors = FALSE, trap = j))
      }
    }
  }
  
  return(list(sim_data = sim, capthist_df = capthist_df))
}

sim_adapt_data_multi <- function(D = NULL, lambda0 = NULL, g0 = NULL,
                                 xlim, ylim, traps_list, detfn = c("hn", 'hhn'),
                                 model = c('h', 'inh'), beta0 = NULL, betas = NULL, covs_session = NULL, alpha0, alphas, 
                                 D_session = NULL) {
  
  detfn <- match.arg(detfn)
  model <- match.arg(model)
  
  if (model == 'h') {
    sim_multi <- sim_scr_data_multi(
      D = D, g0 = g0, sigma = sigma,
      xlim = xlim, ylim = ylim,
      traps_list = traps_list, detfn = detfn
    )
    session.cov <- NULL
  }
  
  if (model == 'inh') {
    sim_multi <- sim_scr_data_multi_covsession(
      beta0 = beta0, betas = betas, covs_session = covs_session, g0 = g0, alpha0 = alpha0, alphas = alphas, 
      xlim = xlim, ylim = ylim,
      detfn = detfn, traps_list = traps_list, D_session = D_session
    )
    session.cov <- as.data.frame(covs_session)
    session.cov$session <- seq_along(traps_list)
  }
  
  all_capthist <- data.frame()
  
  for (s in seq_along(sim_multi$sessions)) {
    capt <- sim_multi$sessions[[s]]$capt
    if (nrow(capt) > 0) {
      for (i in 1:nrow(capt)) {
        for (j in 1:ncol(capt)) {
          if (capt[i, j] > 0) {
            all_capthist <- rbind(
              all_capthist,
              data.frame(
                session = s,
                ID = paste0("s", s, "_", i),
                occasion = 1,
                trap = j,
                stringsAsFactors = FALSE
              )
            )
          }
        }
      }
    }
    
    cat('Completed', s, 'of', length(sim_multi$sessions), 'sims')
  }
  
  return(list(
    sim_data = sim_multi,
    capthist_df = all_capthist,
    session.cov = session.cov
  ))
}

sim_adapt_data_multi_2 <- function(D = NULL, lambda0 = NULL, g0 = NULL,
                                   xlim, ylim, traps_list, detfn = c("hn", 'hhn'),
                                   model = c('h', 'inh'), beta0 = NULL, betas = NULL, 
                                   covs_session = NULL, alpha0, alphas, 
                                   D_session = NULL) {
  
  detfn <- match.arg(detfn)
  model <- match.arg(model)
  
  if(model == 'h'){
    sim_multi <- sim_scr_data_multi(
      D = D, g0 = g0, sigma = sigma,
      xlim = xlim, ylim = ylim,
      traps_list = traps_list, detfn = detfn
    )
    session.cov <- NULL
  }
  
  if(model == 'inh'){
    sim_multi <- sim_scr_data_multi_covsession(
      beta0 = beta0, betas = betas, covs_session = covs_session, 
      g0 = g0, alpha0 = alpha0, alphas = alphas, 
      xlim = xlim, ylim = ylim,
      detfn = detfn, traps_list = traps_list, D_session = D_session
    )
    session.cov <- as.data.frame(covs_session)
    session.cov$session <- seq_along(traps_list)
  }
  
  all_capthist <- do.call(rbind, lapply(seq_along(sim_multi$sessions), function(s){
    capt <- sim_multi$sessions[[s]]$capt
    if(nrow(capt) == 0) return(NULL)
    
    idx <- which(capt > 0, arr.ind = TRUE)
    if(nrow(idx) == 0) return(NULL)
    
    data.frame(
      session = s,
      ID = paste0("s", s, "_", idx[, 1]),
      occasion = 1,
      trap = idx[, 2],
      stringsAsFactors = FALSE
    )
  }))
  
  return(list(
    sim_data = sim_multi,
    capthist_df = all_capthist,
    session.cov = session.cov
  ))
}


# Start simulating multi-session data - list of traps 

# Buffer = 5 * Sigma

# add hazard-half normal detection function simulation - good for acoustics - allow secr & acre to use.

model.sim <- function(n_sims, D, g0 = NULL, lambda0 = NULL, sigma, 
                      xlim, ylim, traps, buffer, detfn, 
                      method = c('secr', 'acre', 'both')) {
  
  
  
  method <- match.arg(method)
  results <- vector("list", n_sims)
  
  
  for (i in 1:n_sims) {
    
    if (is.data.frame(traps) || is.matrix(traps)) {

      sim_data <- sim_adapt_data(D = D, g0 = g0, lambda0 = lambda0, sigma = sigma, 
                                 xlim = xlim, ylim = ylim, traps = traps, 
                                 session_id = 1, detfn = detfn)
      
      secr.traps <- as.data.frame(traps)
      colnames(secr.traps) <- c("x", "y")
     secr_traps <- read.traps(secr.traps, detector = 'proximity')
      
    } else if (is.list(traps)) {

      sim_data <- sim_adapt_data_multi(D = D, g0 = g0, lambda0 = lambda0, sigma = sigma, 
                                       xlim = xlim, ylim = ylim, traps_list = traps, 
                                       detfn = detfn)

      secr.traps <-  lapply(traps, function(df) {
        colnames(df) <- c("x", "y")
        read.traps(data = df, detector = "proximity")
      })
      secr_traps <- secr.traps
    } else {
      stop("traps must be either a data.frame/matrix (single session) or a list of trap arrays (multi-session).")
    }
    sim_results <- list()
    

    if (method %in% c("secr", "both")) {
      
      capthist <- sim_data$capthist_df
      capt_hist <- make.capthist(capthist, secr_traps, fmt = "trapID")
      
      mask <- make.mask(secr_traps, buffer = buffer)
      fit_secr <- secr.fit(capt_hist, mask = mask)
      
      if(is.data.frame(traps)|| is.matrix(traps)){
      pred <- summary(fit_secr)$predicted
    }
      else if(is.list(traps)){
     
      # Without session covariate - homogenous density. 
        
        pred <-summary(fit_secr)$predicted$`session = 1`
      
      }
      
      df_secr <- data.frame(
        method    = "secr",
        parameter = rownames(pred),
        estimate  = pred[, "estimate"],
        se        = pred[, "SE.estimate"]
      )
      sim_results[["secr"]] <- df_secr
    }
    

    if (method %in% c("acre", "both")) {
      captu <- sim_data$capthist_df
      data_acre <- read.acre(captu, traps, control.mask = list(buffer = buffer))
      fit_acre <- fit.acre(data_acre)
      
      coef_acre <- summary(fit_acre)$coefs
      
      se_acre <- summary(fit_acre)$coefs_se |> as.vector()

      df_acre <- data.frame(
        method =  "acre",
        parameter = names(coef_acre), 
        estimate  = as.vector(coef_acre),
        se        = se_acre
      )
      
      sim_results[["acre"]] <- df_acre
    }
    
    if(detfn == "hn") {
      true_vals <- c(D = D, g0 = g0, sigma = sigma)
    } else if(detfn == "hhn") {
      true_vals <- c(D = D, lambda0 = lambda0, sigma = sigma)
    }    
    results[[i]] <- bind_rows(sim_results, .id = "fit_method") |>
      mutate(sim = i,
         true = true_vals[parameter],
         pct_diff = 100 * ((estimate - true)/true))
  }
  
  all_results <- bind_rows(results)
  
  return(list(tidy = all_results, results = results))
}



density_fit <- function(fit, formula, sessioncov, traps_list = NULL, capthist = NULL, family = 'poisson', REML = FALSE){
  
  `%||%` <- function(x, y) if (!is.null(x)) x else y
  ## add stops ie no CL = true 
  
  esa_s <- numeric(length(traps_list))
  n_pars <- NULL
  n <- NULL
  sessioncov <- as.data.frame(sessioncov)
  sessioncov$session <- factor(1:nrow(sessioncov))  # add this
  sessioncov$obs_session <- sessioncov$session
  sessioncov$session <- NULL
  if(class(fit) %in% c('acre_tmb', 'acre')){
    
    
    
    n <- fit$output.tmb$dims$n.IDs
    Z <- fit$esa_partial_derivatives
    esa_s <- acre::esa(fit)$value
    v.theta <- vcov(fit, type = "linked")
    n_pars <- ncol(Z)
    sessioncov$n <- n
    n_sessions <- nrow(sessioncov)
  }
  
  if(class(fit) == 'secr'){
    
    
    for(j in seq_along(traps_list)){
      esa_s[j] <- secr::esa(fit, sessnum = j)[1]
    }
    
    n <- sapply(fit$capthist, function(x) dim(x)[1])
    
    sessioncov$n <- n
    
    n_sessions <- nrow(sessioncov)
    n_pars <- length(coef(fit)[, 1])
    
    f <- function(pars, sessnum){
      log(secr::esa(fit, sessnum = sessnum, beta = pars)[1])
    }
    
    Z <- t(sapply(seq_len(n_sessions), function(s){
      grad(f, coef(fit)[, 1], sessnum = s)
    }))
    
    v.theta <- vcov(fit)
    
  }
    v.log.sd <- log(sqrt(diag(v.theta)))
    v.cor <- cov2cor(v.theta)[lower.tri(v.theta)]
    v.pars <- c(v.log.sd, put_cor(v.cor, input_val = "vec"))
    
    z_terms <- paste0("z", 1:n_pars, collapse = " + ")
    formula_str <- paste("n ~", formula, "+ us(0 +", z_terms, "| dummy)")
    
    # with these
    df_glmm <- data.frame(
      n = n,
      sessioncov[, !colnames(sessioncov) %in% c("n")],  # keep session
      Z,
      dummy = factor(1)
    )
    n_covs <- sum(!colnames(sessioncov) %in% c("session", "n"))
    names(df_glmm)[(n_covs + 2):(n_covs + 1 + n_pars)] <- paste0("z", 1:n_pars)
    
    has_random <- grepl("\\|", formula)
    
    if(has_random){
      v.pars_start <- c(0, v.pars)  
      map_theta <- factor(c(1, rep(NA, length(v.pars))))  
    } else {
      v.pars_start <- v.pars
      map_theta <- factor(rep(NA, length(v.pars)))
    }

    
    fit_varprop <- glmmTMB(
      as.formula(formula_str),
      offset = log(esa_s),
      start = list(theta = v.pars_start),
      map = list(theta = map_theta),
      family = family,
      data = df_glmm, 
      REML = REML
    )
  
    
    return(list(
      fit = fit_varprop,
      coefs = fixef(fit_varprop)$cond,
      se = sqrt(diag(vcov(fit_varprop)$cond)),
      esa_s = esa_s,
      n = n, 
      Z = Z
    ))
  
}


model.sim_inh <- function(n_sims, beta0 = NULL, g0 = NULL, betas = c(), covs_session,
                          lambda0 = NULL, xlim, ylim, traps_list, buffer, detfn = c('hn', 'hhn'), 
                          method = c('secr', 'acre', 'both'), stage = 1, 
                          alpha0, alphas = c(), D_formula = NULL, sigma_formula = NULL, 
                          ncores = NULL, spacing = NULL, varprop = FALSE, sigma_u = NULL, fitted = FALSE,
                          D_session = NULL
){
  
  # Room for improvement: Recording length of each model 
  # - 
  fitted_results <- vector("list", n_sims)
  method <- match.arg(method)
  detfn <- match.arg(detfn)
  results <- vector("list", n_sims)
  
  acre_model <- list()
  if(!is.null(D_formula)) acre_model$D <- as.formula(paste("~", D_formula))
  if(!is.null(sigma_formula)) acre_model$sigma <- as.formula(paste("~", sigma_formula))
  
  secr_model <- list()
  D_formula_fixed <- gsub("\\+?\\s*\\([^)]+\\|[^)]+\\)", "", D_formula)
  D_formula_fixed <- trimws(D_formula_fixed)
  if(!is.null(D_formula)) secr_model[[1]] <- as.formula(paste("D ~", D_formula_fixed))
  if(!is.null(sigma_formula)) secr_model[[2]] <- as.formula(paste("sigma ~", sigma_formula))
  
  secr_model_2 <- list()
  if(!is.null(sigma_formula)) secr_model_2[[1]] <- as.formula(paste("sigma ~", sigma_formula))
  
  
  D_names <- trimws(strsplit(D_formula_fixed, "\\+")[[1]])
  sigma_names <- trimws(strsplit(sigma_formula, "\\+")[[1]])
  
  # handle intercept only
  D_true <- if(length(D_names) == 1 && D_names == "1"){
    setNames(beta0, "D.(Intercept)")
  } else {
    setNames(c(beta0, betas), c("D.(Intercept)", paste0("D.", D_names)))
  }
  
  sigma_true <- if(length(sigma_names) == 1 && sigma_names == "1"){
    setNames(alpha0, "sigma.(Intercept)")
  } else {
    setNames(c(alpha0, alphas), c("sigma.(Intercept)", paste0("sigma.", sigma_names)))
  }
  if(detfn == "hn"){
    true_vals <- c(D_true, sigma_true, g0 = g0)
  } else if(detfn == "hhn"){
    true_vals <- c(D_true, sigma_true, lambda0 = lambda0)
  }
  
  if(!is.null(sigma_u)){
    true_vals <- c(true_vals, sigma_u = sigma_u)
  }
  
  for(i in 1:n_sims){
    
    
    if(!is.null(sigma_u)){
      u_s <- rnorm(length(traps_list), 0, sigma_u)
      D_session_i <- exp(beta0 + as.matrix(covs_session) %*% betas + u_s)
    } else if(!is.null(D_session)){
      D_session_i <- D_session
    } else {
      D_session_i <- NULL
    }
  
  
    sim_data <- sim_adapt_data_multi(beta0 = beta0, lambda0 = lambda0, g0 = g0, 
                                     alpha0 = alpha0, alphas = alphas,
                                     xlim = xlim, ylim = ylim, traps_list = traps_list, 
                                     model = 'inh', betas = betas, 
                                     covs_session = covs_session, detfn = detfn, D_session = D_session_i)
    
    captu <- sim_data$capthist_df
    sessioncov <- sim_data$session.cov
    colnames(sessioncov)[1:ncol(covs_session)] <- colnames(covs_session)
    sim_results <- list()
    
    if(method %in% c('acre', 'both')){
      data_acre <- read.acre(captu, traps_list, control.mask = list(buffer = buffer), 
                             session.cov = sessioncov)
      
      
      fit_acre <- fit.acre(data_acre, model = acre_model)
      coef_acre <- summary(fit_acre)$coefs
      se_acre <- summary(fit_acre)$coefs_se |> as.vector()
      
      sim_results[["acre"]] <- data.frame(
        method    = "acre",
        parameter = names(coef_acre), 
        estimate  = as.vector(coef_acre),
        se        = se_acre
      )
    }
    
    if(method %in% c('secr', 'both')){
      secr.traps <- lapply(traps_list, function(df) {
        colnames(df) <- c("x", "y")
        read.traps(data = df, detector = "proximity")
      })
      
      capt_hist <- make.capthist(sim_data$capthist_df, secr.traps, 
                                 fmt = "trapID")
      mask <- make.mask(secr.traps, buffer = buffer, spacing = spacing)
      
      if(stage %in% c(1, 3)){
        secr_time <- system.time( {fit_secr <- secr.fit(capt_hist, mask = mask, model = secr_model,
                                                        sessioncov = as.data.frame(covs_session), ncores = 
                                                          ncores) })["elapsed"]
        
        
        coef_secr <- summary(fit_secr)$coef
        param_names_secr <- rownames(coef_secr)
        param_names_secr <- gsub("^g0$", "g0", param_names_secr)
        param_names_secr <- gsub("^sigma$", "sigma.(Intercept)", param_names_secr)
        param_names_secr <- gsub("^D$", "D.(Intercept)", param_names_secr)
        
        df_secr_joint <- data.frame(
          method    = "secr_joint",
          parameter = param_names_secr,
          estimate  = coef_secr[, "beta"],
          se        = coef_secr[, "SE.beta"],
          time      = as.numeric(secr_time)
        )
        df_secr_joint$estimate <- ifelse(df_secr_joint$parameter == "g0",
                                         plogis(df_secr_joint$estimate),
                                         df_secr_joint$estimate)
        sim_results[["secr_joint"]] <- df_secr_joint
        
      }
      
      if(stage %in% c(2, 3)){
        esa_s <- numeric(length(traps_list))
        
        
        stage1_time <- system.time({fit_secr2 <- secr.fit(capt_hist, mask = mask, model = secr_model_2,
                                                          sessioncov = as.data.frame(covs_session), 
                                                          CL = TRUE, ncores = ncores)
        
        for(j in seq_along(traps_list)){
          esa_s[j] <- secr::esa(fit_secr2, sessnum = j)[1]
        }
        
        n <- sapply(fit_secr2$capthist, function(x) dim(x)[1])
        sessioncov$n <- n})['elapsed']
        
        glm_time <- system.time({
          df_glm <- as.data.frame(covs_session)
          df_glm$n <- n
          glm_form <- as.formula(paste("n ~", D_formula_fixed))
          glm2 <- glm(glm_form, offset = log(esa_s), family = poisson, data = df_glm)
        })['elapsed']
        
        secr_joint_time <- stage1_time + glm_time 
        
        glm_coefs <- coef(glm2)
        glm_ses <- sqrt(diag(vcov(glm2)))
        
        
        param_names_glm <- ifelse(
          names(glm_coefs) == "(Intercept)",
          "D.(Intercept)",
          paste0("D.", names(glm_coefs))
        )
        
        
        coef_secr2 <- summary(fit_secr2)$coef
        param_names_secr2 <- rownames(coef_secr2)
        param_names_secr2 <- gsub("^g0$", "g0", param_names_secr2)
        param_names_secr2 <- gsub("^sigma$", "sigma.(Intercept)", param_names_secr2)
        param_names_secr2 <- gsub("^sigma\\.", "sigma.", param_names_secr2)
        
        keep_detection <- grepl("^g0|^sigma", param_names_secr2)
        coef_secr2_det <- coef_secr2[keep_detection, ]
        
        df_detection <- data.frame(
          method    = "secr_glm",
          parameter = param_names_secr2[keep_detection],
          estimate  = coef_secr2_det[, "beta"],
          se        = coef_secr2_det[, "SE.beta"]
        ) |>
          mutate(estimate = ifelse(parameter == "g0", plogis(estimate), estimate))
        
        
        df_glm_D <- data.frame(
          method    = "secr_glm",
          parameter = param_names_glm,
          estimate  = as.vector(glm_coefs),
          se        = glm_ses
        )
        
        sim_results[["secr_glm"]] <- rbind(df_detection, df_glm_D) |>
          mutate(time = as.numeric(secr_joint_time))
        
        if(varprop == TRUE){
          time_varprop <- system.time({
            fit_varprop_result <- density_fit(fit_secr2, D_formula, 
                                              covs_session, traps_list)
          })["elapsed"]
          
          fit_varprop <- fit_varprop_result$fit
          time_varprop_2stage <- stage1_time + time_varprop
          
          varprop_coefs <- fixef(fit_varprop)$cond
          varprop_ses <- sqrt(diag(vcov(fit_varprop)$cond))
          
          param_names_varprop <- ifelse(
            names(varprop_coefs) == "(Intercept)",
            "D.(Intercept)",
            paste0("D.", names(varprop_coefs))
          )
          
          df_varprop_D <- data.frame(
            method    = "secr_varprop",
            parameter = param_names_varprop,
            estimate  = as.vector(varprop_coefs),
            se        = varprop_ses
          )
          
          if(!is.null(sigma_u)){
            re_sd <- sqrt(VarCorr(fit_varprop)$cond$obs_session[1,1])
            df_varprop_re <- data.frame(
              method    = "secr_varprop",
              parameter = "sigma_u",
              estimate  = re_sd,
              se        = NA
            )
            df_varprop_D <- rbind(df_varprop_D, df_varprop_re)
          }
          
          df_detection_varprop <- df_detection |> mutate(method = "secr_varprop")
          sim_results[["secr_varprop"]] <- rbind(df_detection_varprop, df_varprop_D) |>
            mutate(time = as.numeric(time_varprop_2stage))
        }
        
        if(fitted == TRUE){
          # GLM fitted density + SE
          glm_pred <- predict(glm2, se.fit = TRUE)
          fitted_D_glm <- exp(glm_pred$fit - log(esa_s))
          # delta method SE on density scale
          se_D_glm <- fitted_D_glm * glm_pred$se.fit
          
          # varprop fitted density + SE
          varprop_pred <- predict(fit_varprop, se.fit = TRUE)
          fitted_D_varprop <- exp(varprop_pred$fit - log(esa_s))
          se_D_varprop <- fitted_D_varprop * varprop_pred$se.fit
          
          # secr joint fitted density + SE
          secr_pred <- predict(fit_secr)
          fitted_D_joint <- sapply(secr_pred, function(x) x["D", "estimate"])
          se_D_joint <- sapply(secr_pred, function(x) x["D", "SE.estimate"])
          
            # true density - use D_session_i if available, otherwise compute from parameters
            true_D_vals <- if(!is.null(D_session_i)){
              D_session_i
            } else {
              exp(beta0 + as.matrix(covs_session) %*% betas)
            }
          
          fitted_df <- data.frame(
            session = 1:length(traps_list),
            true_D = true_D_vals,
            fitted_D_glm = fitted_D_glm,
            se_D_glm = se_D_glm,
            fitted_D_varprop = as.vector(fitted_D_varprop),
            se_D_varprop = as.vector(se_D_varprop),
            fitted_D_joint = fitted_D_joint,
            se_D_joint = se_D_joint,
            sim = i
          )
        
        sim_results[["fitted_density"]] <- fitted_df
        }
      }
    }
    
    results[[i]] <- bind_rows(sim_results[names(sim_results) != "fitted_density"]) |>
      mutate(sim = i,
             true = true_vals[parameter],
             pct_diff = 100 * ((estimate - true) / true))
    
    if(fitted == TRUE){
      fitted_results[[i]] <- sim_results[["fitted_density"]]
    }
    cat('Completed', i, 'of', n_sims, '.')
  }
  
  all_results <- bind_rows(results)
  rownames(all_results) <- 1:nrow(all_results)
  
  
  return(list(
    tidy = all_results, 
    results = results,
    fitted = if(fitted) bind_rows(fitted_results) else NULL
  ))
  
}



# Session - statistically independent capture data 
# Occasion - dependent - animals recognised between time & space - focus on single occasions for now.
# Acre - single occasion.
# Keep traps very separate to be independent. 
# Make grid of densities - changing depending on coordinate/covariate. 
