library(survival)
library(MOSS)
# simulate data
source('./simulate_data_3.R')
# source('./simulate_data_91.R')
library(survtmle)
fit_survtmle <- function(T.tilde, Delta, A, W_df) {
  t_0 <- max(T.tilde)
  fit <- survtmle(ftime = T.tilde,
                  ftype = Delta,
                  trt = A,
                  adjustVars = W_df,
                  SL.trt = c("SL.mean", 'SL.glm', 'SL.gam'),
                  SL.ftime = c("SL.mean", 'SL.glm', 'SL.gam'),
                  SL.ctime = c("SL.mean", 'SL.glm', 'SL.gam'),
                  method = "hazard",
                  returnIC = TRUE,
                  verbose = FALSE
  )
  # extract cumulative incidence at each timepoint
  tpfit <- timepoints(fit, times = seq_len(t_0))
  len_groups <- as.numeric(unique(lapply(lapply(tpfit, FUN = `[[`,
                                                "est"), FUN = length)))
  names_groups <- unique(lapply(lapply(tpfit, FUN = `[[`, "est"),
                                FUN = rownames))[[1]]
  est_only <- t(matrix(unlist(lapply(tpfit, FUN = `[[`, "est")),
                       ncol = len_groups, byrow = TRUE))
  est_only <- as.data.frame(est_only)
  rownames(est_only) <- names_groups
  colnames(est_only) <- paste0("t", seq_len(ncol(est_only)))

  s_0 <- 1 - as.numeric(est_only[1,])
  s_1 <- 1 - as.numeric(est_only[2,])
  return(data.frame(time = 1:t_0, s_0 = s_0, s_1 = s_1))
}
do_once <- function(n_sim = 2e2) {
  simulated <- simulate_data(n_sim = n_sim)
  df <- simulated$dat
  true_surv <- simulated$true_surv1

  sl_lib_g <- c("SL.mean", "SL.glm", "SL.gam")
  sl_lib_censor <- c("SL.mean", "SL.glm", "SL.gam", "SL.earth")
  sl_lib_failure <- c("SL.mean", "SL.glm", "SL.gam", "SL.earth")
  range(df$T.tilde)
  df <- df[df$T.tilde > 0, ]
  k_grid <- 1:max(df$T.tilde)

  message("KM")
  n_sample <- nrow(df)
  km_fit <- survfit(Surv(time = T.tilde, event = Delta) ~ A, data = df)
  surv1_km <- tail(km_fit$surv, km_fit$strata['A=1'])
  time1_km <- tail(km_fit$time, km_fit$strata['A=1'])
  surv0_km <- tail(km_fit$surv, km_fit$strata['A=0'])
  time0_km <- tail(km_fit$time, km_fit$strata['A=0'])
  impute_KM <- function(time, km) {
    surv1_km_final <- rep(NA, max(df$T.tilde))
    surv1_km_final[time] <- km
    surv1_km_final <- na.locf(surv1_km_final, na.rm = FALSE)
    surv1_km_final[is.na(surv1_km_final)] <- 1
    return(surv1_km_final)
  }
  surv1_km_final <- impute_KM(time = time1_km, km = surv1_km)
  surv0_km_final <- impute_KM(time = time0_km, km = surv0_km)
  km_fit_1 <- survival_curve$new(t = k_grid, survival = surv1_km_final)
  km_fit_0 <- survival_curve$new(t = k_grid, survival = surv0_km_final)

  message("SL")
  sl_fit <- initial_sl_fit(
    ftime = df$T.tilde,
    ftype = df$Delta,
    trt = df$A,
    adjustVars = data.frame(df[, c("W", "W1")]),
    t_0 = max(df$T.tilde),
    SL.trt = sl_lib_g,
    SL.ctime = sl_lib_censor,
    SL.ftime = sl_lib_failure
  )
  sl_fit$density_failure_1$hazard_to_survival()
  sl_fit$density_failure_0$hazard_to_survival()
  # WILSON hack no data is t_tilde = 2
  sl_fit$density_failure_1$t <- k_grid
  sl_fit$density_failure_0$t <- k_grid

  # ipcw
  message("ipcw + ee")
  ipcw_fit_1_all <- repeat_t_grid$new(
    method = ipcw,
    A = df$A,
    T_tilde = df$T.tilde,
    Delta = df$Delta,
    density_failure = sl_fit$density_failure_1,
    density_censor = sl_fit$density_censor_1,
    g1W = sl_fit$g1W,
    A_intervene = 1
  )$fit(k_grid = k_grid)
  ipcw_fit_0_all <- repeat_t_grid$new(
    method = ipcw,
    A = df$A,
    T_tilde = df$T.tilde,
    Delta = df$Delta,
    density_failure = sl_fit$density_failure_0,
    density_censor = sl_fit$density_censor_0,
    g1W = sl_fit$g1W,
    A_intervene = 0
  )$fit(k_grid = k_grid)
  ee_fit_1_all <- repeat_t_grid$new(
    method = ee,
    A = df$A,
    T_tilde = df$T.tilde,
    Delta = df$Delta,
    density_failure = sl_fit$density_failure_1,
    density_censor = sl_fit$density_censor_1,
    g1W = sl_fit$g1W,
    A_intervene = 1
  )$fit(k_grid = k_grid)
  ee_fit_0_all <- repeat_t_grid$new(
    method = ee,
    A = df$A,
    T_tilde = df$T.tilde,
    Delta = df$Delta,
    density_failure = sl_fit$density_failure_0,
    density_censor = sl_fit$density_censor_0,
    g1W = sl_fit$g1W,
    A_intervene = 0
  )$fit(k_grid = k_grid)
  ipcw_fit_1 <- survival_curve$new(t = k_grid, survival = ipcw_fit_1_all)
  ipcw_fit_0 <- survival_curve$new(t = k_grid, survival = ipcw_fit_0_all)
  ee_fit_1 <- survival_curve$new(t = k_grid, survival = ee_fit_1_all)
  ee_fit_0 <- survival_curve$new(t = k_grid, survival = ee_fit_0_all)
  library(ggpubr)
  # gg_sl <- ggarrange(
  #   sl_fit$density_failure_1$display(type = "survival", W = df$W),
  #   sl_fit$density_failure_0$display(type = "survival", W = df$W),
  #   ncol = 2,
  #   labels = "AUTO",
  #   common.legend = TRUE,
  #   legend = "bottom"
  # )
  sl_density_failure_1_marginal <- sl_fit$density_failure_1$clone(deep = TRUE)
  sl_density_failure_0_marginal <- sl_fit$density_failure_0$clone(deep = TRUE)
  sl_density_failure_1_marginal$survival <- matrix(colMeans(sl_density_failure_1_marginal$survival), nrow = 1)
  sl_density_failure_0_marginal$survival <- matrix(colMeans(sl_density_failure_0_marginal$survival), nrow = 1)
  # gg_sl2 <- ggarrange(
  #   sl_density_failure_1_marginal$display(type = "survival"),
  #   sl_density_failure_0_marginal$display(type = "survival"),
  #   ncol = 2,
  #   labels = "AUTO",
  #   common.legend = TRUE,
  #   legend = "bottom"
  # )
  # gg_ipcw <- ggarrange(
  #   ipcw_fit_1$display(type = "survival"),
  #   ipcw_fit_0$display(type = "survival"),
  #   ncol = 2,
  #   labels = "AUTO",
  #   common.legend = TRUE,
  #   legend = "bottom"
  # )
  # gg_ee <- ggarrange(
  #   ee_fit_1$display(type = "survival"),
  #   ee_fit_0$display(type = "survival"),
  #   ncol = 2,
  #   labels = "AUTO",
  #   common.legend = TRUE,
  #   legend = "bottom"
  # )
  message("moss")
  moss_fit <- MOSS$new(
    A = df$A,
    T_tilde = df$T.tilde,
    Delta = df$Delta,
    density_failure = sl_fit$density_failure_1,
    density_censor = sl_fit$density_censor_1,
    g1W = sl_fit$g1W,
    A_intervene = 1,
    k_grid = k_grid
  )
  psi_moss_1 <- moss_fit$onestep_curve(
    epsilon = 1e-3,
    # epsilon = 1e-5,
    max_num_interation = 1e2,
    verbose = TRUE
  )
  moss_fit <- MOSS$new(
    A = df$A,
    T_tilde = df$T.tilde,
    Delta = df$Delta,
    density_failure = sl_fit$density_failure_0,
    density_censor = sl_fit$density_censor_0,
    g1W = sl_fit$g1W,
    A_intervene = 0,
    k_grid = k_grid
  )
  psi_moss_0 <- moss_fit$onestep_curve(
    epsilon = 1e-3,
    # epsilon = 1e-5,
    max_num_interation = 1e2,
    verbose = TRUE
  )
  moss_fit_1 <- survival_curve$new(t = k_grid, survival = psi_moss_1)
  moss_fit_0 <- survival_curve$new(t = k_grid, survival = psi_moss_0)
  # gg_moss <- ggarrange(
  #   moss_fit_1$display(type = "survival"),
  #   moss_fit_0$display(type = "survival"),
  #   ncol = 2,
  #   labels = "AUTO",
  #   common.legend = TRUE,
  #   legend = "bottom"
  # )
  # ggarrange(
  #   gg_sl,
  #   gg_sl2,
  #   gg_ipcw,
  #   gg_ee,
  #   gg_moss,
  #   nrow = 5,
  #   common.legend = TRUE,
  #   legend = "bottom"
  # )
  # tmle
  message("tmle")
  tmle_fit <- tryCatch(
    {
      tmle_fit <- fit_survtmle(
        T.tilde = df$T.tilde,
        Delta = df$Delta,
        A = df$A,
        W_df = data.frame(df[, c("W", "W1")])
      )
    },
    error = function(cond) {
      message("tmle error")
      NULL
  })
  if (is.null(tmle_fit)) {
    tmle_fit_1 <- sl_density_failure_1_marginal$clone(deep = TRUE)
    tmle_fit_0 <- sl_density_failure_0_marginal$clone(deep = TRUE)
  } else {
    tmle_fit_1 <- survival_curve$new(t = k_grid, survival = tmle_fit$s_1)
    tmle_fit_0 <- survival_curve$new(t = k_grid, survival = tmle_fit$s_0)
  }
  survival_truth_1 <- survival_curve$new(t = k_grid, survival = simulated$true_surv1(k_grid))
  survival_truth_0 <- survival_curve$new(t = k_grid, survival = simulated$true_surv0(k_grid))

  df_entropy_moss_1 <- evaluate_metric$new(
    survival = moss_fit_1, survival_truth = survival_truth_1
  # )$evaluate_cross_entropy()
  )$evaluate_mse()
  df_entropy_sl_1 <- evaluate_metric$new(
    survival = sl_density_failure_1_marginal, survival_truth = survival_truth_1
  # )$evaluate_cross_entropy()
  )$evaluate_mse()
  df_entropy_ipcw_1 <- evaluate_metric$new(
    survival = ipcw_fit_1, survival_truth = survival_truth_1
  # )$evaluate_cross_entropy()
  )$evaluate_mse()
  df_entropy_ee_1 <- evaluate_metric$new(
    survival = ee_fit_1, survival_truth = survival_truth_1
  # )$evaluate_cross_entropy()
  )$evaluate_mse()
  df_entropy_tmle_1 <- evaluate_metric$new(
    survival = tmle_fit_1, survival_truth = survival_truth_1
  # )$evaluate_cross_entropy()
  )$evaluate_mse()
  df_entropy_km_1 <- evaluate_metric$new(
    survival = km_fit_1, survival_truth = survival_truth_1
  # )$evaluate_cross_entropy()
  )$evaluate_mse()
  df_entropy_moss_1$method <- "MOSS"
  df_entropy_sl_1$method <- "super learner"
  df_entropy_ipcw_1$method <- "IPCW"
  df_entropy_ee_1$method <- "EE"
  df_entropy_tmle_1$method <- "TMLE"
  df_entropy_km_1$method <- "KM"
  df_plot <- rbind(
    df_entropy_moss_1,
    df_entropy_sl_1,
    df_entropy_ipcw_1,
    df_entropy_ee_1,
    df_entropy_tmle_1,
    df_entropy_km_1
  )
  return(df_plot)
}

N_SIMULATION = 1e2
# N_SIMULATION = 8
library(foreach)
# library(Rmpi)
# library(doMPI)
# cl = startMPIcluster()
# registerDoMPI(cl)
# clusterSize(cl) # just to check

library(doSNOW)
library(tcltk)
nw <- parallel:::detectCores()  # number of workers
cl <- makeSOCKcluster(nw)
registerDoSNOW(cl)

# n_sim_grid <- c(50, 1e2)
n_sim_grid <- c(1e2, 3e2)
# n_sim_grid <- c(50, 1e2, 3e2, 5e2)
df_metric <- foreach(
  n_sim = n_sim_grid,
  .combine = rbind,
  .packages = c('R6', 'MOSS', 'survtmle', 'survival'),
  .inorder = FALSE,
  .errorhandling = 'remove',
  .verbose = TRUE
) %:%
  foreach(it2 = 1:N_SIMULATION, .combine = rbind, .errorhandling = 'remove') %dopar% {
  df <- do_once(n_sim = n_sim)
  df$id_mcmc <- it2
  df$n <- n_sim
  return(df)
}
table(df_metric$id_mcmc)

save(df_metric, file = 'df_metric.rda')

# shut down for memory
# closeCluster(cl)
stopCluster(cl)
