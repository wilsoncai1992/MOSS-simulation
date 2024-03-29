assertthat::assert_that(packageVersion("MOSS") >= "1.1.2")
library(survival)
library(MOSS)
library(ggpubr)
library(tidyverse)
library(here)
source(here("fit_survtmle.R"))
# simulate data
source(here("./code_simulation/simulate_data.R"))

do_once <- function(n_sim = 2e2) {
  simulated <- simulate_data(n_sim = n_sim)
  df <- simulated$dat
  true_surv <- simulated$true_surv1
  W_names <- c("W", "W1")

  sl_lib_g <- c("SL.mean", "SL.glm", "SL.gam", "SL.earth")
  sl_lib_censor <- c("SL.mean", "SL.glm", "SL.gam", "SL.earth")
  sl_lib_failure <- c("SL.mean", "SL.glm", "SL.gam", "SL.earth")
  range(df$T.tilde)
  df$T.tilde <- df$T.tilde + 1
  k_grid <- 1:max(df$T.tilde)

  message("KM")
  n_sample <- nrow(df)
  km_fit <- survfit(Surv(time = T.tilde, event = Delta) ~ A, data = df)
  surv1_km <- tail(km_fit$surv, km_fit$strata["A=1"])
  time1_km <- tail(km_fit$time, km_fit$strata["A=1"])
  surv0_km <- tail(km_fit$surv, km_fit$strata["A=0"])
  time0_km <- tail(km_fit$time, km_fit$strata["A=0"])
  library(zoo)
  impute_KM <- function(time, km) {
    surv1_km_final <- rep(NA, max(df$T.tilde))
    surv1_km_final[time] <- km
    surv1_km_final <- na.locf(surv1_km_final, na.rm = FALSE)
    surv1_km_final[is.na(surv1_km_final)] <- 1
    surv1_km_final <- c(1, surv1_km_final)
    surv1_km_final <- surv1_km_final[-length(surv1_km_final)]
    return(surv1_km_final)
  }
  surv1_km_final <- impute_KM(time = time1_km, km = surv1_km)
  surv0_km_final <- impute_KM(time = time0_km, km = surv0_km)
  km_fit_1 <- survival_curve$new(t = k_grid, survival = surv1_km_final)
  km_fit_0 <- survival_curve$new(t = k_grid, survival = surv0_km_final)

  message("SL")
  sl_fit <- initial_sl_fit(
    T_tilde = df$T.tilde,
    Delta = df$Delta,
    A = df$A,
    W = data.frame(df[, W_names]),
    t_max = max(df$T.tilde),
    sl_treatment = sl_lib_g,
    sl_censoring = sl_lib_censor,
    sl_failure = sl_lib_failure
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

  sl_density_failure_1_marginal <- sl_fit$density_failure_1$clone(deep = TRUE)
  sl_density_failure_0_marginal <- sl_fit$density_failure_0$clone(deep = TRUE)
  sl_density_failure_1_marginal$survival <- matrix(colMeans(sl_density_failure_1_marginal$survival), nrow = 1)
  sl_density_failure_0_marginal$survival <- matrix(colMeans(sl_density_failure_0_marginal$survival), nrow = 1)
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
    epsilon = 1e-1 / n_sim,
    max_num_interation = 1e2,
    verbose = FALSE
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
    epsilon = 1e-1 / n_sim,
    max_num_interation = 1e2,
    verbose = FALSE
  )
  moss_fit_1 <- survival_curve$new(t = k_grid, survival = psi_moss_1)
  moss_fit_0 <- survival_curve$new(t = k_grid, survival = psi_moss_0)

  # tmle
  message("tmle")
  tmle_fit <- tryCatch({
    tmle_fit <- fit_survtmle(
      T.tilde = df$T.tilde,
      Delta = df$Delta,
      A = df$A,
      W_df = data.frame(df[, W_names]),
      SL.trt = sl_lib_g,
      SL.ctime = sl_lib_censor,
      SL.ftime = sl_lib_failure
    )
  },
  error = function(cond) {
    message("tmle error")
    NULL
  }
  )
  if (is.null(tmle_fit)) {
    tmle_fit_1 <- sl_density_failure_1_marginal$clone(deep = TRUE)
    tmle_fit_0 <- sl_density_failure_0_marginal$clone(deep = TRUE)
    is_tmle1_converge <- FALSE
  } else {
    s_1 <- c(1, tmle_fit$s_1)
    s_1 <- s_1[-length(s_1)]
    s_0 <- c(1, tmle_fit$s_0)
    s_0 <- s_0[-length(s_0)]
    tmle_fit_1 <- survival_curve$new(t = k_grid, survival = s_1)
    tmle_fit_0 <- survival_curve$new(t = k_grid, survival = s_0)
    is_tmle1_converge <- TRUE
  }
  message("moss with l2 submodel")
  moss_hazard_l2 <- MOSS_hazard$new(
    A = df$A,
    T_tilde = df$T.tilde,
    Delta = df$Delta,
    density_failure = sl_fit$density_failure_1,
    density_censor = sl_fit$density_censor_1,
    g1W = sl_fit$g1W,
    A_intervene = 1,
    k_grid = k_grid
  )
  moss_hazard_l1 <- moss_hazard_l2$clone(deep = TRUE)
  psi_moss_l2_1 <- moss_hazard_l2$iterate_onestep(
    method = "l2", epsilon = 1e-1 / sqrt(n_sim), verbose = FALSE
  )
  moss_hazard_l2_1 <- survival_curve$new(t = k_grid, survival = psi_moss_l2_1)

  psi_moss_hazard_l1_1 <- moss_hazard_l1$iterate_onestep(
    method = "l1", epsilon = 1e-1 / sqrt(n_sim), verbose = FALSE
  )
  moss_hazard_l1_1 <- survival_curve$new(t = k_grid, survival = psi_moss_hazard_l1_1)

  # evaluate against truth
  survival_truth_1 <- survival_curve$new(t = k_grid, survival = simulated$true_surv1(k_grid - 1))
  survival_truth_0 <- survival_curve$new(t = k_grid, survival = simulated$true_surv0(k_grid - 1))

  evaluate_moss <- evaluate_metric$new(
    survival = moss_fit_1, survival_truth = survival_truth_1
  )
  df_entropy_moss_1 <- evaluate_moss$evaluate_cross_entropy()
  df_entropy_moss_1$metric_name <- "cross_entropy"
  df_mse_moss_1 <- evaluate_moss$evaluate_mse()
  df_mse_moss_1$metric_name <- "mse"

  evaluate_moss_l2 <- evaluate_metric$new(
    survival = moss_hazard_l2_1, survival_truth = survival_truth_1
  )
  df_entropy_moss_l2_1 <- evaluate_moss_l2$evaluate_cross_entropy()
  df_entropy_moss_l2_1$metric_name <- "cross_entropy"
  df_mse_moss_l2_1 <- evaluate_moss_l2$evaluate_mse()
  df_mse_moss_l2_1$metric_name <- "mse"

  evaluate_moss_l1 <- evaluate_metric$new(
    survival = moss_hazard_l1_1, survival_truth = survival_truth_1
  )
  df_entropy_moss_l1_1 <- evaluate_moss_l1$evaluate_cross_entropy()
  df_entropy_moss_l1_1$metric_name <- "cross_entropy"
  df_mse_moss_l1_1 <- evaluate_moss_l1$evaluate_mse()
  df_mse_moss_l1_1$metric_name <- "mse"

  evaluate_sl <- evaluate_metric$new(
    survival = sl_density_failure_1_marginal, survival_truth = survival_truth_1
  )
  df_entropy_sl_1 <- evaluate_sl$evaluate_cross_entropy()
  df_entropy_sl_1$metric_name <- "cross_entropy"
  df_mse_sl_1 <- evaluate_sl$evaluate_mse()
  df_mse_sl_1$metric_name <- "mse"
  evaluate_ipcw <- evaluate_metric$new(
    survival = ipcw_fit_1, survival_truth = survival_truth_1
  )
  df_entropy_ipcw_1 <- evaluate_ipcw$evaluate_cross_entropy()
  df_entropy_ipcw_1$metric_name <- "cross_entropy"
  df_mse_ipcw_1 <- evaluate_ipcw$evaluate_mse()
  df_mse_ipcw_1$metric_name <- "mse"
  evaluate_ee <- evaluate_metric$new(
    survival = ee_fit_1, survival_truth = survival_truth_1
  )
  df_entropy_ee_1 <- evaluate_ee$evaluate_cross_entropy()
  df_entropy_ee_1$metric_name <- "cross_entropy"
  df_mse_ee_1 <- evaluate_ee$evaluate_mse()
  df_mse_ee_1$metric_name <- "mse"
  evaluate_tmle <- evaluate_metric$new(
    survival = tmle_fit_1, survival_truth = survival_truth_1
  )
  df_entropy_tmle_1 <- evaluate_tmle$evaluate_cross_entropy()
  df_entropy_tmle_1$metric_name <- "cross_entropy"
  df_mse_tmle_1 <- evaluate_tmle$evaluate_mse()
  df_mse_tmle_1$metric_name <- "mse"
  evaluate_km <- evaluate_metric$new(
    survival = km_fit_1, survival_truth = survival_truth_1
  )
  df_entropy_km_1 <- evaluate_km$evaluate_cross_entropy()
  df_entropy_km_1$metric_name <- "cross_entropy"
  df_mse_km_1 <- evaluate_km$evaluate_mse()
  df_mse_km_1$metric_name <- "mse"

  df_mse_moss_1$method <- "MOSS_classic"
  df_mse_moss_l2_1$method <- "MOSS_l2"
  df_mse_moss_l1_1$method <- "MOSS_l1"
  df_mse_sl_1$method <- "super learner"
  df_mse_ipcw_1$method <- "IPCW"
  df_mse_ee_1$method <- "EE"
  df_mse_tmle_1$method <- "TMLE"
  df_mse_km_1$method <- "KM"
  df_entropy_moss_1$method <- "MOSS_classic"
  df_entropy_moss_l2_1$method <- "MOSS_l2"
  df_entropy_moss_l1_1$method <- "MOSS_l1"
  df_entropy_sl_1$method <- "super learner"
  df_entropy_ipcw_1$method <- "IPCW"
  df_entropy_ee_1$method <- "EE"
  df_entropy_tmle_1$method <- "TMLE"
  df_entropy_km_1$method <- "KM"
  df_plot <- plyr::rbind.fill(
    df_mse_moss_1,
    df_mse_moss_l2_1,
    df_mse_moss_l1_1,
    df_mse_sl_1,
    df_mse_ipcw_1,
    df_mse_ee_1,
    df_mse_tmle_1,
    df_mse_km_1,
    df_entropy_moss_1,
    df_entropy_moss_l2_1,
    df_entropy_moss_l1_1,
    df_entropy_sl_1,
    df_entropy_ipcw_1,
    df_entropy_ee_1,
    df_entropy_tmle_1,
    df_entropy_km_1
  )
  # track if the estimators are monotone
  df_plot$is_monotone_tmle1 <- all(diff(as.numeric(tmle_fit_1$survival)) <= 0)
  df_plot$is_monotone_ee1 <- all(diff(as.numeric(ee_fit_1$survival)) <= 0)
  df_plot$is_monotone_tmle0 <- all(diff(as.numeric(tmle_fit_0$survival)) <= 0)
  df_plot$is_monotone_ee0 <- all(diff(as.numeric(ee_fit_0$survival)) <= 0)
  df_plot$is_tmle1_converge <- is_tmle1_converge
  return(df_plot)
}

N_SIMULATION <- 2
# N_SIMULATION <- 1e3
library(foreach)
# library(Rmpi)
# library(doMPI)
# cl <- startMPIcluster()
# registerDoMPI(cl)
# clusterSize(cl) # just to check

library(doSNOW)
library(tcltk)
nw <- parallel:::detectCores() # number of workers
cl <- makeSOCKcluster(nw)
registerDoSNOW(cl)

n_sim_grid <- c(1e2)
# n_sim_grid <- c(1e3, 1e2)
# n_sim_grid <- c(1e3, 5e2, 1e2)
df_metric <- foreach(
  n_sim = n_sim_grid,
  .combine = rbind,
  .packages = c("R6", "MOSS", "survtmle", "survival"),
  .inorder = FALSE,
  .verbose = TRUE
) %:%
  foreach(
    it2 = 1:N_SIMULATION, .combine = rbind, .errorhandling = "remove"
  ) %dopar% {
    df <- do_once(n_sim = n_sim)
    df$id_mcmc <- it2
    df$n <- n_sim
    return(df)
  }
table(df_metric$id_mcmc)

df_monotone <- df_metric %>%
  select(
    n,
    id_mcmc,
    is_monotone_tmle1,
    is_monotone_ee1,
    is_monotone_tmle0,
    is_monotone_ee0,
    is_tmle1_converge
  )
df_monotone <- df_monotone[!duplicated(df_monotone), ]
df_monotone_summary <- df_monotone %>%
  group_by(n) %>%
  summarise(
    cnt = dplyr::n(),
    is_monotone_tmle1 = sum(is_monotone_tmle1 * is_tmle1_converge) / sum(is_tmle1_converge),
    is_monotone_ee1 = mean(is_monotone_ee1),
    is_monotone_tmle0 = sum(is_monotone_tmle0 * is_tmle1_converge) / sum(is_tmle1_converge),
    is_monotone_ee0 = mean(is_monotone_ee0),
    is_tmle1_converge = mean(is_tmle1_converge)
  )

save(df_metric, df_monotone, df_monotone_summary, file = "df_metric.rda")

# shut down for memory
# closeCluster(cl)
# mpi.quit()
stopCluster(cl)
