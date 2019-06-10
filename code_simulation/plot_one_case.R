library(survival)
library(MOSS)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(survtmle)
source("../fit_survtmle.R")
# simulate data
source("./simulate_data.R")

do_once <- function(n_sim = 2e2) {
  simulated <- simulate_data(n_sim = n_sim)
  df <- simulated$dat
  true_surv <- simulated$true_surv1

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
  tmle_fit <- tryCatch({
    tmle_fit <- fit_survtmle(
      T.tilde = df$T.tilde,
      Delta = df$Delta,
      A = df$A,
      W_df = data.frame(df[, c("W", "W1")]),
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
  } else {
    tmle_fit_1 <- survival_curve$new(t = k_grid, survival = tmle_fit$s_1)
    tmle_fit_0 <- survival_curve$new(t = k_grid, survival = tmle_fit$s_0)
  }
  survival_truth_1 <- survival_curve$new(t = k_grid, survival = simulated$true_surv1(k_grid))
  survival_truth_0 <- survival_curve$new(t = k_grid, survival = simulated$true_surv0(k_grid))

  is_monotone_tmle <- all(diff(as.numeric(tmle_fit_1$survival)) <= 0)
  is_monotone_ipcw <- all(diff(as.numeric(ipcw_fit_1$survival)) <= 0)
  is_monotone_ee <- all(diff(as.numeric(ee_fit_1$survival)) <= 0)
  df_curve_sl1 <- sl_density_failure_1_marginal$create_ggplot_df()
  df_curve_tmle1 <- tmle_fit_1$create_ggplot_df()
  df_curve_moss1 <- moss_fit_1$create_ggplot_df()
  df_curve_km1 <- km_fit_1$create_ggplot_df()
  df_curve_ipcw1 <- ipcw_fit_1$create_ggplot_df()
  df_curve_ee1 <- ee_fit_1$create_ggplot_df()
  df_curve_sl1$method <- "super learner"
  df_curve_tmle1$method <- "TMLE"
  df_curve_moss1$method <- "MOSS"
  df_curve_km1$method <- "KM"
  df_curve_ipcw1$method <- "IPCW"
  df_curve_ee1$method <- "EE"
  df_curve <- rbind(
    df_curve_sl1,
    df_curve_tmle1,
    df_curve_moss1,
    df_curve_km1,
    df_curve_ipcw1,
    df_curve_ee1
  )
  # if (!is_monotone_tmle & !is_monotone_ipcw & !is_monotone_ee) {
  if (!is_monotone_tmle & !is_monotone_ee) {
    # if (!is_monotone_tmle) {
    return(df_curve)
  } else {
    return(NULL)
  }
}

N_SIMULATION <- 2e1
library(foreach)

library(doSNOW)
library(tcltk)
nw <- parallel:::detectCores()
cl <- makeSOCKcluster(nw)
registerDoSNOW(cl)

n_sim_grid <- c(1e2)
# n_sim_grid <- c(1e3)
df_metric <- foreach(
  n_sim = n_sim_grid,
  .combine = rbind,
  .packages = c("R6", "MOSS", "survtmle", "survival"),
  .inorder = FALSE,
  .errorhandling = "remove",
  .verbose = TRUE
) %:%
  foreach(it2 = 1:N_SIMULATION, .combine = rbind, .errorhandling = "remove") %dopar% {
    df <- do_once(n_sim = n_sim)
    if (!is.null(df)) {
      df$id_mcmc <- it2
      return(df)
    } else {
      return(NULL)
    }
  }
unique(df_metric$id_mcmc)
gglist <- list()
for (idx in unique(df_metric$id_mcmc)) {
  gg <- ggplot(df_metric %>% filter(id_mcmc == idx), aes(t, s, color = method)) +
    geom_line() +
    theme_bw() +
    theme(legend.position = "bottom")
  gglist <- c(gglist, list(gg))
}
gg_panel <- ggarrange(plotlist = gglist, legend = "bottom")
ggpubr::ggexport(plotlist = gglist, filename = "panel.pdf", width = 4, height = 4)

# shut down for memory
# closeCluster(cl)
# mpi.quit()
stopCluster(cl)
