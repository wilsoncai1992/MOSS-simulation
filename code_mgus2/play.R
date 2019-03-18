t_censor <- 8
# t_censor <- 12
# t_censor <- 16
source("./preprocess.R")
smooth_2d_function <- function(df_prediction, x = "x", y = "y", z = "z", ...) {
  # loess_fit <- loess(
    # df_prediction[, z] ~ df_prediction[, x] + df_prediction[, y], ...
  # )
  # df_smoothed <- df_prediction
  # df_smoothed[, z] <- predict(loess_fit)
  df_smoothed <- df_prediction
  range_y <- range(df_smoothed[, y])
  breaks_y <- seq(range_y[1], range_y[2], length.out = 10)
  df_smoothed[, y] <- cut(
    df_smoothed[, y], breaks = breaks_y, right = FALSE, include.lowest = TRUE, ordered_result = TRUE
  )
  df_smoothed <- df_smoothed %>%
    group_by_(x, y) %>%
    summarise(s = mean(s))
  return(df_smoothed)
}

do_once <- function(df_train, save_plot = FALSE) {
  # T_tilde <- df_train$futime
  # Delta <- df_train$death
  T_tilde <- df_train$ptime
  Delta <- df_train$pstat
  treatment <- df_train[, A_name]
  range(T_tilde)
  k_grid <- 1:max(T_tilde)

  library(MOSS)
  sl_lib_g <- c("SL.mean", "SL.glm", "SL.gam", "SL.earth", "SL.ranger")
  sl_lib_censor <- c("SL.mean", "SL.glm", "SL.gam", "SL.earth", "SL.ranger")
  sl_lib_failure <- c("SL.mean", "SL.glm", "SL.gam", "SL.earth", "SL.ranger")

  sl_fit <- MOSS::initial_sl_fit(
    ftime = T_tilde,
    ftype = Delta,
    trt = treatment,
    adjustVars = data.frame(df_train[, W_name]),
    t_0 = max(T_tilde),
    SL.trt = sl_lib_g,
    SL.ctime = sl_lib_censor,
    SL.ftime = sl_lib_failure
  )
  sl_fit$density_failure_1$hazard_to_survival()
  sl_fit$density_failure_0$hazard_to_survival()
  # WILSON hack no data is t_tilde = 2
  sl_fit$density_failure_1$t <- k_grid
  sl_fit$density_failure_0$t <- k_grid

  sl_density_failure_1_marginal <- sl_fit$density_failure_1$clone(deep = TRUE)
  sl_density_failure_0_marginal <- sl_fit$density_failure_0$clone(deep = TRUE)
  sl_density_failure_1_marginal$survival <- matrix(colMeans(sl_density_failure_1_marginal$survival), nrow = 1)
  sl_density_failure_0_marginal$survival <- matrix(colMeans(sl_density_failure_0_marginal$survival), nrow = 1)

  g_parametric <- c("SL.mean", "SL.glm")
  censor_parametric <- c("SL.mean", "SL.glm")
  failure_parametric <- c("SL.mean", "SL.glm")
  parametric_fit <- MOSS::initial_sl_fit(
    ftime = T_tilde,
    ftype = Delta,
    trt = treatment,
    adjustVars = data.frame(df_train[, W_name]),
    t_0 = max(T_tilde),
    SL.trt = g_parametric,
    SL.ctime = censor_parametric,
    SL.ftime = failure_parametric
  )
  parametric_fit$density_failure_1$hazard_to_survival()
  parametric_fit$density_failure_0$hazard_to_survival()
  parametric_fit$density_failure_1$t <- k_grid
  parametric_fit$density_failure_0$t <- k_grid
  parametric_failure_1_marginal <- parametric_fit$density_failure_1$clone(deep = TRUE)
  parametric_failure_0_marginal <- parametric_fit$density_failure_0$clone(deep = TRUE)
  parametric_failure_1_marginal$survival <- matrix(
    colMeans(parametric_failure_1_marginal$survival), nrow = 1
  )
  parametric_failure_0_marginal$survival <- matrix(
    colMeans(parametric_failure_0_marginal$survival), nrow = 1
  )

  library(ggplot2)
  df_plots <- list()
  for (W_one in W_name) {
    W_plot <- df_train[, W_one]
    df_s1 <- sl_fit$density_failure_1$create_ggplot_df(W = W_plot)
    df_s0 <- sl_fit$density_failure_0$create_ggplot_df(W = W_plot)
    df_c1 <- sl_fit$density_censor_1$create_ggplot_df(W = W_plot)
    df_c0 <- sl_fit$density_censor_0$create_ggplot_df(W = W_plot)
    df_s1 <- smooth_2d_function(df_s1, x = "t", y = "W", z = "s")
    df_s0 <- smooth_2d_function(df_s0, x = "t", y = "W", z = "s")
    df_c1 <- smooth_2d_function(df_c1, x = "t", y = "W", z = "s")
    df_c0 <- smooth_2d_function(df_c0, x = "t", y = "W", z = "s")
    df_s1$component <- "S_{A=1}(t)"
    df_s0$component <- "S_{A=0}(t)"
    df_c1$component <- "G_{A=1}(t)"
    df_c0$component <- "G_{A=0}(t)"
    df_plot <- rbind(df_s1, df_s0, df_c1, df_c0)
    df_plot$covariate <- W_one
    df_plots <- c(df_plots, list(df_plot))
  }
  df_plots <- do.call(rbind, df_plots)
  # reorder the factors of W, after the rbind of data.frames
  lvls <- levels(factor(df_plots$W))
  left_bound <- as.numeric(sapply(strsplit(substr(lvls, 2, 10), split = ","), function(x) x[1]))
  df_plots$W <- factor(df_plots$W, levels = lvls[order(left_bound)], ordered = TRUE)
  # gg_panel <- ggplot(df_plots, aes(x = t, y = round(W, digits = 1), z = s)) +
  #   geom_raster(aes(fill = s), interpolate = TRUE) +
  #   xlim(c(1, max(df_plots$t))) +
  #   ylab("W") +
  #   theme_bw() +
  #   facet_grid(covariate ~ component, scales = 'free_y')
  gg_panel <- ggplot(df_plots, aes(x = t, y = W, z = s)) +
    geom_tile(aes(fill = s)) +
    xlim(c(1, max(df_plots$t))) +
    ylab("W") +
    theme_bw() +
    facet_grid(covariate ~ component, scales = 'free_y')
  if (save_plot) ggsave(gg_panel, filename = 'panel1.png', width = 8, height = 6)

  denom_1 <- sl_fit$density_censor_1$survival * sl_fit$g1W
  denom_0 <- sl_fit$density_censor_0$survival * (1 - sl_fit$g1W)
  library(stargazer)

  tbl_score <- data.frame(1 / denom_1)
  colnames(tbl_score) <- 1:ncol(tbl_score)
  stargazer(tbl_score)
  tbl_score <- data.frame(1 / denom_0)
  colnames(tbl_score) <- 1:ncol(tbl_score)
  stargazer(tbl_score)
  # summary(1 / denom_1)
  # summary(1 / denom_0)

  density_censor_1_truncated <- sl_fit$density_censor_1$clone(deep = TRUE)
  density_censor_0_truncated <- sl_fit$density_censor_0$clone(deep = TRUE)
  matrix_g1W <- do.call(
    cbind, rep(list(sl_fit$g1W), ncol(sl_fit$density_censor_1$survival))
  )
  # threshold <- 0.02
  # threshold <- 0.01
  threshold <- 1e-5
  density_censor_1_truncated$survival[denom_1 <= threshold] <- (threshold / matrix_g1W)[denom_1 <= threshold]
  density_censor_0_truncated$survival[denom_0 <= threshold] <- (threshold / (1 - matrix_g1W))[denom_0 <= threshold]
  denom_1_truncated <- density_censor_1_truncated$survival * sl_fit$g1W
  denom_0_truncated <- density_censor_0_truncated$survival * (1 - sl_fit$g1W)
  summary(1 / denom_1_truncated)
  summary(1 / denom_0_truncated)

  message("ipcw + ee")
  ipcw_fit_1_all <- repeat_t_grid$new(
    method = ipcw,
    A = treatment,
    T_tilde = T_tilde,
    Delta = Delta,
    density_failure = sl_fit$density_failure_1,
    # density_censor = sl_fit$density_censor_1,
    density_censor = density_censor_1_truncated,
    g1W = sl_fit$g1W,
    A_intervene = 1
  )$fit(k_grid = k_grid)
  ipcw_fit_0_all <- repeat_t_grid$new(
    method = ipcw,
    A = treatment,
    T_tilde = T_tilde,
    Delta = Delta,
    density_failure = sl_fit$density_failure_0,
    # density_censor = sl_fit$density_censor_0,
    density_censor = density_censor_0_truncated,
    g1W = sl_fit$g1W,
    A_intervene = 0
  )$fit(k_grid = k_grid)
  ee_fit_1_all <- repeat_t_grid$new(
    method = ee,
    A = treatment,
    T_tilde = T_tilde,
    Delta = Delta,
    density_failure = sl_fit$density_failure_1,
    # density_censor = sl_fit$density_censor_1,
    density_censor = density_censor_1_truncated,
    g1W = sl_fit$g1W,
    A_intervene = 1
  )$fit(k_grid = k_grid)
  ee_fit_0_all <- repeat_t_grid$new(
    method = ee,
    A = treatment,
    T_tilde = T_tilde,
    Delta = Delta,
    density_failure = sl_fit$density_failure_0,
    # density_censor = sl_fit$density_censor_0,
    density_censor = density_censor_0_truncated,
    g1W = sl_fit$g1W,
    A_intervene = 0
  )$fit(k_grid = k_grid)
  ipcw_fit_1 <- survival_curve$new(t = k_grid, survival = ipcw_fit_1_all)
  ipcw_fit_0 <- survival_curve$new(t = k_grid, survival = ipcw_fit_0_all)
  ee_fit_1 <- survival_curve$new(t = k_grid, survival = ee_fit_1_all)
  ee_fit_0 <- survival_curve$new(t = k_grid, survival = ee_fit_0_all)


  source('../fit_survtmle.R')
  message("tmle")
  tmle_fit <- tryCatch({
      tmle_fit <- fit_survtmle(
        T.tilde = T_tilde,
        Delta = Delta,
        A = treatment,
        W_df = data.frame(df_train[, W_name]),
        SL.trt = g_parametric,
        SL.ctime = censor_parametric,
        SL.ftime = failure_parametric
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
    s_1 <- c(1, tmle_fit$s_1)
    s_1 <- s_1[-length(s_1)]
    s_0 <- c(1, tmle_fit$s_0)
    s_0 <- s_0[-length(s_0)]
    tmle_fit_1 <- survival_curve$new(t = k_grid, survival = s_1)
    tmle_fit_0 <- survival_curve$new(t = k_grid, survival = s_0)
  }

  message("moss with hazard submodel")
  moss_hazard_fit <- MOSS_hazard$new(
    A = treatment,
    T_tilde = T_tilde,
    Delta = Delta,
    density_failure = sl_fit$density_failure_1,
    # density_censor = sl_fit$density_censor_1,
    density_censor = density_censor_1_truncated,
    g1W = sl_fit$g1W,
    A_intervene = 1,
    k_grid = k_grid
  )
  psi_moss_hazard_1 <- moss_hazard_fit$iterate_onestep(
    epsilon = 1e-3, max_num_interation = 3e1, verbose = TRUE
  )
  moss_hazard_fit_1 <- survival_curve$new(t = k_grid, survival = psi_moss_hazard_1)

  moss_hazard_fit_0 <- MOSS_hazard$new(
    A = treatment,
    T_tilde = T_tilde,
    Delta = Delta,
    density_failure = sl_fit$density_failure_0,
    # density_censor = sl_fit$density_censor_0,
    density_censor = density_censor_0_truncated,
    g1W = sl_fit$g1W,
    A_intervene = 0,
    k_grid = k_grid
  )
  psi_moss_hazard_0 <- moss_hazard_fit_0$iterate_onestep(
    epsilon = 1e-3, max_num_interation = 3e1, verbose = TRUE
  )
  moss_hazard_fit_0 <- survival_curve$new(t = k_grid, survival = psi_moss_hazard_0)

  message("moss for ATE")
  moss_hazard_ate_fit <- MOSS_hazard_ate$new(
    A = treatment,
    T_tilde = T_tilde,
    Delta = Delta,
    density_failure = sl_fit$density_failure_1,
    # density_censor = sl_fit$density_censor_1,
    density_censor = density_censor_1_truncated,
    density_failure_0 = sl_fit$density_failure_0,
    # density_censor_0 = sl_fit$density_censor_0,
    density_censor_0 = density_censor_0_truncated,
    g1W = sl_fit$g1W,
    k_grid = k_grid
  )
  psi_moss_ate <- moss_hazard_ate_fit$iterate_onestep(
    epsilon = 5e-2, max_num_interation = 2e1, verbose = TRUE
  )
  moss_ate_fit <- survival_curve$new(t = k_grid, survival = psi_moss_ate)

  # ============================================================================
  # plot the A=1 survival curves

  is_monotone_tmle1 <- all(diff(as.numeric(tmle_fit_1$survival)) <= 0)
  is_monotone_ipcw1 <- all(diff(as.numeric(ipcw_fit_1$survival)) <= 0)
  is_monotone_ee1 <- all(diff(as.numeric(ee_fit_1$survival)) <= 0)
  df_curve_sl1 <- sl_density_failure_1_marginal$create_ggplot_df()
  df_curve_tmle1 <- tmle_fit_1$create_ggplot_df()
  df_curve_moss1 <- moss_hazard_fit_1$create_ggplot_df()
  df_curve_ipcw1 <- ipcw_fit_1$create_ggplot_df()
  df_curve_ee1 <- ee_fit_1$create_ggplot_df()
  df_curve_glm1 <- parametric_failure_1_marginal$create_ggplot_df()
  df_curve_sl1$method <- "super learner"
  df_curve_tmle1$method <- "TMLE"
  df_curve_moss1$method <- "one-step TMLE"
  df_curve_ipcw1$method <- "IPCW"
  df_curve_ee1$method <- "EE"
  df_curve_glm1$method <- "linear model"
  df_curve1 <- rbind(
    df_curve_sl1,
    df_curve_tmle1,
    df_curve_moss1,
    df_curve_ipcw1,
    df_curve_ee1,
    df_curve_glm1
  )
  df_curve1$counterfactual <- "S_{A=1}(t)"

  is_monotone_tmle0 <- all(diff(as.numeric(tmle_fit_0$survival)) <= 0)
  is_monotone_ipcw0 <- all(diff(as.numeric(ipcw_fit_0$survival)) <= 0)
  is_monotone_ee0 <- all(diff(as.numeric(ee_fit_0$survival)) <= 0)
  df_curve_sl0 <- sl_density_failure_0_marginal$create_ggplot_df()
  df_curve_tmle0 <- tmle_fit_0$create_ggplot_df()
  df_curve_moss0 <- moss_hazard_fit_0$create_ggplot_df()
  df_curve_ipcw0 <- ipcw_fit_0$create_ggplot_df()
  df_curve_ee0 <- ee_fit_0$create_ggplot_df()
  df_curve_glm0 <- parametric_failure_0_marginal$create_ggplot_df()
  df_curve_sl0$method <- "super learner"
  df_curve_tmle0$method <- "TMLE"
  df_curve_moss0$method <- "one-step TMLE"
  df_curve_ipcw0$method <- "IPCW"
  df_curve_ee0$method <- "EE"
  df_curve_glm0$method <- "linear model"
  df_curve0 <- rbind(
    df_curve_sl0,
    df_curve_tmle0,
    df_curve_moss0,
    df_curve_ipcw0,
    df_curve_ee0,
    df_curve_glm0
  )
  df_curve0$counterfactual <- "S_{A=0}(t)"

  gg_s <- ggplot(
      data = rbind(df_curve1, df_curve0),
      aes(x = t, y = s, color = method)
    ) +
    geom_line() +
    theme_bw() +
    facet_wrap(. ~ counterfactual, ncol = 2) +
    theme(legend.position = "bottom")

  ate_moss <- moss_ate_fit$create_ggplot_df()
  ate_moss$method <- "one-step TMLE"
  ate_ee <- df_curve_ee1
  ate_ee$s <- ate_ee$s - df_curve_ee0$s
  ate_tmle <- df_curve_tmle1
  ate_tmle$s <- ate_tmle$s - df_curve_tmle0$s
  ate_ipcw <- df_curve_ipcw1
  ate_ipcw$s <- NA
  ate_sl <- df_curve_sl1
  ate_sl$s <- ate_sl$s - df_curve_sl0$s
  ate_glm <- df_curve_glm1
  ate_glm$s <- ate_glm$s - df_curve_glm0$s

  # compute ci
  eic_fit_1 <- eic$new(
    A = treatment,
    T_tilde = T_tilde,
    Delta = Delta,
    density_failure = moss_hazard_ate_fit$density_failure,
    density_censor = moss_hazard_ate_fit$density_censor,
    g1W = moss_hazard_ate_fit$g1W,
    psi = psi_moss_ate,
    A_intervene = 1
  )$all_t(k_grid = k_grid)
  eic_fit_0 <- eic$new(
    A = treatment,
    T_tilde = T_tilde,
    Delta = Delta,
    density_failure = moss_hazard_ate_fit$density_failure_0,
    density_censor = moss_hazard_ate_fit$density_censor_0,
    g1W = moss_hazard_ate_fit$g1W,
    psi = psi_moss_ate,
    A_intervene = 0
  )$all_t(k_grid = k_grid)
  eic_moss <- eic_fit_1 - eic_fit_0
  ci_moss <- MOSS::compute_simultaneous_ci(eic_fit = eic_moss)
  df_ci_moss <- data.frame(
    lower = psi_moss_ate - ci_moss, upper = psi_moss_ate + ci_moss
  )
  ate_moss <- cbind(ate_moss, df_ci_moss)

  eic_fit_1 <- eic$new(
    A = treatment,
    T_tilde = T_tilde,
    Delta = Delta,
    density_failure = sl_fit$density_failure_1,
    density_censor = sl_fit$density_censor_1,
    g1W = sl_fit$g1W,
    psi = ate_ee$s,
    A_intervene = 1
  )$all_t(k_grid = k_grid)
  eic_fit_0 <- eic$new(
    A = treatment,
    T_tilde = T_tilde,
    Delta = Delta,
    density_failure = sl_fit$density_failure_0,
    density_censor = sl_fit$density_censor_0,
    g1W = sl_fit$g1W,
    psi = ate_ee$s,
    A_intervene = 0
  )$all_t(k_grid = k_grid)
  eic_ee <- eic_fit_1 - eic_fit_0
  ci_ee <- MOSS::compute_simultaneous_ci(eic_fit = eic_ee)
  df_ci_ee <- data.frame(
    lower = ate_ee$s - ci_ee, upper = ate_ee$s + ci_ee
  )
  ate_ee <- cbind(ate_ee, df_ci_ee)

  ate_tmle$upper <- ate_tmle$lower <- NA
  ate_ipcw$upper <- ate_ipcw$lower <- NA
  ate_sl$upper <- ate_sl$lower <- NA
  ate_glm$upper <- ate_glm$lower <- NA
  df_ate <- rbind(ate_moss, ate_ee, ate_tmle, ate_ipcw, ate_sl, ate_glm)

  gg_ate <- ggplot(
      data = df_ate,
      aes(x = t, y = s, color = method)
    ) +
    geom_line() +
    geom_line(aes(x = t, y = lower, color = method), lty = 2) +
    geom_line(aes(x = t, y = upper, color = method), lty = 2) +
    geom_hline(yintercept = 0, lty = 4) +
    theme_bw() +
    theme(legend.position = "bottom")

  library(ggpubr)
  panel2 <- ggarrange(
    gg_s,
    gg_ate,
    nrow = 1,
    ncol = 2,
    labels = "AUTO",
    widths = c(2, 1),
    common.legend = TRUE,
    legend = "bottom"
  )
  if (save_plot) ggsave(panel2, filename = 'panel2.png', width = 12, height = 5)
  df_stats <- data.frame(
    is_monotone_tmle1,
    is_monotone_ipcw1,
    is_monotone_ee1,
    is_monotone_tmle0,
    is_monotone_ipcw0,
    is_monotone_ee0
  )
  return(df_stats)
}
do_once(df_train, save_plot = TRUE)


# N_REPEAT <- 1e2
# # N_REPEAT <- 2
# # n_subsample <- nrow(df_train)
# # n_subsample_grid <- c(1e2, 5e2)
# n_subsample_grid <- c(1e2, 5e2, 1e3)

# library(foreach)
# library(Rmpi)
# library(doMPI)
# cl = startMPIcluster()
# registerDoMPI(cl)
# clusterSize(cl) # just to check

# # library(doSNOW)
# # library(tcltk)
# # nw <- parallel:::detectCores() # number of workers
# # cl <- makeSOCKcluster(nw)
# # registerDoSNOW(cl)

# df_stats_list <- foreach(
#   n_subsample = n_subsample_grid,
#   .combine = rbind,
#   .packages = c("R6", "MOSS", "survtmle", "dplyr"),
#   .inorder = FALSE,
#   .errorhandling = "remove",
#   .verbose = TRUERUE
# ) %:%
#   foreach(
#     i = 1:N_REPEAT, .combine = rbind, .errorhandling = "remove"
#   ) %dopar% {
#     df_stats <- do_once(
#       # sample_n(df_train, size = n_subsample, replace = TRUE)
#       sample_n(df_train, size = n_subsample, replace = F)
#     )
#     df_stats$t_censor <- t_censor
#     df_stats$n_subsample <- n_subsample
#     return(df_stats)
# }
# df_summary <- df_stats_list %>%
#   group_by(t_censor, n_subsample) %>%
#   mutate(cnt = dplyr::n()) %>%
#   summarise_all(mean)
# save(df_stats_list, df_summary, file = "df_stats_list.rda")

