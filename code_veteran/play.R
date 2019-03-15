library(survival)
library(dplyr)
library(tidyr)
data("veteran")
df <- veteran
# df <- df %>% filter(time < 800)
# outcome
summary(df$time)
table(df$status)

# baseline
summary(df$celltype)
summary(df$karno)
summary(df$diagtime)
summary(df$age)
summary(df$prior)

# treatment
table(df$trt)

Y_name <- c("time", 'status')
A_name <- c("trt")
# W_name <- c("celltype", "karno", "diagtime", "age", "prior")
W_name <- c("karno", "diagtime", "age", "prior")

df_train <- df %>%
  select(c(Y_name, A_name, W_name)) %>%
  mutate(
    time = floor(time / 50) + 1,
    trt = trt - 1
  ) %>%
  mutate(prior = prior == 10) %>%
  drop_na()
range(df_train$time)
k_grid <- 1:max(df_train$time)

library(MOSS)
sl_lib_g <- c("SL.mean", "SL.glm", "SL.gam")
sl_lib_censor <- c("SL.mean", "SL.glm", "SL.gam", "SL.earth")
sl_lib_failure <- c("SL.mean", "SL.glm", "SL.gam", "SL.earth")

sl_fit <- MOSS::initial_sl_fit(
  ftime = df_train$time,
  ftype = df_train$status,
  trt = df_train$trt,
  adjustVars = data.frame(df_train[, W_name]),
  t_0 = max(df_train$time),
  SL.trt = sl_lib_g,
  SL.ctime = sl_lib_censor,
  SL.ftime = sl_lib_failure
)
sl_fit$density_failure_1$hazard_to_survival()
sl_fit$density_failure_0$hazard_to_survival()
# WILSON hack no data is t_tilde = 2
sl_fit$density_failure_1$t <- k_grid
sl_fit$density_failure_0$t <- k_grid

s1 <- sl_fit$density_failure_1$display(type = 'survival', W = df_train$karno)
s0 <- sl_fit$density_failure_0$display(type = 'survival', W = df_train$karno)
c1 <- sl_fit$density_censor_1$display(type = 'survival', W = df_train$karno)
c0 <- sl_fit$density_censor_0$display(type = 'survival', W = df_train$karno)
library(ggpubr)
ggarrange(s1, s0, c1, c0, ncol = 2, nrow = 2, common.legend = TRUE, legend = 'bottom', labels = "AUTO")
summary(sl_fit$g1W)

message("ipcw + ee")
ipcw_fit_1_all <- repeat_t_grid$new(
  method = ipcw,
  A = df_train$trt,
  T_tilde = df_train$time,
  Delta = df_train$status,
  density_failure = sl_fit$density_failure_1,
  density_censor = sl_fit$density_censor_1,
  g1W = sl_fit$g1W,
  A_intervene = 1
)$fit(k_grid = k_grid)
ipcw_fit_0_all <- repeat_t_grid$new(
  method = ipcw,
  A = df_train$trt,
  T_tilde = df_train$time,
  Delta = df_train$status,
  density_failure = sl_fit$density_failure_0,
  density_censor = sl_fit$density_censor_0,
  g1W = sl_fit$g1W,
  A_intervene = 0
)$fit(k_grid = k_grid)
ee_fit_1_all <- repeat_t_grid$new(
  method = ee,
  A = df_train$trt,
  T_tilde = df_train$time,
  Delta = df_train$status,
  density_failure = sl_fit$density_failure_1,
  density_censor = sl_fit$density_censor_1,
  g1W = sl_fit$g1W,
  A_intervene = 1
)$fit(k_grid = k_grid)
ee_fit_0_all <- repeat_t_grid$new(
  method = ee,
  A = df_train$trt,
  T_tilde = df_train$time,
  Delta = df_train$status,
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


source('./fit_survtmle.R')
message("tmle")
tmle_fit <- tryCatch({
    tmle_fit <- fit_survtmle(
      T.tilde = df_train$time,
      Delta = df_train$status,
      A = df_train$trt,
      W_df = data.frame(df_train[, W_name])
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
  A = df_train$trt,
  T_tilde = df_train$time,
  Delta = df_train$status,
  density_failure = sl_fit$density_failure_1,
  density_censor = sl_fit$density_censor_1,
  g1W = sl_fit$g1W,
  A_intervene = 1,
  k_grid = k_grid
)
psi_moss_hazard_1 <- moss_hazard_fit$iterate_onestep(epsilon = 1e-2, max_num_interation = 5e1, verbose = T)
moss_hazard_fit_1 <- survival_curve$new(t = k_grid, survival = psi_moss_hazard_1)
moss_hazard_fit_1$display(type = 'survival')

message("moss for ATE")
moss_hazard_ate_fit <- MOSS_hazard_ate$new(
  A = df_train$trt,
  T_tilde = df_train$time,
  Delta = df_train$status,
  density_failure = sl_fit$density_failure_1,
  density_censor = sl_fit$density_censor_1,
  density_failure_0 = sl_fit$density_failure_0,
  density_censor_0 = sl_fit$density_censor_0,
  g1W = sl_fit$g1W,
  k_grid = k_grid
)
psi_moss_ate <- moss_hazard_ate_fit$iterate_onestep(epsilon = 1e-2, max_num_interation = 5e1, verbose = T)
moss_ate_fit <- survival_curve$new(t = k_grid, survival = psi_moss_ate)
moss_ate_fit$display(type = 'survival')

is_monotone_tmle <- all(diff(as.numeric(tmle_fit_1$survival)) <= 0)
is_monotone_ipcw <- all(diff(as.numeric(ipcw_fit_1$survival)) <= 0)
is_monotone_ee <- all(diff(as.numeric(ee_fit_1$survival)) <= 0)
df_curve_sl1 <- sl_density_failure_1_marginal$create_ggplot_df()
df_curve_tmle1 <- tmle_fit_1$create_ggplot_df()
df_curve_moss1 <- moss_hazard_fit_1$create_ggplot_df()
df_curve_ipcw1 <- ipcw_fit_1$create_ggplot_df()
df_curve_ee1 <- ee_fit_1$create_ggplot_df()
df_curve_sl1$method <- "super learner"
df_curve_tmle1$method <- "TMLE"
df_curve_moss1$method <- "MOSS"
df_curve_ipcw1$method <- "IPCW"
df_curve_ee1$method <- "EE"
df_curve <- rbind(
  df_curve_sl1,
  df_curve_tmle1,
  df_curve_moss1,
  df_curve_ipcw1,
  df_curve_ee1
)
library(ggplot2)
ggplot(data = df_curve, aes(x = t, y = s, color = method)) + geom_line()
