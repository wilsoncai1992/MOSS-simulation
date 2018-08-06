library(MOSS)
source('./survError.R')
# simulate data
# source('./simulate_data_0.R')
# source('./simulate_data_1.R')
# source('./simulate_data_2.R')
source('./simulate_data_3.R')
library(survtmle)
fit_survtmle <- function(dat, Wname = c('W', 'W1')) {
  dat$T.tilde[dat$T.tilde <= 0] <- 1
  t_0 <- min(max(dat$T.tilde), 20)
  fit <- survtmle(ftime = dat$T.tilde,
                  ftype = dat$Delta,
                  trt = dat$A,
                  adjustVars = data.frame(dat[,Wname]),
                  SL.trt = c("SL.mean", 'SL.glm', 'SL.gam'),
                  SL.ftime = c("SL.mean", 'SL.glm', 'SL.gam'),
                  SL.ctime = c("SL.mean", 'SL.glm', 'SL.gam'),
                  method = "hazard",
                  # tol = 1/(sqrt(length(dat$T.tilde)))/10,
                  returnIC = TRUE,
                  verbose = FALSE
  )
  # browser()
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
  # ts.plot(s_1)
  # ts.plot(s_0)
  return(data.frame(time = 1:t_0, s_0 = s_0, s_1 = s_1))
}
library(survival)
library(MOSS)
do_once <- function(n_sim = 2e2) {
  simulated <- simulate_data(n_sim = n_sim)
  df <- simulated$dat
  true_surv <- simulated$true_surv1
  # KM
  n.data <- nrow(df)
  km.fit <- survfit(Surv(time = T.tilde, event = Delta) ~ A, data = df)
  # SL
  SL_fit <- MOSS::MOSS$new(dat = df, dW = 1, epsilon.step = 1e-1, max.iter = 2e2, verbose = FALSE)
  SL_fit$initial_fit(g.SL.Lib = c("SL.mean","SL.glm",'SL.gam'),
                     Delta.SL.Lib = c("SL.mean","SL.glm",'SL.gam'),
                     ht.SL.Lib = c("SL.mean","SL.glm",'SL.gam'))
  SL_fit$transform_failure_hazard_to_survival()
  # MOSS
  MOSS_fit <- MOSS::MOSS$new(dat = df, dW = 1, epsilon.step = 1e-1, max.iter = 4e2, verbose = FALSE)
  # MOSS_fit <- MOSS::MOSS$new(dat = df, dW = 1, epsilon.step = 1e0, max.iter = 1e2, verbose = FALSE, tol = 1e-1/nrow(df))
  # MOSS_fit <- MOSS::MOSS$new(dat = df, dW = 1, epsilon.step = 1e1, max.iter = 5e1, verbose = FALSE)
  MOSS_fit$onestep_curve(g.SL.Lib = c("SL.mean","SL.glm",'SL.gam'),
                         Delta.SL.Lib = c("SL.mean","SL.glm",'SL.gam'),
                         ht.SL.Lib = c("SL.mean","SL.glm",'SL.gam'))
  # plot(km.fit, col = 4:5)
  # MOSS_fit$plot_onestep_curve(col = 'green', add = T)
  # curve(true_surv(x), from = 0, 50, add = TRUE)
  # curve(simulated$true_surv0(x), from = 0, 50, add = TRUE)

  # compute error
  error_SL <- survError$new(true_surv = true_surv, object = SL_fit, mode = 'SL')
  error_onestep <- survError$new(true_surv = true_surv, object = MOSS_fit, mode = 'onestep')
  error_KM <- survError$new(true_surv = true_surv, object = km.fit)

  # survtmle
  # browser()
  error_survtmle <- tryCatch({
    # survtmle_out <- fit_survtmle(dat = df,Wname = c('W', 'W1', 'W2', 'W3', 'W4', 'W5'))
    survtmle_out <- fit_survtmle(dat = df,Wname = c('W', 'W1'))
    s_0 <- survtmle_out$s_0
    s_1 <- survtmle_out$s_1
    t_survtmle <- survtmle_out$time
    out <- survError$new(true_surv = true_surv, object = survtmle_out)
    # browser()
    # lines(s_1 ~ t_survtmle, col = 'red', lty = 1) #WILSON
    # return(out)
    out
  },error = function(error_condition) {
    out <- error_SL
    # return(out)
    out
  })
  return(list(error_SL = error_SL,
              error_onestep = error_onestep,
              error_KM = error_KM,
              error_survtmle = error_survtmle))
}

# repeat 100 times
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

n_sim <- 1e2
# n_sim <- 5e2
# n_sim <- 1e3
all_CI <- foreach(it2 = 1:N_SIMULATION,
                  .combine = c,
                  .packages = c('R6', 'MOSS', 'survtmle', 'survival'),
                  .inorder = FALSE,
                  .errorhandling = 'pass',
                  .verbose = T) %dopar% {
                  # .verbose = T) %do% {
                    # if(it2%%10 == 0) print(it2)
                    source('./survError.R')
                    do_once(n_sim = n_sim)
                  }
# shut down for memory
# closeCluster(cl)
head(all_CI)
table(names(all_CI))

error_SL_list <- all_CI[names(all_CI) == 'error_SL']
error_onestep_list <- all_CI[names(all_CI) == 'error_onestep']
error_KM_list <- all_CI[names(all_CI) == 'error_KM']
error_survtmle_list <- all_CI[names(all_CI) == 'error_survtmle']

length(error_SL_list)
length(error_onestep_list)
length(error_KM_list)
length(error_survtmle_list)

error_SL <- survError_list$new(list_of_survError = error_SL_list)$compute()
error_onestep <- survError_list$new(list_of_survError = error_onestep_list)$compute()
error_KM <- survError_list$new(list_of_survError = error_KM_list)$compute()
error_survtmle <- survError_list$new(list_of_survError = error_survtmle_list)$compute()

save(error_SL, error_KM, error_onestep, error_survtmle, n_sim, file = 'scenario.rda')
# ================
# onestep_all_t_fit <- onestep_all_t$survival_df
# sfun_onestep_all_t  <- stepfun(onestep_all_t_fit$T.uniq, c(1, onestep_all_t_fit$s_vec) , f = 1, right = TRUE)
