library(MOSS)
source('./survError.R')
# simulate data
simulate_data <- function(n_sim = 2e2) {
  library(simcausal)
  D <- DAG.empty()
  D <- D +
    node("W1", distr = "rbinom", size = 1, prob = .5) +
    node("W", distr = "runif", min = 0, max = 1.5) +
    node("A", distr = "rbinom", size = 1, prob = .15 + .5*as.numeric(W > .75)) +
    node("Trexp", distr = "rexp", rate = 1 + .5*(W)^2 - .8*A) +
    node("Cweib", distr = "rweibull", shape = 1 - .5*W, scale = 75) +
    # node("T", distr = "rconst", const = ceiling(Trexp*.6)) +
    # node("C", distr = "rconst", const = ceiling(Cweib*.6)) +
    node("T", distr = "rconst", const = ceiling(Trexp*1)) +
    node("C", distr = "rconst", const = ceiling(Cweib*1)) +
    # node("C", distr = "rconst", const = min(10,ceiling(Cweib*1))) +
    # Observed random variable (follow-up time):
    node("T.tilde", distr = "rconst", const = ifelse(T <= C , T, C)) +
    # Observed random variable (censoring indicator, 1 - failure event, 0 - censored):
    node("Delta", distr = "rconst", const = ifelse(T <= C , 1, 0))
  setD <- set.DAG(D)
  dat <- sim(setD, n = n_sim)
  # only grab ID, W's, A, T.tilde, Delta
  Wname <- grep('W', colnames(dat), value = TRUE)
  dat <- dat[,c('ID', Wname, 'A', "T.tilde", "Delta")]

  # input: scalar q, W vector. computes for all W, the S(q|A,W)
  true_surv_one <- function(q, W, A = 1) sapply(W, function(w) {1 - pexp(q, rate = 1 + .5*w^2 - .8*A)})
  # input: vector q. mean(S(q|A,W)|A), average out W. loop over q
  true_surv <- function(q_grid, surv_fn, A) {
    W_grid <- seq(0, 1.5, .01)
    survout <- numeric()
    for (q in q_grid) survout <- c(survout, mean(surv_fn(q = q/1, W = W_grid, A = A)))
    return(survout)
  }
  truth_surv <- function(q) true_surv(q_grid = q, surv_fn = true_surv_one, A = 1)
  truth_surv0 <- function(q) true_surv(q_grid = q, surv_fn = true_surv_one, A = 0)
  return(list(dat = dat, true_surv1 = truth_surv, true_surv0 = truth_surv0))
}

fit_survtmle <- function(dat, Wname = c('W', 'W1')) {
  library(survtmle)
  dat$T.tilde[dat$T.tilde <= 0] <- 1
  t_0 = max(dat$T.tilde)
  fit <- survtmle(ftime = dat$T.tilde,
                  ftype = dat$Delta,
                  trt = dat$A,
                  adjustVars = data.frame(dat[,Wname]),
                  SL.trt = c('SL.glm', 'SL.gam'),
                  SL.ftime = c('SL.glm', 'SL.gam'),
                  SL.ctime = c('SL.glm', 'SL.gam'),
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
  # KM
  library(survival)
  n.data <- nrow(df)
  km.fit <- survfit(Surv(time = T.tilde, event = Delta) ~ A, data = df)
  # SL
  library(MOSS)
  SL_fit <- MOSS::MOSS$new(dat = df, dW = 1, epsilon.step = 1e-2, max.iter = 2e2, verbose = FALSE)
  SL_fit$initial_fit(g.SL.Lib = c("SL.mean","SL.glm",'SL.gam'),
                     Delta.SL.Lib = c("SL.mean","SL.glm",'SL.gam'),
                     ht.SL.Lib = c("SL.mean","SL.glm",'SL.gam'))
  SL_fit$transform_failure_hazard_to_survival()
  # MOSS
  # MOSS_fit <- MOSS::MOSS$new(dat = df, dW = 1, epsilon.step = 1e-1, max.iter = 2e2, verbose = FALSE)
  MOSS_fit <- MOSS::MOSS$new(dat = df, dW = 1, epsilon.step = 1e-2, max.iter = 5e1, verbose = FALSE)
  MOSS_fit$onestep_curve(g.SL.Lib = c("SL.mean","SL.glm",'SL.gam'),
                         Delta.SL.Lib = c("SL.mean","SL.glm",'SL.gam'),
                         ht.SL.Lib = c("SL.mean","SL.glm",'SL.gam'))
  # plot(km.fit, col = 4:5)
  # MOSS_fit$plot_onestep_curve(col = 'green', add = T)
  # curve(true_surv(x), from = 0, 50, add = TRUE)
  # curve(simulated$true_surv0(x), from = 0, 50, add = TRUE)
  # survtmle
  survtmle_out <- fit_survtmle(dat = df)
  s_0 <- survtmle_out$s_0
  s_1 <- survtmle_out$s_1
  t_survtmle <- survtmle_out$time

  # lines(s_1 ~ t_survtmle, col = 'red', lty = 1)
  # compute error
  error_SL <- survError$new(true_surv = true_surv, object = SL_fit, mode = 'SL')
  error_onestep <- survError$new(true_surv = true_surv, object = MOSS_fit, mode = 'onestep')
  error_KM <- survError$new(true_surv = true_surv, object = km.fit)
  error_survtmle <- survError$new(true_surv = true_surv, object = survtmle_out)
  return(list(error_SL = error_SL,
              error_onestep = error_onestep,
              error_KM = error_KM,
              error_survtmle = error_survtmle))
}

# repeat 100 times
N_SIMULATION = 1e2
# N_SIMULATION = 8
library(foreach)
library(Rmpi)
library(doMPI)
cl = startMPIcluster()
registerDoMPI(cl)
clusterSize(cl) # just to check

# library(doSNOW)
# library(tcltk)
# nw <- parallel:::detectCores()  # number of workers
# cl <- makeSOCKcluster(nw)
# registerDoSNOW(cl)

n_sim <- 1e3
all_CI <- foreach(it2 = 1:N_SIMULATION,
                  .combine = c,
                  .packages = c('R6', 'MOSS', 'survtmle', 'survival'),
                  .inorder = FALSE,
                  .errorhandling = 'pass',
                  .verbose = T) %dopar% {
                    if(it2%%10 == 0) print(it2)
                    source('./survError.R')
                    do_once(n_sim = n_sim)
                  }
# shut down for memory
closeCluster(cl)
head(all_CI)

error_SL_list <- all_CI[names(all_CI) == 'error_SL']
error_onestep_list <- all_CI[names(all_CI) == 'error_onestep']
error_KM_list <- all_CI[names(all_CI) == 'error_KM']
error_survtmle_list <- all_CI[names(all_CI) == 'error_survtmle']

error_SL <- survError_list$new(list_of_survError = error_SL_list)$compute()
error_onestep <- survError_list$new(list_of_survError = error_onestep_list)$compute()
error_KM <- survError_list$new(list_of_survError = error_KM_list)$compute()
error_survtmle <- survError_list$new(list_of_survError = error_survtmle_list)$compute()

save(error_SL, error_KM, error_onestep, error_survtmle, n_sim, file = 'scenario.rda')
# ================
# onestep_all_t_fit <- onestep_all_t$survival_df
# sfun_onestep_all_t  <- stepfun(onestep_all_t_fit$T.uniq, c(1, onestep_all_t_fit$s_vec) , f = 1, right = TRUE)
