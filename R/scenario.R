library(MOSS)
# simulate data
simulate_data <- function(n_sim = 2e2) {
  library(simcausal)
  D <- DAG.empty()
  
  D <- D +
    node("W1", distr = "rbinom", size = 1, prob = .5) +
    node("W", distr = "runif", min = 0, max = 1.5) +
    node("A", distr = "rbinom", size = 1, prob = .15 + .5*as.numeric(W > .75)) +
    node("Trexp", distr = "rexp", rate = 1 + .7*(W)^2 - .8*A) +
    node("Cweib", distr = "rweibull", shape = 1 - .5*W, scale = 75) +
    node("T", distr = "rconst", const = round(Trexp*20,0)) +
    node("C", distr = "rconst", const = round(Cweib*20, 0)) +
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
  true_surv_one <- function(q, W, A = 1) sapply(W, function(w) {1 - pexp(q, rate = 1 + .7*w^2 - .8*A)})
  # input: vector q. mean(S(q|A,W)|A), average out W. loop over q
  true_surv <- function(q_grid, surv_fn, A) {
    W_grid <- seq(0, 1.5, .01)
    survout <- numeric()
    for (q in q_grid) survout <- c(survout, mean(surv_fn(q = q/20, W = W_grid, A = A)))
    return(survout)
  }
  truth_surv <- function(q) true_surv(q_grid = q, surv_fn = true_surv_one, A = 1)
  # truesurvExp <- true_surv(q_grid = q_grid, foo = true_surv_one, A = 1)
  
  return(list(dat = dat, true_surv = truth_surv))
}

simulated <- simulate_data(n_sim = 2e2)
df <- simulated$dat
true_surv <- simulated$true_surv

# KM
library(survival)
n.data <- nrow(df)
km.fit <- survfit(Surv(time = T.tilde, event = Delta) ~ A, data = df)
# MOSS
library(MOSS)
MOSS_fit <- MOSS::MOSS$new(dat = df, dW = 1, epsilon.step = 1e-2, max.iter = 2e2, verbose = FALSE)
MOSS_fit$onestep_curve()
MOSS_fit$Psi.hat
# survtmle
fit_survtmle <- function(dat, Wname = c('W', 'W1')) {
  library(survtmle)
  dat$T.tilde[dat$T.tilde <= 0] <- 1
  t_0 = max(dat$T.tilde)
  fit <- survtmle(ftime = dat$T.tilde, 
                  ftype = dat$Delta,
                  trt = dat$A, 
                  adjustVars = data.frame(dat[,Wname]),
                  # t0 = max(dat$T.tilde),
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
survtmle_out <- fit_survtmle(dat = df)
s_0 <- survtmle_out$s_0
s_1 <- survtmle_out$s_1
t_survtmle <- survtmle_out$time


# ================
# onestep_all_t_fit <- onestep_all_t$survival_df
# sfun_onestep_all_t  <- stepfun(onestep_all_t_fit$T.uniq, c(1, onestep_all_t_fit$s_vec) , f = 1, right = TRUE)
