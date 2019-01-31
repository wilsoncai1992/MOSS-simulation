simulate_data <- function(n_sim = 2e2) {
  library(simcausal)
  D <- DAG.empty()
  D <- D +
    node("W1", distr = "rbinom", size = 1, prob = .5) +
    node("W", distr = "runif", min = 0, max = 1.5) +
    node("A", distr = "rbinom", size = 1, prob = .4 + .5*as.numeric(W > .75)) +
    node("Trexp", distr = "rlnorm", meanlog = 2 - 1*W + 1*A, sdlog = .01) +
    node("Cweib", distr = "rconst", const = 200) +
    node("T", distr = "rconst", const = round(Trexp*1)) +
    node("C", distr = "rconst", const = round(Cweib*1)) +
    # Observed random variable (follow-up time):
    node("T.tilde", distr = "rconst", const = ifelse(T <= C, T, C)) +
    # Observed random variable (censoring indicator, 1 - failure event, 0 - censored):
    node("Delta", distr = "rconst", const = ifelse(T <= C, 1, 0))
  setD <- set.DAG(D)
  dat <- sim(setD, n = n_sim)
  # only grab ID, W's, A, T.tilde, Delta
  Wname <- grep('W', colnames(dat), value = TRUE)
  dat <- dat[,c('ID', Wname, 'A', "T.tilde", "Delta")]

  # input: scalar q, W vector. computes for all W, the S(q|A,W)

  true_surv_one <- function(q, W, A = 1) sapply(W, function(w) {
    1 - plnorm(q, meanlog = 2 - 1*w + 1*A, sdlog = .01)
  })
  # input: vector q. mean(S(q|A,W)|A), average out W. loop over q
  true_surv <- function(q_grid, surv_fn, A) {
    W_grid <- seq(0, 1.5, .01)
    survout <- numeric()
    for (q in q_grid) survout <- c(survout, mean(surv_fn(q = q/1, W = W_grid, A = A)))
    # for (q in q_grid) survout <- c(survout, mean(surv_fn(q = (q-1)/.2, W = W_grid, A = A)))
    return(survout)
  }
  truth_surv <- function(q) true_surv(q_grid = q, surv_fn = true_surv_one, A = 1)
  truth_surv0 <- function(q) true_surv(q_grid = q, surv_fn = true_surv_one, A = 0)
  return(list(dat = dat, true_surv1 = truth_surv, true_surv0 = truth_surv0))
}
