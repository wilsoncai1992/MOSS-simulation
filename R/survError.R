library(R6)

survError <- R6Class("survError",
  public = list(
    true_surv = NULL,
    surv_fit = NULL,
    T_range = NULL,
    mode = NULL,
    initialize = function(true_surv, object, T_range = 1:1e2, mode = NULL) {
      self$true_surv <- true_surv
      self$mode <- mode
      self$T_range <- T_range
      if ('survfit' %in% class(object)) {
        # KM output
        ha <- summary(object)
        surv1_km <- ha$surv[ha$strata == 'A=1']
        time1_km <- ha$time[ha$strata == 'A=1']
        self$surv_fit  <- stepfun(time1_km, c(1, surv1_km), f = 1, right = TRUE)
      }
      if ('data.frame' %in% class(object)){
        # survtmle output
        self$surv_fit  <- stepfun(object$time, c(1, object$s_1) , f = 1, right = TRUE)
      }
      if ('MOSS' %in% class(object)) {
        # onestep output
        if (mode == 'SL') self$surv_fit  <- stepfun(1:object$T.max, c(1, colMeans(object$Qn.A1.t_full)), f = 1, right = TRUE)
        if (mode == 'onestep') self$surv_fit  <- stepfun(1:object$T.max, c(1, object$Psi.hat), f = 1, right = TRUE)
      }
    },
    display = function(...){
      curve(self$surv_fit(x), from = 0, to = 1e2, ...)
    },
    calc_errors = function(){
      bias <- function(x) self$surv_fit(x) - self$true_surv(x)
      bias_out <- bias(self$T_range)
      mse <- function(x) (self$surv_fit(x) - self$true_surv(x))^2
      mse_out <- mse(self$T_range)
      variance_out <- mse_out - bias_out^2
      return(list(bias = bias_out, mse = mse_out, variance = variance_out))
    }
    )
)

error_SL <- survError$new(true_surv = true_surv, object = MOSS_fit, mode = 'SL')
error_onestep <- survError$new(true_surv = true_surv, object = MOSS_fit, mode = 'onestep')
error_KM <- survError$new(true_surv = true_surv, object = km.fit)
error_survtmle <- survError$new(true_surv = true_surv, object = survtmle_out)

error_survtmle$display()
yi <- error_survtmle$calc_errors()

# survError$new(p_density)