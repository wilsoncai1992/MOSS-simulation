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
MOSS_fit <- MOSS::MOSS$new(dat = df, dW = 1, epsilon.step = 1e-1, max.iter = 5e2, verbose = FALSE)
# MOSS_fit <- MOSS::MOSS$new(dat = df, dW = 1, epsilon.step = 1e0, max.iter = 1e2, verbose = FALSE, tol = 1e-1/nrow(df))
# MOSS_fit <- MOSS::MOSS$new(dat = df, dW = 1, epsilon.step = 1e1, max.iter = 5e1, verbose = FALSE)
MOSS_fit$onestep_curve(g.SL.Lib = c("SL.mean","SL.glm",'SL.gam'),
                       Delta.SL.Lib = c("SL.mean","SL.glm",'SL.gam'),
                       ht.SL.Lib = c("SL.mean","SL.glm",'SL.gam'))
error_SL <- survError$new(true_surv = true_surv, object = SL_fit, mode = 'SL')
error_onestep <- survError$new(true_surv = true_surv, object = MOSS_fit, mode = 'onestep')
error_KM <- survError$new(true_surv = true_surv, object = km.fit)

# survtmle
survtmle_out <- fit_survtmle(dat = df,Wname = c('W', 'W1'))
s_0 <- survtmle_out$s_0
s_1 <- survtmle_out$s_1
t_survtmle <- survtmle_out$time
out <- survError$new(true_surv = true_surv, object = survtmle_out)
# lines(s_1 ~ t_survtmle, col = 'red', lty = 1) #WILSON
out

df1 <- data.frame(s = s_1, time = t_survtmle, method = 'iterative TMLE')
df2 <- data.frame(s = MOSS_fit$Psi.hat, time = 1:MOSS_fit$T.max, method = 'one-step TMLE')
df3 <- data.frame(s = colMeans(SL_fit$Qn.A1.t_full), time = 1:SL_fit$T.max, method = 'initial fit')
df4 <- data.frame(s = tail(km.fit$surv, km.fit$strata['A=1']), time = tail(km.fit$time, km.fit$strata['A=1']), method = 'Kaplan-Meier')
df5 <- data.frame(s = true_surv(1:20), time = 1:20, method = 'truth')
# df_realization <- rbind(df1, df2, df4)
df_realization <- rbind(df1, df2, df3, df4, df5)
library(ggplot2)
ggplot(df_realization, aes(x = time, y = s, color = method)) + 
  geom_line() +
  xlim(c(1, 10)) +
  ylab('survival probability') + 
  theme_bw() + 
  theme(legend.position = 'bottom') + 
  guides(color=guide_legend(nrow=2))
