load('./scenario.rda')
dir.create('./plot/')
# T_max <- 10
T_max <- 20
# create df
df_bias <- data.frame(time = 1:length(error_KM$bias),
                      KM = error_KM$bias,
                      survtmle = error_survtmle$bias,
                      SL = error_SL$bias,
                      onestep = error_onestep$bias)
df_variance <- data.frame(time = 1:length(error_KM$variance),
                          KM = error_KM$variance,
                          survtmle = error_survtmle$variance,
                          SL = error_SL$variance,
                          onestep = error_onestep$variance)
df_mse <- data.frame(time = 1:length(error_KM$mse),
                     KM = error_KM$mse,
                     survtmle = error_survtmle$mse,
                     SL = error_SL$mse,
                     onestep = error_onestep$mse)
library(tidyr)
gather_into_long <- function(df) {
  df2 <- gather(data = df[,2:5], 'method', 'value')
  df2$method[df2$method == 'SL'] <- 'initial fit'
  df2$method[df2$method == 'survtmle'] <- 'iterative TMLE'
  df2$method[df2$method == 'onestep'] <- 'one-step TMLE'
  df2$method[df2$method == 'KM'] <- 'Kaplan-Meier'
  df2$time <- df$time
  return(df2)
}

library(dplyr)
df_relative_efficiency <- df_mse
df_relative_efficiency <- df_relative_efficiency %>% mutate(onestep = onestep/survtmle, SL = SL/survtmle, KM = KM/survtmle) %>% mutate(survtmle = 1)
df_relative_efficiency$onestep[df_mse$survtmle < 1e-5] <- 1
df_relative_efficiency[is.na(df_relative_efficiency)] <- 1
df_relative_efficiency$KM[df_relative_efficiency$KM>3] <- 3
df_relative_efficiency$SL[df_relative_efficiency$SL>3] <- 3
df_relative_efficiency$onestep[df_relative_efficiency$onestep>3] <- 3

df_bias2 <- gather_into_long(df_bias)
df_variance2 <- gather_into_long(df_variance)
df_mse2 <- gather_into_long(df_mse)
df_relative_efficiency2 <- gather_into_long(df_relative_efficiency)
# bias
library(ggplot2)
plot1 <- ggplot(data = df_bias2, aes(x = time, y = value, color = method)) +
  geom_line() +
  geom_hline(yintercept = 0, lty = 2) +
  ylab('bias') +
  xlim(c(0,T_max)) +
  theme_bw() +
  theme(legend.position="bottom") +
  guides(color=guide_legend(nrow=2))
ggsave(filename=file.path('./plot/', "bias.pdf"), plot=plot1, width=4, height=4)

# variance
plot2 <- ggplot(data = df_variance2, aes(x = time, y = value, color = method)) +
  geom_line() +
  geom_hline(yintercept = 0, lty = 2) +
  ylab('variance') +
  xlim(c(0, T_max)) +
  theme_bw() +
  theme(legend.position="bottom") +
  guides(color=guide_legend(nrow=2))
ggsave(filename=file.path("./plot/", "variance.pdf"), plot=plot2, width=4, height=4)
# mse
plot3 <- ggplot(data = df_mse2, aes(x = time, y = value, color = method)) +
  geom_line() +
  geom_hline(yintercept = 0, lty = 2) +
  ylab('MSE') +
  xlim(c(0,T_max)) +
  theme_bw() +
  theme(legend.position="bottom") +
  guides(color=guide_legend(nrow=2))
ggsave(filename=file.path("./plot/", "MSE.pdf"), plot=plot3, width=4, height=4)

# relative efficiency
plot4 <- ggplot(data = df_relative_efficiency2, aes(x = time, y = value, color = method)) +
  geom_line() +
  geom_hline(yintercept = 1, lty = 2) +
  ylab('relative efficiency (compared to iterative TMLE)') +
  xlim(c(0,T_max)) +
  ylim(c(0, 3)) +
  theme_bw() +
  theme(legend.position="bottom") +
  guides(color=guide_legend(nrow=2))
ggsave(filename=file.path("./plot/", "relative_efficiency.pdf"), plot=plot4, width=4, height=4)

