load('./scenario.rda')
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
# bias
library(ggplot2)
plot1 <- ggplot(data = df_bias, aes(x = time)) +
  geom_line(aes(y = KM, colour = 'KM')) +
  geom_line(aes(y = survtmle, colour = 'survtmle')) +
  geom_line(aes(y = SL, colour = 'SL')) +
  geom_line(aes(y = onestep, colour = 'onestep')) + 
  geom_hline(yintercept = 0, lty = 2) +
  ylab('bias') + 
  labs(colour='estimator')+
  xlim(c(0,20))
ggsave(filename="bias.pdf", plot=plot1, width=8, height=6)

# variance
plot2 <- ggplot(data = df_variance, aes(x = time)) +
  geom_line(aes(y = KM, colour = 'KM')) +
  geom_line(aes(y = survtmle, colour = 'survtmle')) +
  geom_line(aes(y = SL, colour = 'SL')) +
  geom_line(aes(y = onestep, colour = 'onestep')) + 
  geom_hline(yintercept = 0, lty = 2) +
  ylab('variance') + 
  labs(colour='estimator') +
  xlim(c(0,20))
ggsave(filename="variance.pdf", plot=plot2, width=8, height=6)
# mse
plot3 <- ggplot(data = df_mse, aes(x = time)) +
  geom_line(aes(y = KM, colour = 'KM')) +
  geom_line(aes(y = survtmle, colour = 'survtmle')) +
  geom_line(aes(y = SL, colour = 'SL')) +
  geom_line(aes(y = onestep, colour = 'onestep')) + 
  geom_hline(yintercept = 0, lty = 2) +
  ylab('MSE') + 
  labs(colour='estimator') +
  xlim(c(0,20))
ggsave(filename="MSE.pdf", plot=plot3, width=8, height=6)
# relative efficiency

