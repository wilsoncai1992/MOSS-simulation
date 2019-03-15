library(dplyr)
library(ggplot2)
library(ggpubr)


# load("./output/df_metric.rda")
# df_metric <- df_metric %>% filter(metric_name == "cross_entropy")
# df_metric$metric[is.na(df_metric$metric)] <- 1
# df_metric$metric[is.infinite(df_metric$metric)] <- 3
# df_agg <- df_metric %>%
#   group_by(t, n, method) %>%
#   summarise(metric = mean(metric))
#   # summarise(metric = median(metric))
# df_agg_tmle <- df_agg %>% filter(method == "TMLE")
# df_agg2 <- dplyr::left_join(df_agg, df_agg_tmle, by = c("t", "n"))
# df_agg2 <- df_agg2 %>% mutate(metric_relative = metric.x / metric.y, method = method.x)
# df_agg2$metric_relative[df_agg2$metric_relative >= 3] <- 3
# gg2 <- ggplot(df_agg2, aes(x = t, y = metric_relative, color = method)) +
#   geom_line() +
#   ylim(c(0, 3)) +
#   facet_grid(n ~ .)
# ggsave(gg2, filename = './output/cross_entropy_panel2.png', width = 12, height = 6)

# load("./output/df_metric.rda")
# df_metric <- df_metric %>% filter(metric_name == "cross_entropy")
# df_metric_tmle <- df_metric %>% filter(method == "TMLE")
# df_metric2 <- dplyr::left_join(df_metric, df_metric_tmle, by = c("id_mcmc", "t", "n"))
# df_metric2 <- df_metric2 %>%
#   mutate(metric_relative = metric.x / metric.y, method = method.x)
#   # mutate(metric_relative = metric.y / metric.x, method = method.x)
# df_metric2$metric_relative[is.na(df_metric2$metric_relative)] <- 1
# df_metric2$metric_relative[is.infinite(df_metric2$metric_relative)] <- 999
# metric_max <- 3
# metric_min <- 1e-3
# df_metric2$metric_relative[df_metric2$metric_relative >= metric_max] <- metric_max
# df_metric2$metric_relative[df_metric2$metric_relative <= metric_min] <- metric_min
# gg1 <- ggplot(df_metric2, aes(x = as.factor(t), y = metric_relative, color = method)) +
#   geom_boxplot() +
#   ylim(0, metric_max) +
#   facet_grid(n ~ .) +
#   theme_bw()
# df_metric3 <- df_metric2 %>%
#   group_by(method, t, n) %>%
#   summarise(
#     # metric_relative = median(metric_relative),
#     metric_relative = mean(metric_relative),
#     cnt = dplyr::n()
#   )
# gg2 <- ggplot(df_metric3, aes(x = t, y = metric_relative, color = method)) +
#   geom_line() +
#   facet_grid(n ~ .) +
#   theme_bw()
# gg3 <- ggplot(df_metric3, aes(x = t, y = cnt, color = method)) +
#   geom_line() +
#   facet_grid(n ~ .) +
#   theme_bw()
# gg_out1 <- ggarrange(
#   gg1,
#   gg2,
#   gg3,
#   ncol = 3,
#   common.legend = TRUE,
#   legend = 'bottom',
#   labels = "AUTO"
# )
# ggsave(gg_out1, filename = './output/cross_entropy_panel1.png', width = 12, height = 6)


load("./output/df_metric.rda")
df_metric <- df_metric %>% filter(metric_name == "mse")
df_mse <- df_metric %>%
  group_by(method, t, n) %>%
  summarise(
    mse = mean(mse),
    bias = mean(bias),
    # truth = mean(truth),
    cnt = dplyr::n()
  ) %>%
  mutate(
    variance = mse - bias ^ 2,
    rmse = sqrt(mse),
    # nrmse = rmse / truth,
    # bias_percent = bias / truth
  )
df_mse_tmle <- df_mse %>% filter(method == "TMLE")
df_mse_joined <- dplyr::left_join(df_mse, df_mse_tmle, by = c("t", "n"))
df_mse_joined <- df_mse_joined %>% mutate(re = mse.y / mse.x, method = method.x)
gg1 <- ggplot(df_mse, aes(x = t, y = bias, color = method)) +
  geom_line() +
  geom_hline(yintercept = 0, lty = 3) +
  ylab("Bias") +
  facet_wrap(n ~ ., nrow = 1) +
  theme_bw()
gg2 <- ggplot(df_mse, aes(x = t, y = variance, color = method)) +
  geom_line() +
  scale_y_log10() +
  ylab("Variance") +
  # ylim(0, 1) +
  facet_wrap(n ~ ., nrow = 1) +
  theme_bw()
gg3 <- ggplot(df_mse, aes(x = t, y = mse, color = method)) +
  geom_line() +
  ylab("MSE") +
  scale_y_log10() +
  # ylim(0, 1) +
  facet_wrap(n ~ ., nrow = 1) +
  theme_bw()

ymax <- quantile(df_mse_joined$re[!is.na(df_mse_joined$re)], 0.99)
gg3_2 <- ggplot(df_mse_joined, aes(x = t, y = re, color = method)) +
  geom_line() +
  ylim(0, ymax) +
  ylab("Relative Efficiency") +
  facet_wrap(n ~ ., nrow = 1) +
  theme_bw()
gg4 <- ggplot(df_mse, aes(x = t, y = cnt, color = method)) +
  geom_line() +
  ylab("Number of samples") +
  facet_wrap(n ~ ., nrow = 1) +
  theme_bw()
gg_out1 <- ggarrange(
  gg1,
  gg2,
  gg3,
  gg3_2,
  gg4,
  nrow = 5,
  common.legend = TRUE,
  legend = 'bottom'
)

vals <- df_metric$bias * sqrt(df_metric$n)
xlims <- quantile(vals[!is.na(vals)], c(0.025, 0.975))
rootroot <- function(x) x^(1/4)
irootroot <- function(x) x^4
gg5 <- ggplot(
  df_metric %>% filter(t %% 10 == 1 & t <= 50),
  aes(x = bias * sqrt(n), color = method)
) +
  geom_density() +
  geom_vline(xintercept = 0, lty = 2) +
  xlim(xlims) +
  # scale_y_continuous(trans=scales::trans_new("sq", rootroot, irootroot)) +
  ylim(c(0, 10)) +
  facet_grid(t ~ n)

ggsave(gg_out1, filename = "./output/mse_panel1.png", width = 10, height = 10)
ggsave(gg5, filename = "./output/mse_panel2.png", width = 10, height = 10)
