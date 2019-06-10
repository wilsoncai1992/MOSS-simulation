library(tidyverse)
library(ggpubr)

load("./output/df_metric.rda")
df_metric <- df_metric %>%
  filter(metric_name == "mse") %>%
  filter(method != "MOSS_classic")
df_metric$method[df_metric$method == "MOSS_hazard"] <- "one-step TMLE"
df_mse <- df_metric %>%
  group_by(method, t, n) %>%
  summarise(
    mse = mean(mse),
    bias = mean(bias),
    truth = mean(truth),
    cnt = dplyr::n()
  ) %>%
  mutate(
    variance = mse - bias ^ 2,
    rmse = sqrt(mse),
    nrmse = rmse / truth,
    bias_percent = bias / truth
  )
df_mse_tmle <- df_mse %>% filter(method == "TMLE")
df_mse_joined <- dplyr::left_join(df_mse, df_mse_tmle, by = c("t", "n"))
df_mse_joined <- df_mse_joined %>% mutate(re = mse.y / mse.x, method = method.x)

gg1 <- ggplot(df_mse, aes(x = t, y = sqrt(n) * bias, color = method)) +
  geom_line() +
  geom_hline(yintercept = 0, lty = 3) +
  ylab("sqrt(n) * bias") +
  ylim(c(NA, 10)) +
  facet_wrap(n ~ ., nrow = 1) +
  theme_bw() +
  guides(colour = guide_legend(nrow = 1))
gg2 <- ggplot(df_mse, aes(x = t, y = n * variance, color = method)) +
  geom_line() +
  scale_y_log10() +
  ylab("n * Variance") +
  facet_wrap(n ~ ., nrow = 1) +
  theme_bw()
gg3 <- ggplot(df_mse, aes(x = t, y = n * mse, color = method)) +
  geom_line() +
  ylab("n * MSE") +
  scale_y_log10() +
  facet_wrap(n ~ ., nrow = 1) +
  theme_bw()

ymax <- quantile(df_mse_joined$re[!is.na(df_mse_joined$re)], 0.99)
gg3_2 <- ggplot(df_mse_joined, aes(x = t, y = re, color = method)) +
  geom_line() +
  ylim(0, ymax) +
  ylab("Relative Efficiency") +
  facet_wrap(n ~ ., nrow = 1) +
  ylim(NA, 3) +
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
  legend = 'bottom',
  align = "v"
)

vals <- df_metric$bias * sqrt(df_metric$n)
xlims <- quantile(vals[!is.na(vals)], c(0.025, 0.975))
gg5 <- ggplot(
  df_metric %>% filter(t %% 10 == 1 & t <= 50),
  aes(x = bias * sqrt(n), color = method)
) +
  geom_density() +
  geom_vline(xintercept = 0, lty = 2) +
  xlim(xlims) +
  ylim(c(0, 10)) +
  facet_grid(t ~ n) +
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 1))

ggsave(gg_out1, filename = "./output/mse_panel1.png", width = 8, height = 8)
ggsave(gg5, filename = "./output/mse_panel2.png", width = 8, height = 8)
