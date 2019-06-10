library(survival)
library(dplyr)
library(tidyr)
data(mgus2)

df <- mgus2
df <- df %>% mutate(is_male = as.numeric(sex == "M"))
# outcome
summary(df$ptime)
table(df$pstat)

# baseline
summary(df$age)
table(df$is_male)

# treatment
df$mspike_high <- as.numeric(df$mspike >= 1.5)
table(df$mspike_high)

Y_name <- c("ptime", 'pstat')
A_name <- c("mspike_high")
W_name <- c("age", "is_male", "hgb", "creat")

df_train <- df %>%
  select(c(Y_name, A_name, W_name)) %>%
  mutate(ptime = floor(ptime / 20) + 1) %>%
  drop_na()

# don't analyze t with high censoring probability
# t_censor <- 8
is_censored <- df_train$ptime >= t_censor
df_train$pstat[is_censored] <- 0
df_train$ptime[is_censored] <- t_censor
