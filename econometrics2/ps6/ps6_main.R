
# Setup -------------------------------------------------------------------

try(setwd("ps6"), silent = TRUE)

library(glue)
library(lfe)
library(furrr)
library(tidyverse)

plan(multisession, workers = 7)


simulate_y <- function(a) {
  etas <- runif(n_obs, 1, 5)
  nus <- rnorm(n_obs*(n_time - 1), 0, 1)
  
  list(
    y1 = y1 <- rcauchy(n_obs, 0, 1),
    y2 = y2 <- a*y1 + etas + nus[1:n_obs],
    y3 = a*y2 + etas + nus[(n_obs+1):length(nus)]
  )
}

n_mc <- 10000
n_obs <- 10000
n_time <- 3



# Question 1 --------------------------------------------------------------

set.seed(220206)

a_q1 <- 0.2

results_q1 <- future_map_dfr(1:n_mc, function(rep) {
  y <- simulate_y(a_q1)
  
  data_1 <- tibble(
    y_ols = c(y$y2, y$y3), x_ols = c(y$y1, y$y2),
    y_with = y_ols - rowMeans(cbind(y$y2, y$y3)),
    x_with = x_ols  - rowMeans(cbind(y$y1, y$y2))
  )
  data_2 <- tibble(
    y_diff = y$y3 - y$y2, x_diff = y$y2 - y$y1,
    z_ah = y$y1
  )
  
  mods <- list(
    ols = felm(y_ols ~ x_ols - 1, data = data_1),
    with = felm(y_with ~ x_with - 1, data = data_1),
    diff = felm(y_diff ~ x_diff - 1, data = data_2),
    ah = felm(y_diff ~ -1 | 0 | (x_diff ~ z_ah), data = data_2)
  )
  
  map_dbl(mods, ~.x$coefficients)
},
  .options = furrr_options(seed = TRUE),
  .progress = TRUE
)

# Average bias:
colMeans(results_q1 - a_q1) %>% round(7)
#       ols       with       diff         ah 
# 0.0015171 -0.0002359 -0.0002359  0.0000003

# Proportion of WG estimator greater than OLS:
mean(results_q1$with - results_q1$ols > 0) %>% round(5)
# 0.09



# Question 2 --------------------------------------------------------------

set.seed(220206)

a_q2 = 0.9999

results_q2 <- future_map_dfr(1:n_mc, function(rep) {
  y <- simulate_y(a_q2)
  
  data <- tibble(
    y_diff = y$y3 - y$y2, x_diff = y$y2 - y$y1,
    z_ah = y$y1
  )
  
  mod <- felm(y_diff ~ -1 | 0 | (x_diff ~ z_ah), data = data)
  
  c(mod$stage1$tval^2, mod$coefficients) %>% set_names(c("tval", "coef"))
},
  .options = furrr_options(seed = TRUE),
  .progress = TRUE
)

# Deciles of the 1st stage F-stat:
quantile(results_q2$tval, probs = seq(0.1, 0.9, 0.1)) %>% round(3)
#   10%   20%   30%   40%   50%   60%   70%   80%   90% 
# 0.024 0.102 0.264 0.455 0.799 1.315 2.032 3.110 6.674

# Average bias:
mean(results_q2$coef - a_q2) %>% round(3)
# 0.027
