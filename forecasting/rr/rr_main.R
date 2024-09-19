
# Packages and Functions --------------------------------------------------

# For models:
library(HDeconometrics) #for ic.glmnet, boosting, bagging
library(vars) #for VAR, VARselect
library(randomForest) #for randomForest
library(forecast) #for auto.arima

# Aesthetics:
library(varr)
library(gt)

# Basics:
library(tidyverse)
library(furrr)
library(lubridate)

plan(multisession, workers = 7)
theme_set(theme_bw())

source("rr/models.R")
source("rr/results.R")

pca <- function(data, n_lags, k = 1:(n_lags - 1)) {
  names <- pmap_chr(expand_grid("comp", k, "_lag", 0:(n_lags - 1)), paste0)
  pca <- princomp(scale(data, scale = FALSE)) #scale = TRUE
  
  pca$scores[, k] %>%
    embed(n_lags) %>%
    `colnames<-`(names)
}



# Data --------------------------------------------------------------------

utils <- list()

# Loading data:
#load("rr/rawdata.rda")
#load("rr/rawdata.Rdata")
data_raw <- read_delim("rr/dataset_fred_m.csv", delim = ";")

# Removal of incomplete variables:
utils$cols_rmv <- map_dbl(data_raw, ~ sum(is.na(.x))) %>%
  .[order(.)] %>%
  .[. >= 0.05 * nrow(data_raw)] %>%
  names()

# Transformation plans for each variable:
utils$trans <- data_raw[1,-1] %>%
  unlist() %>%
  split(names(.), .) %>%
  set_names(~ paste0("tcode", .x)) %>%
  map(~ setdiff(.x, utils$cols_rmv))

# Chosen y and training window size:
y_name <- "CPIAUCSL"
n_train <- 132


# Data processing:
data_base <- data_raw %>%
  #dealing with NA's:
  slice(-(1:12)) %>% #first is "transform" row, 2:13 are from Permit's columns NA's
  select(-all_of(utils$cols_rmv)) %>%
  fill(-sasdate, .direction = "down") %>% 
  #processing variables:
  mutate(
    #dates and dummy of post november 2009:
    sasdate = as.Date(sasdate, "%m/%d/%Y"),
    POST2009 = as.integer(sasdate >= as.Date("2009-11-01")),
    #variables transformations:
    across(utils$trans$tcode2, ~ c(NA, diff(.x))),
    across(utils$trans$tcode4, ~ log(.x)),
    across(utils$trans$tcode5, ~ c(NA, diff(log(.x)))),
    across(utils$trans$tcode6, ~ c(NA, NA, diff(log(.x), 2))),
    across(utils$trans$tcode7, ~ .x/lag(.x) - 1)
  ) %>%
  #removing NAs created when transforming, and dropping 11 variables as in their paper
  slice(-(1:2)) %>%
  #relocating y candidates:
  relocate(all_of(c("CPIAUCSL", "RPI", "INDPRO")), .after = 1) #PCEPI

data_split <- list(
  full = data_base,
  post_pandemic = filter(data_base, sasdate >= as.Date("2019-12-01") %m-% months(n_train))
) %>%
  pluck("post_pandemic")

# Creating X and Y data subsets:
n_lags <- 5
n_lags_x <- 5
utils$adj_rows <- if (n_lags - n_lags_x > 0) -(1:(n_lags - n_lags_x)) else TRUE

utils$data_x <- select(data_split, -all_of(c("sasdate", y_name)))

data <- list()

data$y <- data_split %>%
  slice(-(1:(n_lags - 1))) %>%
  pull(y_name)

data$std <- cbind(
  embed(data_split[[y_name]], n_lags)[, -1] %>%
    `colnames<-`(paste0(y_name, 1:(n_lags - 1))),
  slice(utils$data_x, -(1:(n_lags - 1)))
) %>%
  as.matrix()

data$lag <- cbind(
  data$std[, 1:(n_lags - 1)],
  embed(as.matrix(utils$data_x), n_lags_x)[utils$adj_rows,]
) %>%
  `colnames<-`(c(
    paste0(y_name, 1:(n_lags - 1)),
    expand_grid(names(utils$data_x), 1:n_lags_x) %>% pmap_chr(paste0)
  )) %>%
  as.matrix()

data$pca <- pca(data_split[-1], n_lags)


# Utility values for models:
utils$lasso_relevant <- model$shrink(data)$adalasso %>%
  coef() %>%
  .[-c(1, which(grepl(y_name, names(.), fixed = TRUE)))] %>%
  .[. > 0] %>%
  names() %>%
  #gsub("^(.+)[0-9]$", "\\1", .) %>%
  unique()

utils$target_select <- HDeconometrics:::baggit(as.matrix(cbind(data$y, data$std)),
  pre.testing = "individual", fixed.controls = 2:5
) %>%
  {which(.[-(1:5)] != 0)}

data$pca_target <- pca(data_split[-1][utils$target_select], n_lags)



# Main Loop ---------------------------------------------------------------

horizons <- set_names(1:12)
window_type <- "rolling" #rolling or fixed
windows <- seq(1, length(data$y) - n_train - max(horizons) - 1, 1)

predictions <- map(horizons, function(h) {
  h_preds <- future_map_dfr(windows, .options = furrr_options(seed = TRUE), function(w) {
    utils <- utils; ic.glmnet <- ic.glmnet;  randomForest <- randomForest #to solve furrr bug
    start <- switch(window_type, "fixed" = 1, "rolling" = w)
    
    data_adj <- imap(data, function(x, name) {
      if (name == "y") {
        x[(h + 1):length(x)][start:(n_train + w)]
      } else {
        x[1:(nrow(x) - h), ][start:(w + n_train), ]
      }
    })
    
    data_pred <- imap(data_adj, function(x, name) {
      if (name == "y") x[length(x), drop = FALSE] else x[nrow(x), , drop = FALSE]
    })
    
    mods <- c(
      model$bench(data_adj),
      model$shrink(data_adj),
      model$factor(data_adj),
      model$ensemble(data_adj)#, model$var(data_adj)
    )
    
    map_dbl(set_names(names(mods)), \(name) pred[[name]](mods[[name]], data_pred))
  })
  
  h_preds %>%
    mutate(
      mean = rowMeans(.),
      tmean = apply(., 1, mean, trim = 0.05),
      median = apply(., 1, median)
    )
})

errors <- imap_dfr(pred_cpi_fixed, function(h_preds, h) {
  h_preds %>%
    mutate(across(-(1:2), \(x) x - data_split[[y_name]][windows])) %>%
    imap_dfr(function(x, name) {
      
      list(
        Horizon = as.integer(h),
        Model = name,
        RMSE = sqrt(mean(x^2)),
        MAE = mean(abs(x))#, MAD = median(abs(x - median(x)))
      )
    }) %>%
    mutate(across(-c(Horizon, Model), \(x) x / x[1]))
})


# Presentation ------------------------------------------------------------

utils$model_labs <- set_names(
  c(
    "RW", "AR", "PCA→ARDL",
    "LASSO", "ELNET", "adaLASSO", "adaELNET",
    "targetPCA",
    "CSR", "Forest",
    "Mean", "T. Mean", "Median"
  ),
  c(
    "rw", "ar", "ardlpca",
    "lasso", "elnet", "adalasso", "adaelnet",
    "target", "csr",
    "forest",
    "mean", "tmean", "median"
  )
)

table_errors_avg(errors, utils$model_labs, latex = TRUE)

table_errors_h(errors, utils$model_labs, latex = TRUE)


table_errors_avg(errors, utils$model_labs)

g = predictions$`1` %>%
  mutate(
    Window = windows,
    True = data_split[[y_name]][windows],
    .before = 1
  ) %>%
  pivot_longer(-c(Window, True))

ggplot(g, aes(Window, value, colors = name)) +
  geom_line(linetype = "dashed") +
  geom_line(aes(y = True))


# Exploratory Analisys ----------------------------------------------------

ggvar_history(data_base, "CPIAUCSL", index = data_base$sasdate)

ggvar_acf(data_base, "CPIAUCSL", lag.max = 700, graph_type = "area")
ggvar_acf(data_base, "CPIAUCSL", lag.max = 700, type = "partial", graph_type = "area")

aTSA::adf.test(data_base$CPIAUCSL, output = FALSE)[[2]]

ggvar_distribution(data_base, "CPIAUCSL", plot_normal = FALSE)

strucchange::breakpoints(data$y ~ data$pca) %>%
  {
    print(summary(.))
    plot(.)
  }



# Diagnostics -------------------------------------------------------------

auto.arima(data$y, max.q = 0, d = 0, seasonal = FALSE, ic = "bic")

VARselect(cbind(data$y, data$std[, utils$lasso_relevant]), n_lags) %>%
  ggvar_select()

VARselect(cbind(data$y, data$pca[, grepl("0$", colnames(data$pca))]), n_lags,) %>%
  ggvar_select()
