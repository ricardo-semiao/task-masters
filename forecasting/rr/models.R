
# Models ------------------------------------------------------------------

model <- list()

# Benchmark:
model$bench <- function(data, select = utils$lasso_relevant) {
  res <- list()
  res$rw <- arima(data$y, c(0, 1, 0))
  res$ar <- lm(data$y ~ ., as.data.frame(data$std[,1:4]))
  #res$ardllasso <- lm(data$y ~ .,
  # as.data.frame(cbind(data$std[,1:4], data$std[, select]))
  #)
  res$ardlpca <- lm(data$y ~ .,
    as.data.frame(cbind(data$std[,1:4], data$pca[, grepl("0$", colnames(data$pca))]))
  )
  res
}

# Shrinkage:
model$shrink <- function(data, x_which = "std", run_ridge = FALSE, elnet_alpha = 0.5) {
  ridge <- if (run_ridge) ic.glmnet(data$std, data$y, alpha = 0)
  
  std <- c(lasso = 1, elnet = elnet_alpha) %>%
    map(function(alpha) {
      ic.glmnet(data[[x_which]], data$y, alpha = alpha)
    })
  
  ada <- c(adalasso = 1, adaelnet = elnet_alpha) %>%
    map2(std[-1], function(alpha, mod) {
      penalty <- (abs(coef(mod)[-1]) + 1 / sqrt(length(data$y)))^(-1)
      ic.glmnet(data[[x_which]], data$y, penalty.factor = penalty, alpha = alpha)
    })
  
  c(ridge = ridge, std, ada)
}

# Factor:
model$factor <- function(data, run_boost = FALSE) {
  res <- list()
  res$target <- lm(data$y ~ ., as.data.frame(data$pca_target))
  if (run_boost) res$boosting <- boosting(data$pca, data$y)
  res
}

# Ensemble:
model$ensemble <- function(data, run_bag = FALSE) {
  res <- list()
  if (run_bag) {
    res$bagging <- bagging(as.matrix(data$pca), data$y, l = 5, pre.testing = "group-joint")
  }
  res$csr <- csr(data$pca, data$y, K = 10, fixed.controls = seq(1, ncol(data$pca), 5))
  res$forest <- randomForest(data$pca, data$y, importance = TRUE)
  res
}

# VAR:
model$var <- function(data, ps = c(3, 3), select = utils$lasso_relevant) {
  res <- list()
  res$lasso <- VAR(cbind(data$y, data$std[, select]), p = ps[1])
  res$varpca <- VAR(cbind(y = data$y, data$pca[, grepl("0$", colnames(data$pca))]), p = ps[2])
  res
}



# Predictions -------------------------------------------------------------


pred <- list()

pred$rw <- \(x, data) predict(x, 1)$pred
pred$ar <- \(x, data) predict(x, as.data.frame(as.list(data$std[, 1:4])))
pred$ardllasso <- \(x, data, select = utils$lasso_relevant) {
  predict(x,
    as.data.frame(as.list(c(data$std[,1:4], data$std[, select])))
  )
}
pred$ardlpca <- \(x, data) {
  predict(x,
    as.data.frame(as.list(c(data$std[,1:4], data$pca[, grepl("0$", colnames(data$pca))])))
  )
}

#pred$ridge
pred$lasso <- \(x, data) predict(x, data$std)[1,1]
pred$elnet <- pred$adalasso <- pred$adaelnet <- pred$lasso

pred$target <- \(x, data) predict(x, as.data.frame(data$pca_target))
#pred$boosting

#pred$bagging
pred$csr <- \(x, data) predict(x, as.matrix(data$pca))
pred$forest <- \(x, data) predict(x, as.matrix(data$pca))

#pred$varlasso <- \(x, data) {
#  endog <- cbind(data$y, data$std[utils$lasso_relevant])
#  
#  cbind(lag(endog), lag(endog, 2), lag(endog, 3), 1) %>%
#    na.omit() %>%
#    `%*%`(., coef(x$varresult$y))
#}
#pred$varpca <- \(x, data) {
#  endog <- c(data$y, data$pca[, grepl("0$", colnames(data$pca))])
#  
#  cbind(lag(endog), lag(endog, 2), lag(endog, 3), 1) %>%
#    na.omit() %>%
#    `%*%`(., coef(x$varresult$y))
#}

