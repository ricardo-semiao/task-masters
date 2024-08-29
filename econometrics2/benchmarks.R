
# PS6 - recursive ts creation ---------------------------------------------

bench::mark(
  {
    y1 <- rcauchy(n_obs, 0, 1)
    etas <- rep(runif(n_obs, 1, 5), 2)
    nus <- rnorm(n_obs*(n_time - 1), 0, 1)
    rests <- split_at(etas + nus, (1:(n_time - 1)) * n_obs)
    
    y <- accumulate(rests, \(y, r) a*y + r, .init = y1) %>%
      set_names(glue("y{1:3}"))
  },
  {
    y1 <- rcauchy(n_obs, 0, 1)
    etas <- runif(n_obs, 1, 5)
    nus <- rnorm(n_obs*(n_time - 1), 0, 1)
    
    y <- list(
      y1 = y1,
      y2 = y2 <- a*y1 + etas + nus[1:n_obs],
      y3 = a*y2 + etas + nus[(n_obs+1):2*n_obs]
    )
  },
  {
    eta <- runif(n = n_obs, min = 0, max = 5)
    nu2 <- rnorm(n = n_obs, mean = 0, sd = 1)
    nu3 <- rnorm(n = n_obs, mean = 0, sd = 1)
    y <- list(
      y1 = Y1 <- rcauchy(n = n_obs, location = 0, scale = 1),
      y2 = Y2 <- alpha * Y1 + eta + nu2,
      y3 = Y3 <- alpha * Y2 + eta + nu3
    )
  },
  check = FALSE
)