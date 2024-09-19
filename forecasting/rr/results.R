# Utility functions:
sort_names <- function(horizons) {
  map(1:length(horizons), ~ .x + c(1, max(horizons) + 1)) %>%
    list_c() %>%
    c(1, .)
}

get_spanners <- function(horizons, ba = c("(", ")")) {
  map(horizons, function(h) {
    \(x) tab_spanner(x, paste0(ba[1], h, ba[2]), paste0(c("RMSE_", "MAE_"), h))
  })
}

get_row_groups <- function(indexes = NULL) {
  indexes <- indexes %||% list(c(11:13),c(8:10), c(4:7), c(1:3))
  map2(indexes, 1:4, function(ind, id) {
    \(x) tab_row_group(x, "", ind, id = id)
  })
}

get_labels <- function(horizons) {
  labels <- rep(c("RMSE", "MAE"), length(horizons))
  names(labels) <- map(horizons, ~ paste0(c("RMSE", "MAE"), "_", .x)) %>% list_c()
  as.list(labels)
}


# Table creation:
table_errors_avg <- function(errors, mod_labs, latex = FALSE) {
  acum_funs <- list(
    `Acum. 3m` = \(x) round(sum(x[1:3]) / 3, 2),
    `Acum. 6m` = \(x) round(sum(x[1:6]) / 6, 2),
    `Acum. 12m` = \(x) round(sum(x[1:12]) / 12, 2)
  )
  
  table <- errors %>%
    mutate(Model = fct_relabel(fct(Model), ~ mod_labs[.x])) %>%
    group_by(Model) %>%
    summarize(across(c(RMSE, MAE), acum_funs, .names = "{.col}_{.fn}")) %>%
    gt() %>%
    reduce(get_spanners(names(acum_funs), ba = c("", "")), \(x, f) f(x), .init = .) %>%
    cols_label(.list = get_labels(names(acum_funs))) %>%
    reduce(get_row_groups(), \(x, f) f(x), .init = .)
  
  if (latex) cat(as_latex(table)[[1]]) else table
}


table_errors_h <- function(errors, mod_labs, latex = FALSE) {
  table <- errors %>%
    mutate(Model = fct_relabel(fct(Model), ~ mod_labs[.x])) %>%
    pivot_wider(names_from = Horizon, values_from = c(RMSE, MAE)) %>%
    .[sort_names(horizons)] %>%
    mutate(across(-Model, ~ round(.x, 2))) %>%
    gt() %>%
    tab_spanner("Horizon", Model) %>%
    reduce(get_spanners(horizons), \(x, f) f(x), .init = .) %>%
    cols_label(.list = get_labels(horizons)) %>%
    reduce(get_row_groups(), \(x, f) f(x), .init = .)
  
  if (latex) cat(as_latex(table)[[1]]) else table
}
