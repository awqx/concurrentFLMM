#' Generate a long table of true fixed effects
#'
#' To evaluate coverage and visually inspect model fits, it's helpful to have
#' the true underlying fixed coefficients handy. This function converts the data
#' generating functions (the scripts in `syn/`) to long tables with the location
#' on the functional domain, the coefficient (`beta_0` or `beta_1`) and the true
#' value.
#'
#' Current functionality only lists `beta_0` and `beta_1`.
#'
#' @param syn_fun Character, the name of the file (with or without `".R"`)
#' containing the scripts of interest.

truth_from_syn <- function(syn_fun) {
  # Check the file name
  syn_fun <- ifelse(
    grepl("\\.R$", syn_fun),
    syn_fun,
    paste0(syn_fun, ".R")
  )

  # Source the fixed effects
  source(paste0("fun/", syn_fun))

  # Generate values and bind them in a long table
  truth_df <- rbind(
    data.frame(
      cov = "beta_0",
      truth_val = beta_0(),
      s = 1:L
    ),
    data.frame(
      cov = "beta_1",
      truth_val = beta_1(),
      s = 1:L
    )
  ) %>%
    mutate(syn_function = syn_fun)
}
