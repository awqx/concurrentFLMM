#' Summarize across the iterations of simulated fits
#'
#' Given a parameter ID (i.e., a label corresponding to a file in `params/`), we
#' can pull the iterations of the simulation and generate summary statistics.
#'
#' There are three main groups of summaries to return. Columns passed to `joint`
#' are those where we want the proportion of rows that exceed a `cutoff` value.
#' Columns passed to `to_avg` only need the average report (e.g., pointwise
#' coverage or RMSE). The last category, `to_list`, are read directly from the
#' JSON of parameter settings.
#'
#' @param id Character, ID of parameters, passed without the file extension
#' (`".json"`)
#' @param joint Character vector, column names of parameters that should be
#' reported as the proportion above `cutoff`
#' @param to_avg Character vector, column names of parameters that should be
#' reported as their mean
#' @param to_list Character vector, column names of parameters that should be
#' directly read from the parameters JSON
#' @param cutoff Numeric, cutoff for `joint` statistics
#' @param params_dir Character, name of the parameters directory
#' @param read_dir Character, name of the directory containing the results to be
#' read

summarize_sims <- function(
  id,
  joint = c("beta_CI_joint", "beta_CI_joint_noIntrcpt"),
  to_avg = c("beta_CI_naive", "beta_CI_naive_noIntrcpt"),
  to_list = c("n_ids", "snr_sigma", "syn_functions"),
  cutoff = 1,
  params_dir = "params/",
  read_dir = "results/"
) {
  d <- read.csv(paste0(read_dir, id, ".csv"))
  json <- jsonlite::read_json(
    paste0(params_dir, id, ".json"),
    simplifyVector = T
  )
  res_params <- c(joint, to_avg, to_list)
  res <- vector("list", length(res_params))
  names(res) <- res_params

  dummy <- lapply(
    joint,
    function(x) res[[x]] <<- mean(d[[x]] >= cutoff)
  )
  dummy <- lapply(
    to_avg,
    function(y) res[[y]] <<- mean(d[[y]])
  )
  dummy <- lapply(
    to_list,
    function(z) res[[z]] <<- json[[z]]
  )

  res <- data.frame(res) %>%
    mutate(params_id = id)
}
