# A simple helper that converts a fastFMM model object to plottable elements
plottable_fmm <- function(fmm, cov_names = NULL) {
  n_covs <- nrow(fmm$betaHat)

  # Get covariate names
  if (is.null(cov_names))
    cov_names <- rownames(fmm$betaHat)
  if (length(cov_names) != n_covs)
    stop(
      "Dimension mismatch:", "\n",
      "Provided covariate names are length ", length(cov_names), "; ",
      "model covariates are length ", n_covs
    )

  # Extract the coefficients
  res <- fmm$betaHat %>% t() %>% data.frame()
  names(res) <- cov_names
  res$s <- fmm$argvals
  res <- tidyr::pivot_longer(
    res,
    -s,
    names_to = "cov",
    values_to = "val"
  )

  # Extract the spacing for ci and joint_ci
  ci <- do.call(
    rbind,
    lapply(
      1:length(cov_names),
      function(x) {
        data.frame(
          s = fmm$argvals,
          cov = cov_names[x],
          ci = 1.96 * sqrt(diag(fmm$betaHat_var[, , x])),
          joint_ci = fmm$qn[x] * sqrt(diag(fmm$betaHat_var[, , x]))
        )
      }
    )
  )

  # Bind the CI df to the results and return it
  dplyr::left_join(res, ci, by = c("s", "cov"))
}

# A simple helper that converts a fastFMM model object to plottable elements
plottable_fmm_unsmooth <- function(fmm, cov_names = NULL) {
  n_covs <- nrow(fmm$betaTilde)

  # Get covariate names
  if (is.null(cov_names))
    cov_names <- rownames(fmm$betaTilde)
  if (length(cov_names) != n_covs)
    stop(
      "Dimension mismatch:", "\n",
      "Provided covariate names are length ", length(cov_names), "; ",
      "model covariates are length ", n_covs
    )

  # Extract the coefficients
  res <- fmm$betaTilde %>% t() %>% data.frame()
  names(res) <- cov_names
  res$s <- fmm$argvals
  res <- tidyr::pivot_longer(
    res,
    -s,
    names_to = "cov",
    values_to = "val"
  )

  # Extract the spacing for ci and joint_ci
  ci <- do.call(
    rbind,
    lapply(
      1:length(cov_names),
      function(x) {
        data.frame(
          s = fmm$argvals,
          cov = cov_names[x],
          ci = 1.96 * sqrt(diag(fmm$betaHat_var[, , x])),
          joint_ci = fmm$qn[x] * sqrt(diag(fmm$betaHat_var[, , x]))
        )
      }
    )
  )

  # Bind the CI df to the results and return it
  dplyr::left_join(res, ci, by = c("s", "cov"))
}
