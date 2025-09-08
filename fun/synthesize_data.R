# fLME data simulation #########################################################

#' @param I (int) Number of clusters.
#' @param J Mean number of trials for each cluster.
#' @param snr_B Fixed to random effects ratio
#' @param snr_sigma Signal-noise ratio
#' @param ptwise_snr Boolean, `== TRUE` when `snr_B, snr_sigma` is applied to
#' each point on the functional domain separately (i.e., noise is independent
#' but not identically distributed).
#' @param x_name Name of functional covariate
#' @param return_fixed Whether to return the fixed effects matrix
#' @param return_noise whether to return the noise
#' @param asis_cols Whether to return functional covariates in matrix columns

synthesize_data <- function(
  I,
  J,
  snr_B,
  snr_sigma,
  ptwise_snr = T,
  y_name = "y",
  x_name = "x",
  return_fixed = F,
  return_noise = F,
  asis_cols = F,
  fix_J = T
) {
  if (fix_J) {
    J_list <- rep(J, I)
  } else {
    J_list <- rpois(I, J)
  }

  n <- sum(J_list)
  L <- length(beta_0())

  # Functional covariate generation
  X <- replicate(n, fun_cov())

  # Fixed effects
  fixed_fx <- X * beta_1()
  fixed_fx <- fixed_fx + beta_0()

  # Random effects
  g_0 <- do.call(
    cbind,
    lapply(
      1:I,
      function(i)
        replicate(J_list[i], gamma_0())
    )
  )
  g_1 <- do.call(
    cbind,
    lapply(
      1:I,
      function(i) {
        replicate(J_list[i], gamma_1())
      }
    )
  )
  if (!is.numeric(g_0)) {
    rand_fx <- as.numeric(g_0) + X * as.numeric(g_1)
    rand_fx <- matrix(rand_fx, ncol = n)
  } else {
    rand_fx <- g_0 + X * g_1
  }

  # Scale rand_fx
  if (ptwise_snr) {
    rand_scales <- sqrt(Rfast::rowVars(fixed_fx)) / sqrt(Rfast::rowVars(rand_fx)) / snr_B
    rand_fx <- rand_fx * rand_scales
  } else {
    rand_fx <- sd(fixed_fx) / sd(rand_fx) / snr_B * rand_fx
  }

  # Calculate noise
  # Add together signal and noise
  eta <- rand_fx + fixed_fx
  if (ptwise_snr) {
    # Replace eta with fixed_fx for sd calculation
    noise_sds <- sqrt(Rfast::rowVars(fixed_fx)) / snr_sigma
    noise <- replicate(n, rnorm(n = L, mean = 0, sd = noise_sds))
  } else {
    noise <- rnorm(length(eta), sd = sd(fixed_fx) / snr_sigma)
  }
  synth <- eta + noise

  if (asis_cols) {
    df <- data.frame(matrix(NA, nrow = ncol(X), ncol = 2))
    colnames(df) <- c(x_name, y_name)
    df[x_name] <- I(t(X))
    df[y_name] <- I(t(synth))
  } else {
    # Transpose, make a data.frame, add columns for trial and cluster
    df <- data.frame(t(X), t(synth))
    colnames(df) <- c(
      paste0(x_name, "_", 1:L),
      paste0(y_name, "_", 1:L)
    )
  }

  df$id <- as.factor(rep(1:I, J_list))
  df$trial <- sequence(J_list)

  # Return the fixed effects and the synthetic data
  if (!return_fixed) fixed_fx <- NULL
  if (!return_noise) noise <- NULL
  list(
    data = df,
    beta = rbind(beta_0(), beta_1()),
    fixed_fx = fixed_fx,
    noise = noise
  )
}
