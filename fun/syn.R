# Read a ground truth model, simulate a draw, fit a model, then report

# 0 Load and setup #############################################################

# Rstudio mode allows for testing without the command line
rstudio_mode <- F
rstudio_id <- 1
rstudio_params <- "syn_256"

# 0.0 Packages =================================================================

suppressMessages(library(dplyr))

# 0.1 Command-line args ========================================================

args <- commandArgs(TRUE)
iter <- as.integer(as.numeric(args[1]))

if (rstudio_mode) {
  iter <- rstudio_id
}

if (is.na(iter)) {
  stop(
    "Argument `iter` not detected.", "\n",
    "Please ensure iter and params_id are given as command-line args."
  )
}

params_id <- args[2]

if (rstudio_mode) params_id <- rstudio_params

# 0.2 Directories ==============================================================

# End directories with slash
base_dir <- "~/"
params_dir <- paste0(base_dir, "params/elife/")
fun_dir <- paste0(base_dir, "fun/")
syn_dir <- paste0(base_dir, "syn/")
save_dir <- paste0(base_dir, "raw/outs/", params_id, '/')
coefs_dir <- paste0(base_dir, "raw/coefs/", params_id, "/")

if (!dir.exists(save_dir)) dir.create(save_dir, recursive = T)
if (!dir.exists(coefs_dir)) dir.create(coefs_dir, recursive = T)

# 1 Simulation parameters ######################################################

# 1.1 Consistent across runs ==================================================+

# Testing for existence
params_json <- tryCatch(
  jsonlite::read_json(
    paste0(params_dir, params_id, ".json"),
    simplifyVector = T
  ),
  error = function(e) {
    message("Params JSON not found in the specified path.")
    message(paste("Check", params_id, "is in directory", params_dir))
  }
)

res_names <- c(
  "params_id",
  "synth_time",
  "fLME_time",
  "beta_CI_joint",
  "beta_CI_joint_noIntrcpt",
  "beta_CI_joint0",
  "beta_CI_naive",
  "beta_CI_naive_noIntrcpt",
  "fLME_beta_rmse",
  "fLME_beta_rmse_noIntrcpt",
  "fLME_bias"
  # "n_trials",
  # "n_ids",
  # "knots_div",
  # "smooth_method",
  # "iter"
)

res <- data.frame(matrix(NA, nrow = 1, ncol = length(res_names)))
colnames(res) <- res_names
res$params_id <- params_id

# 2 Data generation ############################################################

# 2.1 Load functions ===========================================================

# Check for the correct file extension for the sourced functions
syn_fun <- params_json$syn_functions
syn_fun <- ifelse(
  grepl("\\.R$", syn_fun),
  syn_fun,
  paste0(syn_fun, ".R")
)

sim_fun <- params_json$sim_functions
sim_fun <- ifelse(
  grepl("\\.R$", sim_fun),
  sim_fun,
  paste0(sim_fun, ".R")
)

source(paste0(syn_dir, syn_fun))
source(paste0(fun_dir, sim_fun))
message(
  "Sourcing functions from ",
  params_json$syn_functions, " and ",
  params_json$sim_functions
)

# 2.2 Data simulation ==========================================================

set.seed(iter)

timeStart <- Sys.time()

fixed_noise <- params_json$fixed_noise
if (is.null(fixed_noise)) fixed_noise <- F

if (fixed_noise) {
  dat_sim <- synthesize_data(
    I = params_json$n_ids,
    J = params_json$mean_trials,
    rand_scale = params_json$rand_scale,
    sigma_eps = params_json$sigma_eps
  )
} else {
  dat_sim <- synthesize_data(
    I = params_json$n_ids,
    J = params_json$mean_trials,
    snr_B = params_json$snr_B,
    snr_sigma = params_json$snr_sigma
  )
}

timeEnd <- Sys.time()
message('Finished synthesizing data.')

res$synth_time <- as.numeric(difftime(timeEnd, timeStart, units = "mins"))

# 3 Fit fastFMM model #########################################################

# 3.1 Model fitting ============================================================

# number of knots
L <- length(beta_0())
nknots_min <- round(L / params_json$knots_div)
# nknots_min_cov <- round(L / 4)
n <- length(unique(dat_sim$data$id))

# fit model and get model CIs
timeStart <- Sys.time()
fit <- fastFMMconc::fui(
  formula = as.formula(params_json$formula),
  data = dat_sim$data,
  family = params_json$family,
  var = params_json$var,
  analytic = params_json$analytic,
  parallel = params_json$parallel,
  silent = params_json$silent,
  smooth_method = params_json$smooth_method,
  seed = params_json$seed,
  n_cores = params_json$n_cores,
  non_neg = params_json$non_neg,
  MoM = params_json$MoM,
  concurrent = params_json$concurrent,
  nknots_min_cov = params_json$nknots_min_cov
  # unsmooth = F
)
timeEnd <- Sys.time()
message('Finished fitting new models.')

res$fLME_time <- as.numeric(difftime(timeEnd, timeStart, units = "mins"))
message(
  "Finished fitting fastFMM models. Time (min): ", signif(res$fLME_time, 3)
)

# 3.2 Analysis =================================================================

beta_idx <- 2
# bias of cue term
res$fLME_bias <- mean(
  dat_sim$beta[beta_idx, ] - fit$betaHat[beta_idx, ]
)

# sample size of sims (i.e., number of trials analyzed)
# res$n_trials <- nrow(dat_sim$data)

# 3.2.1 Generate 95% CIs -------------------------------------------------------

# initialize CIs
lower <- upper <- lower_joint <- upper_joint <- matrix(
  NA, nrow = nrow(fit$betaHat), ncol = ncol(fit$betaHat)
)
joint_incl <- naive_incl <- lower # initialize inclusion indicator matrix
# AX: Check inconsistency with the ground truth
p <- nrow(dat_sim$beta)
if (p != nrow(fit$betaHat)) {
  stop(
    "Fitted design matrix has different structure than the ground truth.", "\n",
    "If the formula is consistent, this may be caused by a factor ",
    "having a different number of levels."
  )
}

# iterate across regression coefficients, calculate PIs and determine inclusion across fn. domain
for (r in 1:p) {
  # joint CIs
  lower_joint[r, ] <- fit$betaHat[r, ] -
    fit$qn[r] * sqrt(diag(fit$betaHat_var[, , r]))
  upper_joint[r, ] <- fit$betaHat[r, ] +
    fit$qn[r] * sqrt(diag(fit$betaHat_var[, , r]))

  # check whether estimate in CI
  joint_incl[r, ] <- dat_sim$beta[r, ] <= upper_joint[r, ] &
    dat_sim$beta[r, ] >= lower_joint[r, ]

  # naive CIs
  lower[r, ] = fit$betaHat[r, ] - 1.96*sqrt(diag(fit$betaHat_var[, , r]))
  upper[r, ] = fit$betaHat[r, ] + 1.96*sqrt(diag(fit$betaHat_var[, , r]))

  # check whether estimate in CI
  naive_incl[r, ] <- dat_sim$beta[r, ] <= upper[r, ] &
    dat_sim$beta[r, ] >= lower[r, ]
}

# calculate regression coefficient inclusion in 95% CI (across functional domain)
joint_incl <- rowSums(joint_incl) / ncol(joint_incl)
naive_incl <- rowSums(naive_incl) / ncol(naive_incl)

# average CI coverage across regression coefficients
res$beta_CI_joint <- mean(joint_incl)
res$beta_CI_naive <- mean(naive_incl)
res$fLME_beta_rmse <- sqrt(mean((dat_sim$beta - fit$betaHat)^2))

res$beta_CI_joint_noIntrcpt <- mean(joint_incl[-1])
res$beta_CI_joint0 <- mean(joint_incl[1])
res$beta_CI_naive_noIntrcpt <- mean(naive_incl[-1])
res$fLME_beta_rmse_noIntrcpt <- sqrt(
  mean((dat_sim$beta[-1, ] - fit$betaHat[-1, ])^2)
)


# 4 Save results and models ####################################################

# res$run_id <- run_id

# Write run details
# res$iter <- iter
# res$n_ids <- params_json$n_ids
# res$knots_div <- params_json$knots_div
# res$smooth_method <- params_json$smooth_method

write.csv(
  do.call(cbind, res),
  # paste0(save_dir, sprintf("%02d_%04d", run_id, iter), '.csv'),
  paste0(save_dir, sprintf("%04d", iter), '.csv'),
  row.names = F # col.names = F
)

if (params_json$save_coefs) {
  # Remove unnecessary and massive components
  fit$H <-
    fit$R <-
    fit$G <-
    fit$GHat <-
    fit$Z <- NULL

  saveRDS(fit, paste0(coefs_dir, "iter_", sprintf("%03d", iter), ".RDS"))
}

rm(fit)

# 5 pffr models ################################################################

if (params_json$run_pffr) {
  # Same analyses as fastFMM

  set.seed(iter)

  if (fixed_noise) {
    dat_sim <- synthesize_data(
      I = params_json$n_ids,
      J = params_json$mean_trials,
      rand_scale = params_json$rand_scale,
      sigma_eps = params_json$sigma_eps,
      asis_cols = T
    )
  } else {
    dat_sim <- synthesize_data(
      I = params_json$n_ids,
      J = params_json$mean_trials,
      snr_B = params_json$snr_B,
      snr_sigma = params_json$snr_sigma,
      asis_cols = T
    )
  }

  timeStart <- Sys.time()
  fit_pffr <- refund::pffr(
    formula = as.formula(params_json$refund_formula),
    algorithm = "bam",
    method = "GCV.Cp",
    family = "gaussian",
    bs.yindex = list(bs = "tp", k = nknots_min, m = c(2, 1)),
    data = dat_sim$data
  )
  timeEnd <- Sys.time()

  res$pffr_time <- as.numeric(difftime(timeEnd, timeStart, units = "mins"))
  message(
    "Finished fitting refund models. Time (min): ", signif(res$pffr_time, 3)
  )

  # 4.1 95% CIs ================================================================

  coef_pffr <- coef(fit_pffr, n1 = L)
  betaHat_pffr <- betaHat_pffr_var <- dat_sim$beta
  betaHat_pffr[1, ] <- as.vector(coef_pffr$smterms$`Intercept(yindex)`$value) +
    fit_pffr$coefficients[1]
  betaHat_pffr[2, ] <- as.vector(coef_pffr$smterms$`x(yindex)`$value)

  betaHat_pffr_var[1, ] <- as.vector(coef_pffr$smterms$`Intercept(yindex)`$se^2)
  betaHat_pffr_var[2, ] <- as.vector(coef_pffr$smterms$`x(yindex)`$se^2)

  # initialize CIs
  lower <- upper <- lower_joint <- upper_joint <- matrix(
    NA,
    nrow = nrow(dat_sim$beta),
    ncol = ncol(dat_sim$beta)
  )
  # initialize inclusion indicator matrix
  joint_incl <- naive_incl <- lower

  # iterate across regression coefficients
  # calculate PIs and determine inclusion across fn. domain
  for (r in 1:p) {
    # naive CIs
    lower[r, ] = betaHat_pffr[r, ] - 1.96 * sqrt(betaHat_pffr_var[r, ])
    upper[r, ] = betaHat_pffr[r, ] + 1.96 * sqrt(betaHat_pffr_var[r, ])

    # check whether estimate in CI
    naive_incl[r, ] <- dat_sim$beta[r, ] <= upper[r, ] &
      dat_sim$beta[r, ] >= lower[r, ]
  }

  # 4.2 pffr results ===========================================================

  # calculate regression coefficient inclusion in 95% CI (across functional domain)
  naive_incl <- rowSums(naive_incl) / ncol(naive_incl)

  # average CI coverage across regression coefficients
  res$pffr_CI_joint <- mean(naive_incl)
  # rmse of model fit
  res$pffr_beta_rmse <- sqrt(mean((dat_sim$beta - betaHat_pffr)^2))

  # Coverage of the intercept
  res$pffr_CI_joint0 <- mean(naive_incl[1])

  # don't evaluate performance on intercept
  res$pffr_CI_joint_noIntrcpt <- mean(naive_incl[-1])
  res$pffr_beta_rmse_noIntrcpt <- sqrt(
    mean((dat_sim$beta[-1, ] - betaHat_pffr[-1, ])^2)
  )
  res$pffr_bias_noIntrcpt <- mean(
    dat_sim$beta[beta_idx, ] - betaHat_pffr[beta_idx, ]
  )

  # 4.3 Save results and coefficients ##########################################

  write.csv(
    do.call(cbind, res),
    # paste0(save_dir, sprintf("%02d_%04d", run_id, iter), '.csv'),
    paste0(save_dir, sprintf("%04d", iter), '.csv'),
    row.names = F # col.names = F
  )

  if (params_json$save_pffr) {

    # Predict on dummy variables
    dummy <- data.frame(
      x = 0,
      id = params_json$n_ids + 1,
      yindex.vec = 1:L,
      type = "link",
      se.fit = T
    )
    preds <- suppressWarnings(
      predict(fit_pffr, dummy, type = "terms", se.fit = T)
    )

    pffr_res <- list(
      betaHat_0 = fit_pffr$coefficients[1],
      betaHat = betaHat_pffr,
      betaHat_var = betaHat_pffr_var,
      preds = as.numeric(preds$fit[[1]][1, ]),
      preds_se = as.numeric(preds$se.fit[[1]][1, ])
    )

    pffr_dir <- paste0(base_dir, "pffr/", params_id, "/")
    if (!dir.exists(pffr_dir)) dir.create(pffr_dir)
    saveRDS(pffr_res, paste0(pffr_dir, "iter_", sprintf("%03d", iter), ".RDS"))
  }
}
