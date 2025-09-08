# Functions for generating synthetic data for concurrent models

# Certain helper functions are not called during simulation
# This allow for flexibility in specifying the functions

# Basic parameters #############################################################

# Reminder: these are only set within these files
# This allows for some documentation of what settings were used for simulations

S <- seq(0, 5, by = 0.05)
L <- length(S)

# Functional covariates ########################################################

# Radial basis function kernel and i.i.d. multivariate Normal draws

rbfkernel <- function(x, y, sigma = 1, l = 10) {
  sigma^2 * exp(-(x - y)^2 / (2 * l^2))
}

fun_cov_sig <- matrix(0, nrow = L, ncol = L)

for (i in 1:L) {
  for (j in i:L) {
    fun_cov_sig[i, j] <- fun_cov_sig[j, i] <- rbfkernel(i, j)
  }
}

# Draws a single observation x_ij
fun_cov <- function(noisy = T) {
  x <- MASS::mvrnorm(n = 1, mu = rep(0, L), Sigma = fun_cov_sig)
  if (noisy)
    x <- x + rnorm(n = length(x), mean = 0, sd = 0.2)
  x
}

# Fixed coefficients ###########################################################

beta_0 <- function() {
  -.75 - 0.5 * sin(0.5 * pi * S) - 0.5 * cos(0.5 * pi * S)
}

beta_1 <- function() {
  1 * dnorm((S - 0.6) / 0.1) +
    2 * dnorm((S - 0.2) / 1) +
    1 * dnorm((S - 4) / 1.5)
}

# Random coefficients ##########################################################

psi_0 <- function() {
  1.5 * sin(0.5 * pi * S) - cos(0.5 * pi * S)
}

psi_1 <- function() {
  sin(pi * S)
}

# Make the outputs orthonormal
psi_mat <- rbind(
  psi_0() / sqrt(sum(psi_0()^2)),
  psi_1() / sqrt(sum(psi_1()^2))
)

# Output fixed coefficients as a matrix
gamma_0 <- function() {
  c_i <- MASS::mvrnorm(n = 1, mu = rep(0, 2), Sigma = diag(c(3, 1.5)))
  as.numeric(c_i %*% psi_mat)
}

gamma_1 <- function() {
  c_i <- MASS::mvrnorm(n = 1, mu = rep(0, 2), Sigma = diag(c(1, 2)))
  as.numeric(c_i %*% psi_mat)
}
