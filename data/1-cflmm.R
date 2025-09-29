# Data-generating functions
source("fun/fig1a.R")
# Function to create random draws from data-generating functions
source("fun/synthesize_data.R")

# Synthesize the data ##########################################################

# S comes from the synthetic data function
L <- length(S)
s_pt <- S[25]

set.seed(20250605)
dat_list <- synthesize_data(
  I = 40,
  J = 50,
  snr_B = 0.5,
  snr_sigma = 0.5,
  ptwise_snr = T,
  y_name = "y",
  x_name = "x",
  return_fixed = T,
  return_noise = T,
  asis_cols = F,
  fix_J = T
)

dat <- dat_list$data
# Write the full synthetic data frame
write.csv(dat, file = "data/1-cflmm.csv", row.names = F)


