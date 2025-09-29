df <- read.csv("data-raw/3-d2pvt/D2PVT_FLMM_input_HSMSvsHH2O.csv") %>%
  select(-Trial_no) # deprecated trial number coding

# Remove "X."
times <- as.numeric(gsub("^X(\\.)*", "", colnames(df)[7:ncol(df)]))
# Convert to negative
times[1:24] <- -times[1:24]

colnames(df) <- c(
  "id",
  "session",
  "outcome",
  "latency",
  "trial",
  paste0("photometry_", 1:(ncol(df) - 5))
)

# Convert session from character to numeric
df$session <- as.numeric(gsub("^S", "", df$session))

# Clean up ID by removing the date and L#
# There should only be 5 different mice
df$id <- gsub("^[[:digit:]]{8}\\_", "", df$id)
df$id <- gsub("\\_L[[:digit:]]+$", "", df$id)

# Make outcome a factor
df$outcome <- as.factor(df$outcome)
# Scale trial
# df$trial <- as.numeric(scale(df$trial))

latencies <- df$latency
reward <- do.call(
  rbind,
  lapply(
    latencies,
    function(x)
      as.numeric(times > x)
  )
)
reward <- data.frame(reward)
colnames(reward) <- paste0("reward_", 1:length(times))

# Return data frame of photometry and functional covariate `approaching`
dat <- cbind(df, reward)

# Center and scale
dat <- dat %>%
  mutate(
    session = session - 1,
    trial = trial - 1
  ) %>%
  mutate(
    trial = trial / max(trial),
    session = session / max(session)
  )

# Remove columns of all zeros
zero_var <- sapply(
  select(dat, reward_1:reward_120),
  function(x)
    var(x) <= 0.01
)

zero_var_cols <- c(
  paste0("photometry_", which(zero_var)),
  paste0("reward_", which(zero_var))
)

dat <- dat[, !colnames(dat) %in% zero_var_cols]

reward_cols <- grep("^reward", colnames(dat))
colnames(dat) <- c(
  "id",
  "session",
  "outcome",
  "latency",
  "trial",
  paste0("photometry_", 1:length(reward_cols)),
  paste0("reward_", 1:length(reward_cols))
)


write.csv(dat, "data/3-d2pvt.csv", row.names = F)
