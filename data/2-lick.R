source("fun/jeong_3.R")
source("fun/time_align.R")

pp <- jsonlite::read_json(
  paste0(params_dir, "jeong_3-reward-", pre_rwd_len, ".json"),
  simplifyVector = T
)

class(pp) <- pp$align_to

ids <- c(
  'HJ-FP-M2',
  'HJ-FP-M3',
  'HJ-FP-M4',
  'HJ-FP-F1',
  'HJ-FP-F2',
  'HJ-FP-M6',
  'HJ-FP-M7',
  'HJ-FP-M8'
)

session_min <- pp$session_min
session_max <- pp$session_max

iri_cutoff <- 5

silent <- T

# I won't save these in the github repo since it seems inappropriate to
# reupload an existing dataset wholesale

# They can be found at https://dandiarchive.org/dandiset/000351/draft
reward_list <- mapply(
  function(id, s) {
    f_name <- paste0(raw_dir, "sub-", id, "_ses-Day", s, ".nwb")
    # Check for file existence and return early if not found
    if (!file.exists(f_name)) return()

    # Read and align ttls and photometry data
    dat <- read_matlab(f_name)

    # Stop function if IRI is too low
    iri_length <- mean(
      diff(dat$ttls$time[dat$ttls$event == pp$label_rwd])
    )

    if (!silent)
      message(id, s, " mean(IRI): ", round(mean(iri_length), 2), "s.")
    if (iri_length < iri_cutoff)
      return()

    # Filter trials
    keepers <- filter_trials(pp, dat)
    if (sum(keepers) < 1)
      return()

    # Align trials and label the photometry data
    trial_idx <- index_trials(pp, dat, keepers, downsample = T)
    photo_df <- apply(trial_idx, 1, function(x) dat$dat_photo$data[x]) %>%
      as.data.frame()
    colnames(photo_df) <- paste0("photometry_", 1:ncol(photo_df))

    # Get the functional covariate of lick
    downsample_by <- round(pp$Hz / pp$target_Hz)
    lick_df <- get_lick_fun(pp, dat, trial_idx, downsample_by)
    colnames(lick_df) <- paste0("lick_", 1:ncol(photo_df))

    # Get IRI
    iri <- c(
      dat$ttls$time[dat$ttls$event == pp$label_rwd][1],
      diff(dat$ttls$time[dat$ttls$event == pp$label_rwd])
    )
    # Get lick times, totals, and probability
    sess_info <- data.frame(
      id = id,
      session = s,
      trial = (1:length(keepers))[keepers],
      lick_time = get_lick_time(pp, dat)[keepers],
      lick_prob_050 = get_lick_prob(pp, dat, T, 0.5)[keepers],
      lick_prob_100 = get_lick_prob(pp, dat, T, 1.0)[keepers],
      lick_prob_150 = get_lick_prob(pp, dat, T, 1.5)[keepers],
      lick_prob_200 = get_lick_prob(pp, dat, T, 2.0)[keepers],
      # lick_I_050 = get_lick_I(pp, dat, 0.5)[keepers],
      # lick_I_100 = get_lick_I(pp, dat, 1.0)[keepers],
      # lick_I_150 = get_lick_I(pp, dat, 1.5)[keepers],
      # lick_I_200 = get_lick_I(pp, dat, 2.0)[keepers],
      licks = get_lick_total(pp, dat)[keepers],
      iri = iri[keepers]
    )

    return(cbind(sess_info, photo_df, lick_df))
  },
  rep(ids, session_max - session_min + 1),
  rep(session_min:session_max, each = length(ids))
)

reward <- do.call(rbind, reward_list)
rm(reward_list)
