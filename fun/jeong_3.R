# Reading files ################################################################

#' Read MATLAB files
#' 
#' Read the MATLAB output from experiment 3 and do basic data frame formatting.
#'
#' @param f_name name of the MATLAB file to read
#' @output A list with two data frames: 
#' 1) time and event, 
#' 2) time and photometry.

read_matlab <- function(f_name) {
  dat <- rhdf5::h5read(file = f_name, name = "/acquisition")
  
  # event time ttls
  ttls <- as.data.frame(
    cbind(
      dat$eventlog$eventtime, 
      dat$eventlog$eventindex
    )
  )
  colnames(ttls) <- c("time", "event")
  
  dat_photo <- rhdf5::h5read(
    file = f_name, 
    name = "/processing/photometry/dff"
  )
  
  dat_photo <- as.data.frame(do.call(cbind, dat_photo))
  dat_photo <- dat_photo[!is.na(dat_photo$data), ] 
  
  ttls <- time_align(
    time_truth = dat_photo$timestamps, 
    data = ttls, 
    name = "time"
  )
  
  return(
    list(dat_photo = dat_photo, ttls = ttls)
  )
}

# Filter trials ################################################################

#' Filter trials based on a target alignment parameter.
#'
#' @param pp  preprocessing object with parameters
#' @param ... additional arguments

filter_trials <- function(pp, ...) {
  UseMethod("filter_trials")
}

#' Filter trials based on length
#' 
#' Trials are extracted based on the first lick after a reward. Trials shorter
#' than the preset minimum trial length are removed.
#' 
#' @param dat list of data frames output from `read_matlab`

filter_trials.first_lick <- function(pp, dat) {
  ttls <- dat$ttls
  dat_photo <- dat$dat_photo
  
  reward_times <- ttls$time[which(ttls$event == pp$label_rwd)]
  keep_trials <- rep(T, length(reward_times))
  
  # Remove trials that are shorter than min_tm
  short_trials <- unique(which(diff(reward_times) < pp$min_tm) + 1)
  keep_trials[short_trials] <- F
  
  # Remove trials that start too early for sufficient photometry data
  pre_min_tm <- pp$pre_reward_length + pp$pre_lick_reward_period 
  photo_start <- min(dat_photo$timestamps)
  early_trials <- I(reward_times - photo_start <= pre_min_tm)
  keep_trials[early_trials] <- F
  
  # Remove trials that don't have a lick after the reward
  lick_times <- ttls$time[which(ttls$event == pp$label_lick)]
  lick_totals <- sapply(
    reward_times, 
    function(x)  
      sum(lick_times > x & lick_times <= (x + pp$min_tm))
  ) 
  keep_trials <- I(keep_trials * I(lick_totals > 0) == 1) 
  
  keep_trials
}

filter_trials.lick <- function(pp, dat) {
  ttls <- dat$ttls
  dat_photo <- dat$dat_photo
  
  lick_times <- ttls$time[which(ttls$event == pp$label_lick)]
  keep_trials <- rep(T, length(lick_times))
  
  # Remove trials that are shorter than min_tm
  short_trials <- unique(which(diff(lick_times) < pp$min_tm) + 1)
  keep_trials[short_trials] <- F
  
  # Remove trials that start too early for sufficient photometry data
  pre_min_tm <- pp$pre_reward_length + pp$pre_lick_reward_period 
  photo_start <- min(dat_photo$timestamps)
  early_trials <- I(lick_times - photo_start <= pre_min_tm)
  keep_trials[early_trials] <- F
  
  keep_trials  
}

filter_trials.reward <- function(pp, dat) {
  ttls <- dat$ttls
  dat_photo <- dat$dat_photo
  
  reward_times <- ttls$time[which(ttls$event == pp$label_rwd)]
  keep_trials <- rep(T, length(reward_times))
  
  # Remove trials that are shorter than min_tm
  short_trials <- unique(which(diff(reward_times) < pp$min_tm) + 1)
  keep_trials[short_trials] <- F
  
  # Remove trials that start too early for sufficient photometry data
  pre_min_tm <- pp$pre_reward_length + pp$pre_lick_reward_period 
  photo_start <- min(dat_photo$timestamps)
  early_trials <- I(reward_times - photo_start <= pre_min_tm)
  keep_trials[early_trials] <- F
  
  keep_trials
}

# Reindex trials ###############################################################

#' Downsample indices
#' 
#' @param pre_samps integer number of samples before the trial aligner
#' @param L integer size of the functional domain
#' @param downsample_by integer of how many entries to skip (-1)

downsample_indices <- function(pre_samps, L, downsample_by) {
  pre_idx <- sort(-seq(-pre_samps, -1, by = downsample_by))
  post_idx <- seq(
    pre_samps + downsample_by, L, by = downsample_by
  )
  unique(c(pre_idx, post_idx)) 
}

#' Select indices of trial events and photometry 
#' 
#' @param pp  preprocessing object with parameters
#' @param ... additional arguments

index_trials <- function(pp, ...) {
  UseMethod("index_trials")
}

#' @param pp preprocessing object with parameters
#' @param dat list of `ttls` and `dat_photo` data frames
#' @param keepers boolean list of which trials to keep in analysis
#' @param downsample boolean of whether to downsample observations
#' 
#' @out matrix of indices with each row corresponding to a trial

index_trials.first_lick <- function(pp, dat, keepers, downsample = T, ...) {
  # Setup variables, filtering by keepers
  reward_times <- dat$ttls$time[which(dat$ttls$event == pp$label_rwd)][keepers]
  lick_times   <- dat$ttls$time[which(dat$ttls$event == pp$label_lick)]
  
  pre_min_tm <- pp$pre_reward_length + pp$pre_lick_reward_period
  pre_samps <- pre_min_tm * pp$Hz
  
  post_min_tm <- pp$post_reward_length + pp$post_lick_reward_period
  post_samps <- post_min_tm * pp$Hz
  
  # Get trial indices
  trial_idx <- sapply(
    reward_times, 
    function(x) {
      # Get the first lick index
      first_lick_when <- lick_times[
        I(lick_times > x & lick_times <= (x + pp$min_tm))
      ][1]
      first_lick_idx <- which.min(
        abs(first_lick_when - dat$dat_photo$timestamps)
      )
      
      seq(
        first_lick_idx - pre_samps, 
        first_lick_idx + post_samps
      )
    }
  )
  
  # Return as-is if no downsampling required
  if (!downsample)
    return(trial_idx)
  
  downsample_by <- round(pp$Hz / pp$target_Hz)
  downsample_idx <- downsample_indices(
    pre_samps, nrow(trial_idx), downsample_by
  )
  
  # Return downsampled trials
  trial_idx[downsample_idx, ]
}

index_trials.reward <- function(pp, dat, keepers, downsample = T, ...) {
  # Setup variables, filtering by keepers
  reward_times <- dat$ttls$time[which(dat$ttls$event == pp$label_rwd)][keepers]
  
  pre_min_tm <- pp$pre_reward_length + pp$pre_lick_reward_period
  pre_samps <- pre_min_tm * pp$Hz
  
  post_min_tm <- pp$post_reward_length + pp$post_lick_reward_period
  post_samps <- post_min_tm * pp$Hz
  
  # Get trial indices
  trial_idx <- sapply(
    reward_times, 
    function(x) {
      reward_idx <- which.min(abs(x - dat$dat_photo$timestamps))
      seq(
        reward_idx - pre_samps, 
        reward_idx + post_samps
      )
    }
  )
  
  # Return as-is if no downsampling required
  if (!downsample)
    return(trial_idx)
  
  downsample_by <- round(pp$Hz / pp$target_Hz)
  downsample_idx <- downsample_indices(
    pre_samps, nrow(trial_idx), downsample_by
  )
  
  # Return downsampled trials
  trial_idx[downsample_idx, ]
}

# Generate covariates ##########################################################

stepwise_floor <- function(x, steps) {
  steps_ <- c(0, steps)
  max(steps_[steps_ < x])
}

# PDF of an exponential

get_reward_dexp <- function(pp, dat) {
  # Get pdf of exponential 
  reward_times <- dat$ttls$time[which(dat$ttls$event == pp$label_rwd)]
  reward_dexp <- dexp(
    x = dat$dat_photo$timestamps - 
      sapply(
        dat$dat_photo$timestamps, 
        stepwise_floor, 
        steps = reward_times
      ),
    rate = 1 / 12
  )
}

# Total number of licks in reward period

get_lick_total <- function(pp, dat) {
  reward_times <- dat$ttls$time[which(dat$ttls$event == pp$label_rwd)]
  lick_times <- dat$ttls$time[which(dat$ttls$event == pp$label_lick)] 
  
  sapply(
    reward_times, 
    function(x)  
      sum(lick_times > x & lick_times <= (x + pp$min_tm))
  ) 
}

get_lick_time <- function(pp, dat) {
  reward_times <- dat$ttls$time[which(dat$ttls$event == pp$label_rwd)]
  lick_times <- dat$ttls$time[which(dat$ttls$event == pp$label_lick)] 
  
  sapply(
    reward_times, 
    function(x)  
      lick_times[I(lick_times > x & lick_times <= (x + pp$min_tm) )][1] - x 
  )  
}

# Averaged number of licks in reward period

get_lick_prob <- function(
  pp, dat, downsample = T, reward_period_length = NULL
) {
  reward_times <- dat$ttls$time[which(dat$ttls$event == pp$label_rwd)]
  lick_times <- dat$ttls$time[which(dat$ttls$event == pp$label_lick)] 
  if (is.null(reward_period_length)) 
    reward_period_length <- pp$reward_period_length
  
  sapply(
    reward_times, 
    function(x) {
      lick_tot <- sum(
        lick_times > x & lick_times <= (x + reward_period_length)
      )
      
      if (downsample) {
        lick_tot / round(pp$target_Hz * reward_period_length)
      } else {
        lick_tot / round(pp$Hz * reward_period_length)
      }
    }
  )
}

# Indicator of licks in reward period

get_lick_I <- function(
    pp, dat, reward_period_length = NULL
) {
  reward_times <- dat$ttls$time[which(dat$ttls$event == pp$label_rwd)]
  lick_times <- dat$ttls$time[which(dat$ttls$event == pp$label_lick)] 
  if (is.null(reward_period_length)) 
    reward_period_length <- pp$reward_period_length
  
  sapply(
    reward_times, 
    function(x) {
      lick_tot <- sum(
        lick_times > x & lick_times <= (x + reward_period_length)
      )
      as.numeric(lick_tot > 0)
    }
  )
}

# Get lick as a functional covariate
get_lick_fun <- function(
  pp, 
  dat, 
  trial_idx, 
  downsample = T
) {
  # Get the lick observations
  lick_times <- dat$ttls$time[which(dat$ttls$event == pp$label_lick)] 
  lick_idx <- sapply(
    lick_times, 
    function(x) 
      which.min(abs(x - dat$dat_photo$timestamps))
  )
  # Creating a new vector of the lick observations
  lick_obs <- rep(0, length(dat$dat_photo$timestamps))
  lick_obs[lick_idx] <- 1
  
  if (!downsample) {
    as.data.frame( 
      apply(trial_idx, 1, function(x) lick_obs[x])
    )
  }
  
  # Return "1" if a lick occurs in the region of the downsample
  downsample_by <- round(pp$Hz / pp$target_Hz)
  apply(
    trial_idx,
    1, 
    function(row_idx) {
      sapply(
        row_idx, 
        function(x) 
          as.numeric(sum(lick_obs[(x - downsample_by + 1):x]) > 0)
      )
    }
  ) %>%
    as.data.frame() 
}