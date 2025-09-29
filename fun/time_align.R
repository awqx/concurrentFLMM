# written by Gabriel Loewinger


# time aligning
time_align <- function(time_truth, data, name = "timestamps", save_time = FALSE){
  # time_truth - is a vector of timestamps for the (usually photometry) timepoints to align to (ground truth)
  # data - is the dataset with timestamps to change to align to the time_truth
  # "name"  -  is the column name in "data" variable
  # save_time  - whether to save original time variable
  
  data_new <- data.table::as.data.table( data[name] ) 
  tm <- data.table::data.table(time_temp = as.numeric(time_truth) ) # vector of times we want to align with (photometry times) -- ground truth time
  tm[, time_aligned := time_truth]
  
  data.table::setkeyv(data_new, name) # column name in data file that we want to align to be consistent with photometry (the ground truth time)
  data.table::setkeyv(tm, c('time_temp'))
  
  data_new <- as.data.frame( tm[data, roll='nearest'] )
  
  # delete original time variable
  if(!save_time){
    data_new <- subset(data_new, select = -time_temp )
    colnames(data_new)[ colnames(data_new) == "time_aligned" ] <- name
  }   
  
  return(data_new)
  
}