concat_results <- function(
    id,
    read_dir = "raw/outs/",
    niter = 1:200,
    silent = F,
    target_cols = 22
) {
  rpath <- paste0(read_dir, id, "/")

  do.call(
    rbind,
    lapply(
      paste0(rpath, sprintf("%04d", niter), ".csv"),
      function(f) {
        if (file.exists(f)) {
          temp <- read.csv(f)
          if (ncol(temp) != target_cols) {
            message(f, " has ", ncol(temp), " columns.")
            return()
          } else {
            temp
          }
        } else {
          if (!silent) message(f, " not found.")
        }
      }
    )
  )
}
