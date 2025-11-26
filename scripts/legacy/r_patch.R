# Tame vroom on macOS/ARM (optional but helps stability)
Sys.setenv(VROOM_NUM_THREADS = "1", VROOM_SHOW_PROGRESS = "false")

# Safe base reader for the tiny progress CSV
read_base_progress <- function(pf) {
  if (file.exists(pf) && file.info(pf)$size > 0) {
    utils::read.csv(pf, stringsAsFactors = FALSE)
  } else {
    data.frame(stage=character(), pid=character(), status=character(), t=as.POSIXct(character()))
  }
}
