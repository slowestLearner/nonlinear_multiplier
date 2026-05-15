library(this.path)
setwd(this.path::this.dir())
source("../../utilities/runmefirst.R")
source("../../utilities/regressions.R")
options(width = 200)

# number of jkp characteristics?
from_dir <- "../../tmp/raw_data/jkp_chars_not_lagged/1_unif/"
files <- list.files(from_dir, pattern = "RDS")

out <- data.table()
for (this_file in files) {
  tic(this_file)
  tmp <- readRDS(paste0(from_dir, this_file))
  out <- rbind(out, unique(tmp[, .(file = this_file, var)]))
  toc()
}

# save locally
out[, file_num := as.integer(gsub("part_|.RDS", "", file))]
out <- out[order(file_num, var)]

to_file <- "tmp/raw_files/list_of_jkp_controls.RDS"
dir.create(dirname(to_file), recursive = T, showWarnings = F)
saveRDS(out, to_file)
