library(this.path)
setwd(this.path::this.dir())
source("../../utilities/runmefirst.R")
source("../../utilities/regressions.R")
library(latex2exp)
options(width = 120)

# number of jkp characteristics?
from_dir <- "../../tmp/raw_data/jkp_chars_not_lagged/1_unif/"
files <- list.files(from_dir, pattern = "RDS")

out <- data.table()
for (this_file in files) {
  tic(this_file)
  tmp <- readRDS(paste0(from_dir, this_file))
  out <- rbind(out, unique(tmp[, .(var)]))
  toc()
}

# a total of X variables
