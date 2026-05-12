# --- industry dummy controls
library(this.path)
setwd(this.path::this.dir())
source("../utilities/runmefirst.R")

# sic codes
data <- readRDS("../../../../data/stocks/prices/quarterly_return.RDS")[, .(yyyymm, permno, siccd)]

# go through different industry defs
from_dir <- "../../../../data/stocks/controls/industry_definitions/"
ind_nums <- c(5, 12, 17, 30, 49)
for (ind_num in ind_nums) {
    # ind_num <- 5
    print(paste0("Processing industry def for ", ind_num, " industries"))
    tic()
    tmp <- readRDS(paste0(from_dir, ind_num, "ind.RDS"))
    for (i in 1:nrow(tmp)) {
        data[siccd %in% c(tmp[i, sicStart]:tmp[i, sicEnd]), paste0("ind", ind_num) := tmp[i, ind]]
    }
    toc()
}

# fill in other
data[is.na(data)] <- "other"
data[, siccd := NULL]

# save locally
to_file <- "../../../../data/stocks/controls/ff_refined_industry_dummies.RDS"
dir.create(dirname(to_file), showWarnings = FALSE, recursive = TRUE)
saveRDS(data, to_file)
