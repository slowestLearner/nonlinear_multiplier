# --- Let's check the cross-sectional correlation of dispersion measures
library(this.path)
setwd(this.path::this.dir())
source("../../utilities/runmefirst.R")
options(width = 200)

data <- readRDS("tmp/raw_files/ofi_dispersion.RDS")[, c("obs", "ofi_sum", "abs_ofi_sum") := NULL]

data_list <- split(data, by = "yyyymm")

p.get_one_period <- function(this_data) {
    # this_data <- data_list[[1]]
    return(cor(this_data[, -c(1:2), with = F], method = "spearman"))
}

p.get_one_period(data_list[[100]])
