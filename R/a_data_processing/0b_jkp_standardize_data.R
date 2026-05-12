# --- standardize into unif[-0.5, 0.5]. Still have to it by file.
library(this.path)
setwd(this.path::this.dir())
source("../utilities/runmefirst.R")

# read in data and put into long format
from_dir <- "../tmp/raw_data/jkp_chars_not_lagged/0_raw/"
files <- list.files(from_dir, pattern = "*.RDS")

to_dir <- "../tmp/raw_data/jkp_chars_not_lagged/1_unif/"
dir.create(to_dir, recursive = T, showWarnings = F)

# restrict attention to the stocks we care about
stock_data <- readRDS("../../../../data/stocks/prices/quarterly_return.RDS")[, .(yyyymm, permno)] %>% unique()
stock_data <- stock_data[yyyymm %in% 199212:202303] # include 1 quarter before and after the sample

# process one file
p.get_one <- function(this_file) {
    # this_file <- files[1]
    data <- readRDS(paste0(from_dir, this_file))
    data <- merge(data, stock_data, by = c("yyyymm", "permno"), all.y = T)

    # first turn into numeric
    cols_to_melt <- names(data)[!names(data) %in% c("yyyymm", "permno")]
    data[, (cols_to_melt) := lapply(.SD, as.numeric), .SDcols = cols_to_melt]

    # into long table
    data <- melt(data, id.vars = c("yyyymm", "permno"), variable.name = "var", value.name = "char", variable.factor = FALSE) %>% setDT()

    # check availability. Require having at least 20% in all periods
    tmp <- data[, .(frac_obs = sum(!is.na(char)) / .N), .(yyyymm, var)]
    tmp[, min_frac_obs := min(frac_obs), var]
    tmp <- tmp[min_frac_obs > .2]
    tmp <- tmp[, .(yyyymm, var)]
    data <- merge(data, tmp, by = c("yyyymm", "var"))

    # turn into unif[-0.5, 0.5], fill zeros
    data[!is.na(char), char := rank(char, ties.method = "first") / sum(!is.na(char)) - .5, .(yyyymm, var)]
    data[is.na(char), char := 0]

    to_file <- paste0(to_dir, this_file)
    saveRDS(data, to_file)
    rm(data)
    gc()
}

# process all
tic()
mclapply(files, p.get_one, mc.cores = 2)
toc()
