# --- Compute stock-specific volatility using rolling windows (after removing market factor)
library(zoo)
library(this.path)
setwd(this.path::this.dir())
nc <- detectCores() - 2
source("../../utilities/runmefirst.R")
options(width = 120)

# load daily stock returns
tic("loading daily stock returns")
from_dir <- "../../../../../data/stocks/prices/daily_return/"
files <- list.files(from_dir, pattern = ".RDS")
p.read_one <- function(this_file) {
    data <- readRDS(paste0(from_dir, this_file))[, .(date, permno, ret)]
    return(data)
}

data <- rbindlist(mclapply(files, p.read_one, mc.cores = nc)) %>% na.omit()
rm(files, p.read_one, from_dir)
toc()

# merge with market factor return
tmp <- readRDS("../../../../../data/factors/daily_ff.RDS")[, .(date, rf, mktrf)]
data <- merge(data, tmp, by = "date")
data[, ret_rf := ret - rf]
data[, c("rf", "ret") := NULL]

# estimate by quarter
tmp <- unique(data[, .(date)]) %>%
    mutate(yyyymm = 100 * year(date) + 3 * quarter(date)) %>%
    setDT()
tmp[, idx := frank(yyyymm, ties.method = "dense")]
data <- merge(data, tmp, by = "date")
data_all <- copy(data)
rm(tmp, data)

# in each quarter, use 12q to estimate volatility
min_obs <- 252 * 3 / 2 # require observing data for half of the sample
lookback <- 12

# worker function
p.get_one_idx <- function(this_idx) {
    # this_idx <- target_idx[1]
    data <- data_all[idx %in% c(this_idx - c(0:(lookback - 1)))]

    # keep track of stocks with enough obs
    data[, obs := .N, permno]
    data <- data[obs >= min_obs][, obs := NULL]

    # regress out market factor
    data[, beta := cov(ret_rf, mktrf) / var(mktrf), permno]
    data[, ret_resid := ret_rf - beta * mktrf]

    # compute volatility
    data <- data[, .(idx = this_idx, vol_total = sd(ret_rf) * sqrt(252), vol_idio = sd(ret_resid) * sqrt(252)), permno]
    return(data)
}

# estimate
tmp <- unique(data_all[, .(idx)])[order(idx)]
tmp <- tmp[idx >= (min(idx) + lookback - 1)]

block_size <- 12
tmp[, block_idx := ceiling(.I / block_size)]

# on a 6 core machine, takes around one min
out <- data.table()
for (i in unique(tmp[, block_idx])) {
    tic(i / max(tmp[, block_idx]))
    out <- rbind(out, rbindlist(mclapply(tmp[block_idx == i, idx], p.get_one_idx, mc.cores = nc)))
    gc()
    toc()
}

# get yyyymm
tmp <- unique(data_all[, .(idx, yyyymm)])
out <- merge(out, tmp, by = "idx")[, idx := NULL]
setcolorder(out, c("yyyymm", "permno", "vol_total", "vol_idio"))

to_file <- "tmp/vols/stock_vol_rolling_12q.RDS"
dir.create(dirname(to_file), showWarnings = FALSE, recursive = TRUE)
saveRDS(out, to_file)

# # --- sanity check: compare with earlier

# # computed earlier
# data <- readRDS("../../../../../../../../Desktop/J-Leaves/data/stockprices/processed/ivol/byMonth_ff3_usingDailyData.RDS")
# data[, obs := NULL]
# data <- data[yyyymm >= 199003]

# # rolling 12q
# data[, total_vol := frollapply(total_vol, 36, mean) * sqrt(252), permno]
# data[, ivol := frollapply(ivol, 36, mean) * sqrt(252), permno]
# names(data)[3] <- "tvol"
# data <- melt(data, id.vars = c("yyyymm", "permno"), variable.name = "type", value.name = "old") %>%
#     mutate(type = as.character(type)) %>%
#     na.omit() %>%
#     setDT()

# # get new
# tmp <- readRDS("tmp/vols/stock_vol_rolling_12q.RDS")
# names(tmp)[3:4] <- c("tvol", "ivol")
# tmp <- melt(tmp, id.vars = c("yyyymm", "permno"), variable.name = "type", value.name = "new") %>%
#     mutate(type = as.character(type)) %>%
#     na.omit() %>%
#     setDT()

# data <- merge(data, tmp, by = c("yyyymm", "permno", "type"))
# rm(tmp)

# dcast(data[, cor(old, new), .(yyyymm, type)], yyyymm ~ type, value.var = "V1")
