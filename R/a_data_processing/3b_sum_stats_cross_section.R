# === produce summary statistics as time-series of cross-sectional moments
library(this.path)
setwd(this.path::this.dir())
source("../utilities/runmefirst.R")

data <- readRDS("../tmp/raw_data/reg_inputs/all_ofi_and_ret.RDS")
data <- data[yyyymm >= 199306]
data[, ofi := 100 * ofi]
data[, ret := 100 * ret]

# convert to wide format
data <- dcast(data, yyyymm + permno + ret ~ type, value.var = "ofi", fun.aggregate = mean)

# get lagged market cap
tmp <- readRDS("../../../data/stocks/prices/quarterly_return.RDS")
tmp <- tmp[, .(yyyymm, permno, me_1)]
data <- merge(data, tmp, by = c("yyyymm", "permno"), all.x = T)

# into long format
data <- data.table(melt(data, id.vars = c("yyyymm", "permno")))
names(data)[3:4] <- c("var", "xx")
data <- data[!is.na(xx)]

# compute sum-stats for each (yyyymm, var)
p.get_one_period <- function(this_data) {
  # this_data <- data_list[[1]]
  tmp <- this_data[, list(
    var = this_data[1, var],
    yyyymm = this_data[1, yyyymm],
    obs = nrow(this_data),
    Mean = mean(xx),
    StDev = sd(xx),
    p01 = quantile(xx, .01),
    p05 = quantile(xx, .05),
    p25 = quantile(xx, .25),
    p50 = quantile(xx, .50),
    p75 = quantile(xx, .75),
    p95 = quantile(xx, .95),
    p99 = quantile(xx, .99)
  )]
  return(tmp)
}

data_list <- split(data, by = c("yyyymm", "var"))
data_list <- data_list[sapply(data_list, nrow) > 0]
out <- rbindlist(lapply(data_list, p.get_one_period))

# get time-series averages
tmp <- out[, .(obs = mean(obs)), var]
out[, obs := NULL]
out <- melt(out, id.vars = c("var", "yyyymm"))
out <- out[, .(value = mean(value)), .(var, variable)]
out <- dcast(out, var ~ variable, value.var = "value")
out <- merge(tmp, out, by = "var")
rm(tmp)

to_dir <- "../tmp/raw_data/sum_stats/cross_section/"
dir.create(to_dir, recursive = T, showWarnings = F)
saveRDS(out, paste0(to_dir, "all.RDS"))

# # === SANITY CHECK: compare with old results

# # bmi
# new = readRDS('tmp/raw_data/sum_stats/bmi.RDS')
# old = readRDS('../20250117_quarterly/tmp/raw_data/sum_stats/bmi.RDS')
# mean(new == old, na.rm = T)

# # fit and ofi
# new = readRDS('tmp/raw_data/sum_stats/fit_and_ofi.RDS')[, var := as.character(var)]
# old = readRDS('../20250117_quarterly/tmp/raw_data/sum_stats/fit_and_ofi.RDS')[, var := as.character(var)]
# old = old[var != 'OFI']
# old[var == 'OFI_resid', var := 'OFI']

# new = melt(new, id.vars = 'var')
# old = melt(old, id.vars = 'var')
# compare = merge(new, old, by = c('var','variable'), all = T); rm(new, old)
# compare[, cor(value.x, value.y)]
