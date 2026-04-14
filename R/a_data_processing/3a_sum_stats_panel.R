# === script to produce summary statistics for demand measures
library(this.path)
setwd(this.path::this.dir())
source("../utilities/runmefirst.R")

# ======= 1) BMI

# read data
data <- readRDS("../tmp/raw_data/reg_inputs/all_ofi_and_ret.RDS")
data <- data[type == "BMI"][, type := NULL]
setnames(data, "ofi", "d_bmi")
data[, d_bmi := 100 * d_bmi] # convert to percentage

# compute summary statistics
tmp <- data.table(obs = nrow(data) / data[, length(unique(yyyymm))]) # num of observations per year
tmp <- cbind(tmp, data[, list(
  Mean = mean(d_bmi), StDev = sd(d_bmi),
  p01 = quantile(d_bmi, .01),
  p05 = quantile(d_bmi, .05),
  p25 = quantile(d_bmi, .25),
  p50 = quantile(d_bmi, .50),
  p75 = quantile(d_bmi, .75),
  p95 = quantile(d_bmi, .95),
  p99 = quantile(d_bmi, .99)
)])
to_dir <- "../tmp/raw_data/sum_stats/panel/"
dir.create(to_dir, recursive = T, showWarnings = F)
saveRDS(tmp, paste0(to_dir, "bmi.RDS"))

# ======= 2) OFI and FIT, as well as other variables
rm(list = ls())

data <- readRDS("../tmp/raw_data/reg_inputs/all_ofi_and_ret.RDS")[type != "BMI"]
data <- data[yyyymm >= 199306]
data[, ofi := 100 * ofi]
data[, ret := 100 * ret]

# data_ret = data[, .(ret = last(ret)), .(yyyymm, permno)]
# data[, ret := NULL]
# data = dcast(data, yyyymm + permno ~ type, value.var = 'ofi')
# data = merge(data, data_ret, by = c('yyyymm','permno')); rm(data_ret)

# convert to wide format
data <- dcast(data, yyyymm + permno + ret ~ type, value.var = "ofi")

# get lagged market cap
tmp <- readRDS("../../../data/stocks/prices/quarterly_return.RDS")
tmp <- tmp[, .(yyyymm, permno, me_1)]
data <- merge(data, tmp, by = c("yyyymm", "permno"))
rm(tmp)

# into long format
data <- data.table(melt(data, id.vars = c("yyyymm", "permno")))
names(data)[3:4] <- c("var", "xx")

# helper to compute summary statistics for one variable
p.getOne <- function(this_v) {
  tmp <- copy(data[var == this_v])
  tmp <- tmp[!is.na(xx)]
  return(tmp[, list(
    var = this_v,
    obs = round(tmp[, length(permno), yyyymm][, mean(V1)]),
    x_mean = mean(xx), x_sd = sd(xx),
    p01 = quantile(xx, .01),
    p05 = quantile(xx, .05),
    p25 = quantile(xx, .25),
    p50 = quantile(xx, .50),
    p75 = quantile(xx, .75),
    p95 = quantile(xx, .95),
    p99 = quantile(xx, .99)
  )])
}

all_vars <- unique(data[, var])
data <- Reduce(rbind, lapply(all_vars, p.getOne))
rm(all_vars, p.getOne)

options(width = 200)

to_dir <- "../tmp/raw_data/sum_stats/panel/"
dir.create(to_dir, recursive = T, showWarnings = F)
saveRDS(data, paste0(to_dir, "fit_and_ofi.RDS"))

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
