# --- produce summary statistics as time-series of cross-sectional moments
library(this.path)
setwd(this.path::this.dir())
source("../utilities/runmefirst.R")

data <- readRDS("../tmp/raw_data/reg_inputs/all_ofi_and_ret.RDS")[yyyymm %in% 199306:202212]


# turn into percentages
data[, ofi := 100 * ofi]
data[, ret := 100 * ret]

# convert to wide format
data <- dcast(data, yyyymm + permno + ret ~ type, value.var = "ofi", fun.aggregate = mean)

# merge with market cap
tmp <- readRDS("../../../../data/stocks/prices/quarterly_return.RDS")[, .(yyyymm, permno, me_1)]
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

to_file <- "../tmp/raw_data/sum_stats/all.RDS"
dir.create(dirname(to_file), showWarnings = FALSE, recursive = TRUE)
saveRDS(out, to_file)

# # === SANITY CHECK: sufficiently similar to old results
# new <- readRDS("../tmp/raw_data/sum_stats/all.RDS")
# old <- readRDS("../tmp/raw_data/sum_stats_todel/cross_section/all.RDS")
