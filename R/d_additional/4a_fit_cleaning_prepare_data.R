# code adapted from /Users/slowlearner/Dropbox/SpeculativeIdeas/Nonlinear_multipliers/formal_tests/20250117_quarterly/5_demand_cleaning.Rmd
# Put together a regression table.
library(this.path)
setwd(this.path::this.dir())
source("../utilities/runmefirst.R")

# get all trades in S12
files <- list.files("../../../../data/institutional/s12_quarterly_trades/", full.names = T)
data <- rbindlist(mclapply(files, readRDS, mc.cores = nc))
rm(files)
gc()

# clean data slightly
setnames(data, c("prc_tfn_1", "shrout_tfn_1"), c("prc_1", "shrout_1"))
data[, ncusip := NULL]
data <- data[!is.na(permno)]
data <- data[!is.na(prc_1) & !is.na(shrout_1)]

# restrict universe to common stocks
tmp <- readRDS("../../../../data/stocks/prices/quarterly_return.RDS")[, .(yyyymm, permno)]
data <- merge(data, tmp, by = c("yyyymm", "permno"))
rm(tmp)
gc()

# require the trades to be between quarter ends
tmp <- unique(data[, .(yyyymm, wficn, rdate_1, rdate)])
tmp <- merge(tmp, tmp[, .(rdate_1_max = max(rdate_1), rdate_max = max(rdate)), yyyymm], by = "yyyymm")
tmp <- tmp[(rdate_1 == rdate_1_max) & (rdate == rdate_max), .(yyyymm, wficn)]
data <- merge(data, tmp, by = c("yyyymm", "wficn"))
rm(tmp)
data[, c("rdate", "rdate_1") := NULL]

# deal with the rare non-unique entries
data[, obs := .N, .(yyyymm, wficn, permno)]
data <- data[obs == 1][, obs := NULL]
setkey(data, yyyymm, wficn, permno)

# sanity checks on position sizes
cut <- .2
data[shares_1 / shrout_1 > cut, shares_1 := shrout_1 * cut]
data[shares / shrout_1 > cut, shares := shrout_1 * cut]
rm(cut)

# merge with quarterly domestic equity fund flows
tmp <- readRDS("../../../../data/institutional/fund_flows/quarterly_flow_by_wficn.RDS")

# holdings-implied fund AUM cannot differ by too much (check match quality)
tt <- data[, .(aum_1 = sum(shares_1 * prc_1)), .(yyyymm, wficn)]
tmp <- merge(tmp, tt, by = c("yyyymm", "wficn"))
rm(tt)
tmp[, size_ratio := aum_1 / tna_1]
tmp <- tmp[(size_ratio >= .8) & (size_ratio <= 1.2)][, .(yyyymm, wficn, flow)]
data <- merge(data, tmp, by = c("yyyymm", "wficn"))
rm(tmp)

# fund can't have too few obs or be too small in AUM
tmp <- data[, list(obs = sum(shares_1 > 0), aum_1 = sum(shares_1 * prc_1)), list(yyyymm, wficn)]
tmp <- tmp[(obs >= 20) & (aum_1 >= 10)][, obs := NULL]
data <- merge(data, tmp, by = c("yyyymm", "wficn"))
rm(tmp)

# get portfolio and market-implied weights
data[shares_1 > 0, w_1 := prc_1 * shares_1 / aum_1]
data <- merge(data, data[shares_1 > 0, .(xx = sum(prc_1 * shrout_1)), .(yyyymm, wficn)], by = c("yyyymm", "wficn"))
data[shares_1 > 0, w_mkt_1 := prc_1 * shrout_1 / xx]
data[, xx := NULL]

# data qualit check: one stock cannot dominate the whole portfolio
cut <- .4
tmp <- data[w_1 > 0, list(max_w_1 = max(w_1), max_w_mkt_1 = max(w_mkt_1)), list(yyyymm, wficn)]
tmp <- tmp[(max_w_1 < cut) & (max_w_mkt_1 < cut)]
tmp[, c("max_w_1", "max_w_mkt_1") := NULL]
data <- merge(data, tmp, by = c("yyyymm", "wficn"))
rm(tmp, cut)
gc()

# sort into bins within each (yyyymm, wficn)
data[w_1 > 0, bin := ntile(w_1, 20), .(yyyymm, wficn)]
data[w_1 > 0, bin_to_mkt := ntile(w_1 - w_mkt_1, 20), .(yyyymm, wficn)]
data[w_1 == 0, c("bin", "bin_to_mkt") := NA]

stopifnot(nrow(out) == nrow(data))
data <- copy(out)
rm(out, tmp, blocks, i, p.get_one)
data[, idx := NULL]

# put into easier-to-use formats
data[, dw := (shares - shares_1) * prc_1 / aum_1]
data[, dw := Winsorize(dw, val = quantile(dw, probs = c(.0001, .9999)))]
data[, c("shares", "shares_1", "prc_1", "shrout_1") := NULL]

# save
to_file <- "../tmp/additional/clean_fit/flow_to_trade/regression_table.RDS"
dir.create(dirname(to_file), showWarnings = FALSE, recursive = TRUE)
saveRDS(data, to_file)

# # --- SANITY CHECKS: compare with earlier. Have more data, and is sufficiently close

# data <- readRDS("../tmp/additional/clean_fit/flow_to_trade/regression_table.RDS")
# data_old <- readRDS("../../../20250117_quarterly/tmp/raw_data/cleaning/fit/trade_and_flow_regression_table.RDS")

# dim(data)
# dim(data_old)

# data <- merge(data, data_old, by = c('yyyymm','wficn','permno'), all = T)
# rm(data_old)
# data[, cor(dw.x, dw.y, use = "pairwise.complete.obs")]
# data[, cor(bin.x, bin.y, use = "pairwise.complete.obs")]
