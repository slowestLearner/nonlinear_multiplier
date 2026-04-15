# Estimate heterogeneous trading response to flows.
library(this.path)
setwd(this.path::this.dir())
source("../utilities/runmefirst.R")

# regression table
data <- readRDS("../tmp/additional/clean_fit/flow_to_trade/regression_table.RDS")

# sort flows into bins
n_bins <- 20
cut <- .01
tmp <- data[, .(flow = last(flow)), .(yyyymm, wficn)]
tmp[abs(flow) < cut, bin_flow := 0]
tmp[flow > cut, bin_flow := ntile(flow, n_bins)]
tmp[flow < -cut, bin_flow := ntile(flow, n_bins) - (n_bins + 1)]
tmp[, flow := NULL]
data <- merge(data[, bin_flow := NULL], tmp, by = c("yyyymm", "wficn"))
rm(tmp)

# get flow-size bins
data[, bin_composite := paste0(bin, "_", bin_flow)]
gc()

# main LHS variable (trade)
data[w_1 > 0, dw_to_w_1 := dw / w_1]
data[, dw_to_w_1 := Winsorize(dw_to_w_1, val = quantile(dw_to_w_1, probs = c(.001, .999), na.rm = T))]

# - first step: regress out the effect of pre-existing position bins
ols <- feols(dw_to_w_1 ~ 1 | bin, data[w_1 > 0], weights = data[w_1 > 0, w_1], cluster = c("yyyymm", "wficn"))
data[(w_1 > 0) & !is.na(dw_to_w_1), dw_to_w_1_resid := ols$residuals]

# - second step, keep track of the residual trade response to flow by flow-size bins
ols <- feols(dw_to_w_1_resid ~ 1 | bin_composite, data[w_1 > 0], weights = data[w_1 > 0, w_1], cluster = c("yyyymm", "wficn"))
ff <- fixef(ols)

out <- data.table(vv = names(ff$bin_composite), coef = ff$bin_composite)
out[, loc := sapply(vv, function(x) {
  gregexpr("_", x)[[1]][1]
})]
out[, bin := as.integer(substr(vv, 1, loc - 1))]
out[, bin_flow := as.integer(substr(vv, loc + 1, nchar(vv)))]
out[, c("vv", "loc") := NULL]
out[, lab := paste0("w_1 bin = ", bin)]
rm(ff, ols)

# merge with flow sizes and position sizes
tt <- data[, .(flow = mean(flow), w_1 = mean(w_1)), .(bin, bin_flow)]
out <- merge(out, tt, by = c("bin", "bin_flow"))
out <- merge(out, out[, .(w_1_mean = mean(w_1)), bin], by = "bin")
out[, lab := paste0("w_1 = ", round(w_1_mean * 100, 2), "%")]
rm(tt)

to_file <- "../tmp/additional/clean_fit/flow_to_trade/heterogeneous_trade_to_flow.RDS"
dir.create(dirname(to_file), showWarnings = FALSE, recursive = TRUE)
saveRDS(out, to_file)

# # --- SANITY CHECK: sufficiently close to earlier

# out <- readRDS("../tmp/additional/clean_fit/flow_to_trade/heterogeneous_trade_to_flow.RDS")
# out_old <- readRDS("../../../20250117_quarterly/tmp/raw_data/cleaning/fit/heterogeneous_trade_to_flow.RDS")

# dim(out) == dim(out_old)
# out <- merge(out, out_old, by = c("bin", "bin_flow"), all = T)
# out[, cor(w_1_mean.x, w_1_mean.y, use = "pairwise.complete.obs")]
# out[, cor(coef.x, coef.y, use = "pairwise.complete.obs")]
