# --- compute stdev and (modified) HHI for quarterly OFI
library(this.path)
setwd(this.path::this.dir())
source("../../utilities/runmefirst.R")

# daily OFI
tic("read data")
data <- readRDS("../../../../../data/demand_shocks/ofi/daily.RDS")
tmp <- unique(data[, .(date)])
tmp[, yyyymm := 100 * year(date) + 3 * quarter(date)]
data <- merge(data, tmp, by = "date")
rm(tmp)

# # get rid of outliers
# data[, ofi := Winsorize(ofi, quantile(ofi, probs = c(.005, .995)))]

# compute HHI in the same direction
data[, sum_ofi := sum(ofi), by = .(yyyymm, permno)]
data[, abs_ofi_same_dir := ifelse(sign(ofi) == sign(sum_ofi), abs(ofi), 0)]
data[, denom := sum(abs_ofi_same_dir), .(yyyymm, permno)]
data[, w := abs_ofi_same_dir / denom]
data <- data[denom > .001]

# make sure weights sum to 1
stopifnot(min(data[, sum(w), .(yyyymm, permno)][, V1]) > .999999)
stopifnot(max(data[, sum(w), .(yyyymm, permno)][, V1]) < 1.000001)
toc()

tic("summarize and save")
# summarize by (yyyymm, permno)
out <- data[, .(
  ofi_sum = sum(ofi),
  abs_ofi_sum = sum(abs(ofi)),
  obs = .N,
  hhi_pow1.5 = sum(w^1.5)^(1 / (1.5 - 1)),
  hhi_pow2 = sum(w^2)^(1 / (2 - 1)),
  hhi_pow3 = sum(w^3)^(1 / (3 - 1)),
  hhi_pow5 = sum(w^5)^(1 / (5 - 1)),
  hhi_pow10 = sum(w^10)^(1 / (10 - 1)),
  ofi_sd = sd(ofi, na.rm = T)
), .(yyyymm, permno)] %>%
  na.omit() %>%
  setDT()

# just keep the subset with quarterly OFI
tmp <- readRDS("../../tmp/raw_data/reg_inputs/reg_table_static.RDS")[type == "OFI", .(yyyymm, permno)]
out <- merge(out, tmp, by = c("yyyymm", "permno"))
rm(tmp)

to_file <- "tmp/raw_files/ofi_dispersion.RDS"
dir.create(dirname(to_file), showWarnings = FALSE, recursive = TRUE)
saveRDS(out, to_file)
toc()
