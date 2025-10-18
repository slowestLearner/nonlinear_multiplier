# Put together data for dynamic regressions
source("utilities/runmefirst.R")

# Start with static demand for FIT and OFI
data <- readRDS("tmp/raw_data/reg_inputs/reg_table_static.RDS") %>% filter(type != "BMI")

# get rid of the static nonlinear variables
data[, c("ofi_absofi", "ofi_bin1", "ofi_bin2", "ofi_bin3", "bin") := NULL]

# get lagged demand over the previous 4 quarters
data[, idx := frank(yyyymm, ties.method = "dense")]
for (i in 1:4) {
  data <- merge(data, data[, .(idx = idx + i, permno, type, xx = ofi)], by = c("idx", "permno", "type"), all.x = T)
  setnames(data, "xx", paste0("ofi_", i))
}
data[, idx := NULL]

# get cumulative past demand
data[, cumofi_1 := ofi_1]
for (i in 2:4) {
  setnames(data, paste0("cumofi_", (i - 1)), "xx")
  setnames(data, paste0("ofi_", i), "yy")
  data[, zz := xx + yy]
  setnames(data, "xx", paste0("cumofi_", (i - 1)))
  setnames(data, "yy", paste0("ofi_", i))
  setnames(data, "zz", paste0("cumofi_", i))
}
data[, paste0("ofi_", 1:4) := NULL]

saveRDS(data, "tmp/raw_data/reg_inputs/reg_table_dynamic.RDS")
