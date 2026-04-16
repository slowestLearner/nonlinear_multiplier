# ----- Put together panel data for dynamic regressions
library(this.path)
setwd(this.path::this.dir())
source("../utilities/runmefirst.R")

# static demand for FIT and OFI
data <- readRDS("../tmp/raw_data/reg_inputs/reg_table_static.RDS") %>% filter(type != "BMI")

# get rid of the static nonlinear variables
data[, c("ofi_absofi", "ofi_bin1", "ofi_bin2", "ofi_bin3", "bin") := NULL]

# get lagged demand over the previous 4 quarters
data[, idx := frank(yyyymm, ties.method = "dense")]
for (i in 1:4) {
  data <- merge(data, data[, .(idx = idx + i, permno, type, xx = ofi)], by = c("idx", "permno", "type"), all.x = T)
  setnames(data, "xx", paste0("ofi_", i))
}
data[, idx := NULL]

# compute cumulative past demand
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

to_dir <- "../tmp/raw_data/reg_inputs/"
dir.create(to_dir, recursive = T, showWarnings = F)
saveRDS(data, paste0(to_dir, "reg_table_dynamic.RDS"))

# # --- SANITY CHECK: sufficiently similar except OFI-pre-whitened

# new <- readRDS('../tmp/raw_data/reg_inputs/reg_table_dynamic.RDS')
# old <- readRDS('../tmp/raw_data/reg_inputs_todel/reg_table_dynamic.RDS')
# names(new) == names(old)

# out <- merge(new, old, by = c("yyyymm", "permno", "type"))
# out[, cor(ofi.x, ofi.y, use = "complete.obs"), type]
# out[, cor(cumofi_1.x, cumofi_1.y, use = "complete.obs"), type]
# out[, cor(cumofi_2.x, cumofi_2.y, use = "complete.obs"), type]
# out[, cor(cumofi_3.x, cumofi_3.y, use = "complete.obs"), type]
# out[, cor(cumofi_4.x, cumofi_4.y, use = "complete.obs"), type]
