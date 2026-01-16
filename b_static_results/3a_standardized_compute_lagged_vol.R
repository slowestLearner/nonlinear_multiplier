# Compare rolling standard deviations of FIT and OFI
library(this.path)
setwd(this.path::this.dir())
source("../utilities/runmefirst.R")

# get FIT and OFI data
data <- readRDS("../tmp/raw_data/reg_inputs/all_ofi_and_ret.RDS")
data <- data[type != "BMI"]

# use up to 12q lags
data[, idx := frank(yyyymm, ties.method = "dense")]
for (i in 1:12) {
  data <- merge(data, data[, .(idx = idx + i, type, permno, xx = ofi)], by = c("idx", "type", "permno"), all.x = T)
  setnames(data, "xx", paste0("ofi_", i))
}

# compute standard errors
ofi_matrix <- as.matrix(data[, .SD, .SDcols = paste0("ofi_", 1:12)])

# 4 lags
X <- copy(ofi_matrix[, 1:4])
idx <- rowSums(is.na(X)) == 0
sd_values <- apply(X[idx, ], 1, sd)
data[idx, sd_ofi_4q := sd_values]

# 8 lags
X <- copy(ofi_matrix[, 1:8])
idx <- rowSums(is.na(X)) == 0
sd_values <- apply(X[idx, ], 1, sd)
data[idx, sd_ofi_8q := sd_values]

# 12 lags
X <- copy(ofi_matrix[, 1:12])
idx <- rowSums(is.na(X)) == 0
sd_values <- apply(X[idx, ], 1, sd)
data[idx, sd_ofi_12q := sd_values]

rm(X, ofi_matrix, idx, sd_values)

# save
data <- data[, .(yyyymm, permno, type, sd_ofi_4q, sd_ofi_8q, sd_ofi_12q)]
data <- data.table(melt(data, id.vars = c("yyyymm", "permno", "type")))
names(data)[4:5] <- c("spec", "std_ofi")
data[, spec := as.character(spec)]

# require having some data
data <- merge(data, unique(data[!is.na(std_ofi), .(yyyymm, permno)]), by = c("yyyymm", "permno"))
data <- data[yyyymm >= 199303]

to_dir <- "../tmp/raw_data/reg_inputs/"
dir.create(to_dir, recursive = T, showWarnings = F)
saveRDS(data, paste0(to_dir, "quarterly_fit_and_ofi_lagged_rolling_stdev.RDS"))


# # === SANITY

# new <- readRDS("../tmp/raw_data/reg_inputs/quarterly_fit_and_ofi_lagged_rolling_stdev.RDS")
# old <- readRDS("../../20250117_quarterly/tmp/raw_data/reg_inputs/quarterly_fit_and_ofi_lagged_rolling_stdev.RDS")
# old[type == "OFI_resid", type := "OFI_pre_whitened"]

# new <- new[0 == rowSums(is.na(new))]
# old <- old[0 == rowSums(is.na(old))]
# compare <- merge(new, old, by = c("yyyymm", "permno", "type", "spec"))
# compare[, cor(std_ofi.x, std_ofi.y), type]
