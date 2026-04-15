# Put together data for static regressions
library(this.path)
setwd(this.path::this.dir())
# source("utilities/runmefirst.R")

# All demand and returns
data <- readRDS("../tmp/raw_data/reg_inputs/all_ofi_and_ret.RDS")
data <- data[yyyymm >= 199306] # liquidity-characteristics are not available for 199303

# add controls that are specific to BMI
tmp <- readRDS("../tmp/raw_data/controls/controls_for_BMI.RDS")
controls_bmi <- setdiff(names(tmp), c("yyyymm", "permno")) # variable names for BMI-specific controls
data <- merge(data, tmp, by = c("yyyymm", "permno"), all.x = T)
rm(tmp)

# merge in other controls
tmp <- readRDS("../tmp/raw_data/controls/quarterly_controls_lagged.RDS")
data <- merge(data, tmp, by = c("yyyymm", "permno"), all.x = T)
rm(tmp)

# get the names of other controls
cdata <- readRDS("../tmp/raw_data/controls/controls_classification.RDS")
controls_liq <- cdata[control_type == "liquidity", var]
controls_char <- cdata[control_type == "return-predictor", var]
controls_ind <- cdata[control_type == "industry", var]
rm(cdata)

# don't use industries (TODO: streamline code later)
data[, c(controls_ind) := NULL]
rm(controls_ind)

# Re-center characteristics within each time period and demand type
for (this_v in c(controls_liq, controls_char)) {
  print(paste0("recentering control = ", this_v))
  setnames(data, this_v, "xx")
  data[, xx := xx - mean(xx, na.rm = T), .(yyyymm, type)]
  setnames(data, "xx", this_v)
}
rm(this_v)
data[is.na(data)] <- 0 # fill zeros if needed (for these characteristics)
controls_list <- c(controls_char, controls_liq)

# get nonlinear x variables
data_pos <- copy(data)

# keep positive OFI only
data_pos[, ofi := fifelse(ofi > 0, ofi, 0)]

# nonlinear terms
data_pos[, ofi_absofi := ofi * abs(ofi)]
data_pos[, sd_ofi := sd(ofi), .(yyyymm, type)]

data_pos[, bin := 1L]
data_pos[abs(ofi) > sd_ofi,        bin := 2L]
data_pos[abs(ofi) > 2 * sd_ofi,    bin := 3L]

data_pos[, ofi_bin1 := fifelse(bin == 1L, ofi, 0)]
data_pos[, ofi_bin2 := fifelse(bin == 2L, ofi, 0)]
data_pos[, ofi_bin3 := fifelse(bin == 3L, ofi, 0)]

data_pos[, sd_ofi := NULL]


data_neg <- copy(data)

# keep negative OFI only (preserve sign)
data_neg[, ofi := fifelse(ofi < 0, ofi, 0)]

# nonlinear terms
data_neg[, ofi_absofi := ofi * abs(ofi)]
data_neg[, sd_ofi := sd(ofi), .(yyyymm, type)]

data_neg[, bin := 1L]
data_neg[abs(ofi) > sd_ofi,        bin := 2L]
data_neg[abs(ofi) > 2 * sd_ofi,    bin := 3L]

data_neg[, ofi_bin1 := fifelse(bin == 1L, ofi, 0)]
data_neg[, ofi_bin2 := fifelse(bin == 2L, ofi, 0)]
data_neg[, ofi_bin3 := fifelse(bin == 3L, ofi, 0)]

data_neg[, sd_ofi := NULL]


to_dir <- "../tmp/raw_data/reg_inputs/"
dir.create(to_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(data_pos, file.path(to_dir, "reg_table_static_pos_ofi.RDS"))
saveRDS(data_neg, file.path(to_dir, "reg_table_static_neg_ofi.RDS"))
