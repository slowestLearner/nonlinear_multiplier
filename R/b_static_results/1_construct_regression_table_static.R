# ---- Put together a panel dataset for static regressions
library(this.path)
setwd(this.path::this.dir())
source("../utilities/runmefirst.R")

# All demand and returns
data <- readRDS("../tmp/raw_data/reg_inputs/all_ofi_and_ret.RDS")[yyyymm %in% 199306:202212]

# add controls that are specific to BMI
tmp <- readRDS("../tmp/raw_data/controls/controls_for_BMI.RDS")
controls_bmi <- setdiff(names(tmp), c("yyyymm", "permno")) # variable names for BMI-specific controls
data <- merge(data, tmp, by = c("yyyymm", "permno"), all.x = T)
rm(tmp)

# get the names of other controls
cdata <- readRDS("../tmp/raw_data/controls/controls_classification.RDS")
controls_liq <- cdata[control_type == "liquidity", var]
controls_char <- cdata[control_type == "return-predictor", var]
rm(cdata)

# merge in other controls
tmp <- readRDS("../tmp/raw_data/controls/quarterly_controls_lagged.RDS")
tmp <- tmp[, c("yyyymm", "permno", controls_liq, controls_char), with = F]
data <- merge(data, tmp, by = c("yyyymm", "permno"), all.x = T)
rm(tmp)

# Re-center characteristics within each time period and demand type
for (this_v in c(controls_liq, controls_char)) {
  print(paste0("recentering control = ", this_v))
  setnames(data, this_v, "xx")
  data[, xx := xx - mean(xx, na.rm = T), .(yyyymm, type)]
  setnames(data, "xx", this_v)
}
rm(this_v)
data[is.na(data)] <- 0 # fill zeros if needed for characteristics
controls_list <- c(controls_char, controls_liq)

# create nonlinear demand variables
data[, ofi_absofi := ofi * abs(ofi)]
data[, sd_ofi := sd(ofi), .(yyyymm, type)]
data[, bin := 1]
data[abs(ofi) > sd_ofi, bin := 2]
data[abs(ofi) > 2 * sd_ofi, bin := 3]
data[, ofi_bin1 := ifelse(bin == 1, ofi, 0)]
data[, ofi_bin2 := ifelse(bin == 2, ofi, 0)]
data[, ofi_bin3 := ifelse(bin == 3, ofi, 0)]
data[, sd_ofi := NULL]

# save
setkey(data, yyyymm, permno, type)

to_file <- "../tmp/raw_data/reg_inputs/reg_table_static.RDS"
dir.create(dirname(to_file), recursive = T, showWarnings = F)
saveRDS(data, to_file)

# # --- SANITY CHECK: identical except for pre-whitened OFI
# tmp <- readRDS(to_file)
# setkey(data, yyyymm, permno, type)
# setkey(tmp, yyyymm, permno, type)
# mean(data == tmp, na.rm = T)
# mean(is.na(data) == is.na(tmp))

# out <- merge(data, tmp, by = c("yyyymm", "permno", "type"))
# out[, cor(ofi.x, ofi.y, use = "complete.obs"), type]
# out[, cor(bin.x, bin.y, use = "complete.obs"), type]
# out[, cor(mom.x, mom.y, use = "complete.obs"), type]
# out[, cor(realized_vol.x, realized_vol.y, use = "complete.obs"), type] # some directions flipped but it does not matter for results
