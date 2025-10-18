# Put together data for static regressions
source("utilities/runmefirst.R")

# All demand and returns
data <- readRDS("tmp/raw_data/reg_inputs/all_ofi_and_ret.RDS")
data <- data[yyyymm >= 199306] # liquidity-characteristics are not available for 199303

# add controls that are specific to BMI
tmp <- readRDS("tmp/raw_data/controls/controls_for_BMI.RDS")
controls_bmi <- setdiff(names(tmp), c("yyyymm", "permno")) # variable names for BMI-specific controls
data <- merge(data, tmp, by = c("yyyymm", "permno"), all.x = T)
rm(tmp)

# merge in other controls
tmp <- readRDS("tmp/raw_data/controls/quarterly_controls_lagged.RDS")
data <- merge(data, tmp, by = c("yyyymm", "permno"), all.x = T)
rm(tmp)

# get the names of other controls
cdata <- readRDS("tmp/raw_data/controls/controls_classification.RDS")
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
data[, ofi_absofi := ofi * abs(ofi)]
data[, sd_ofi := sd(ofi), .(yyyymm, type)]
data[, bin := 1]
data[abs(ofi) > sd_ofi, bin := 2]
data[abs(ofi) > 2 * sd_ofi, bin := 3]
data[, ofi_bin1 := ifelse(bin == 1, ofi, 0)]
data[, ofi_bin2 := ifelse(bin == 2, ofi, 0)]
data[, ofi_bin3 := ifelse(bin == 3, ofi, 0)]
data[, sd_ofi := NULL]

saveRDS(data, "tmp/raw_data/reg_inputs/reg_table_static.RDS")
