# Contemporaneous regression of returns on demand, for alternative FIT constructions
library(this.path)
setwd(this.path::this.dir())
source("../utilities/runmefirst.R")
source("../utilities/regressions.R")

# static regression table: swap FIT in these with new data
data <- readRDS("../tmp/raw_data/reg_inputs/reg_table_static.RDS")[type == "FIT"][, type := NULL]

# remove existing demand variables
vv <- names(data)
vv <- vv[grepl("ofi", vv)]
data[, c(vv, "bin") := NULL]

# append new FIT data
# tmp <- readRDS("../../../data/demand_shocks/j_fit/quarterly_residuals.RDS")[, .(yyyymm, permno, origin, type = flow_type, ofi = fit2shrout_cut01)]
tmp <- readRDS("../tmp/additional/clean_fit/flow_residuals/2_fit.RDS")[, .(yyyymm, permno, type = flow_type, ofi = fit2shrout_cut01)]
data <- merge(data, tmp, by = c("yyyymm", "permno"), allow.cartesian = T)
rm(tmp)

# get control variable names
cdata <- readRDS("../tmp/raw_data/controls/controls_classification.RDS")
controls_char <- cdata[control_type == "return-predictor", var]
controls_liq <- cdata[control_type == "liquidity", var]
controls_list <- c(controls_char, controls_liq)
rm(cdata)

# get names for bmi controls
tmp <- readRDS("../tmp/raw_data/controls/controls_for_BMI.RDS")
controls_bmi <- setdiff(names(tmp), c("yyyymm", "permno"))
rm(tmp)

# create bins again based on new FIT data
data[, sd_ofi := sd(ofi), .(yyyymm, type)]
data[, bin := 1]
data[abs(ofi) > sd_ofi, bin := 2]
data[abs(ofi) > 2 * sd_ofi, bin := 3]
data[, ofi_bin1 := ifelse(bin == 1, ofi, 0)]
data[, ofi_bin2 := ifelse(bin == 2, ofi, 0)]
data[, ofi_bin3 := ifelse(bin == 3, ofi, 0)]
data[, sd_ofi := NULL]
data_all <- copy(data)
rm(data)

# ###########################################################
# Fama-MacBeth with nonlinear or stdev-based specification
# ###########################################################

# control variables for different specifications
control_formulas <- c(
  "1",
  paste0(c(controls_char), collapse = "+"),
  paste0(c(controls_char, controls_liq), collapse = "+")
)

# process one type of data. reg_spec = "nonlinear" or "stdev"
p.process_one_type <- function(data, reg_spec = "nonlinear") {
  this_type <- data[1, type]

  out_all <- data.table() # save all results here

  # introduce direct controls
  for (spec_idx in 1:length(control_formulas)) {
    ff_no_ofi <- paste0("ret ~ ", control_formulas[spec_idx])
    if (reg_spec == "nonlinear") {
      ff <- paste0("ret ~ ", control_formulas[spec_idx], " + ofi + ofi_absofi")
    } else if (reg_spec == "stdev") {
      ff <- paste0("ret ~ ", control_formulas[spec_idx], " + ofi_bin1 + ofi_bin2 + ofi_bin3")
    }
    if (this_type == "BMI") {
      ff <- paste0(ff, "+", paste0(controls_bmi, collapse = "+"))
      ff_no_ofi <- paste0(ff_no_ofi, "+", paste0(controls_bmi, collapse = "+"))
    }
    out <- p.fama_macbeth(data, ff, compare_coefs = ifelse(reg_spec == "nonlinear", FALSE, TRUE))
    out_no_ofi <- p.fama_macbeth(data, ff_no_ofi, compare_coefs = FALSE)

    out[, r2_no_ofi := out_no_ofi[1, r2]]
    out[, spec_idx := spec_idx]
    out_all <- rbind(out_all, out)
  }

  # then add interactions with ofi
  for (this_v in c(controls_char, controls_liq)) {
    setnames(data, this_v, "xx")
    data[, yy := xx * ofi]
    setnames(data, c("xx", "yy"), c(this_v, paste0("ofi_", this_v)))
    ff <- paste0(ff, " + ofi_", this_v)
    ff_no_ofi <- paste0(ff_no_ofi, " + ofi_", this_v)
    spec_idx <- spec_idx + 1

    out <- p.fama_macbeth(data, ff, compare_coefs = ifelse(reg_spec == "nonlinear", FALSE, TRUE))
    out_no_ofi <- p.fama_macbeth(data, ff_no_ofi, compare_coefs = FALSE)
    out[, r2_no_ofi := out_no_ofi[1, r2]]
    out[, spec_idx := spec_idx]
    out_all <- rbind(out_all, out)
  }

  # name the control variables being added
  tmp <- data.table(
    spec_idx = 1:(length(controls_list) + 4),
    var_added = c("none_init", "controls_char", "controls_char+controls_liq", "none", controls_list),
    var_type = c(
      rep("", 4), rep("return-predicting chars", length(controls_char)),
      rep("liquidity", length(controls_liq))
    )
  )
  out_all <- merge(out_all, tmp, by = "spec_idx", all.x = T)
  out_all[, type := this_type]
  # out_all[, origin := this_origin]

  return(out_all)
}


# # nonlinear specification---
# tic("static fm: nonlinear")
# out_nonlinear <- rbindlist(mclapply(split(data_all, by = "type"), function(x) {
#   p.process_one_type(x, reg_spec = "nonlinear")
# }, mc.cores = nc))
# toc()

# to_dir <- "../tmp/price_impact/regression_contemp/"
# dir.create(to_dir, recursive = T, showWarnings = F)
# saveRDS(out_nonlinear, paste0(to_dir, "fm_nonlinear.RDS"))

# stdev-based specification---
# tic("static fm: stdev")
# out_stdev <- rbindlist(mclapply(split(data_all, by = c("type", "origin")), function(x) {
#   p.process_one_type(x, reg_spec = "stdev")
# }, mc.cores = nc))
# toc()

# takes min mins on a 6 core computer
tic("parallel processing: stdev specification")
plan(multisession, workers = detectCores() - 2)
out_stdev <- rbindlist(future_lapply(split(data_all, by = c("type")), function(x) {
  p.process_one_type(x, reg_spec = "stdev")
}, future.packages = c("data.table"), future.seed = 123))
plan(sequential)
toc()

to_dir <- "../tmp/price_impact/regression_contemp/alternative_fit/"
dir.create(to_dir, recursive = T, showWarnings = F)
saveRDS(out_stdev, paste0(to_dir, "fm_stdev.RDS"))

# # === SANITY check with earlier

# new <- readRDS("../tmp/price_impact/regression_contemp/alternative_fit/fm_stdev.RDS")
# old <- readRDS("../tmp/price_impact/regression_contemp/alternative_fit_todel/fm_stdev.RDS")[origin == "flow"][, origin := NULL]
# dim(new) == dim(old)

# out <- merge(new, old, by = c("type", "var", "spec_idx"), all = T)
# out[, cor(coef.x, coef.y, use = "complete.obs")]
# out[, cor(se.x, se.y, use = "complete.obs")]
