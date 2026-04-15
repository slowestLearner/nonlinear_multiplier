# Contemporaneous regression of returns on demand: panel regressions
library(this.path)
setwd(this.path::this.dir())
source("../utilities/runmefirst.R")
source("../utilities/regressions.R")

# load regression data
tic("loading data")
data_all <- readRDS("../tmp/raw_data/reg_inputs/reg_table_static.RDS")

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


# ###########################################################
# Estimation, with nonlinear or stdev-based specification
# ###########################################################

# control variables for different specifications
control_formulas <- c(
  "1",
  paste0(c(controls_char), collapse = "+"),
  paste0(c(controls_char, controls_liq), collapse = "+")
)

# process one type of data. reg_spec = "nonlinear" or "stdev"
# TODO: delete the nonlinear specification in the final version
p.process_one_type <- function(data, reg_spec = "nonlinear") {
  this_type <- data[1, type]

  out_all <- data.table() # save all results here

  # introduce direct controls
  for (spec_idx in 1:length(control_formulas)) {
    # print(paste0("this_type = ", this_type, " spec_idx = ", spec_idx))
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
    out <- p.panel_regression(data, ff, compare_coefs = ifelse(reg_spec == "nonlinear", FALSE, TRUE))
    out_no_ofi <- p.panel_regression(data, ff_no_ofi, compare_coefs = FALSE)
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
    # print(paste0("this_type = ", this_type, " spec_idx = ", spec_idx))

    out <- p.panel_regression(data, ff, compare_coefs = ifelse(reg_spec == "nonlinear", FALSE, TRUE))
    out_no_ofi <- p.panel_regression(data, ff_no_ofi, compare_coefs = FALSE)
    out[, r2_no_ofi := out_no_ofi[1, r2]]
    out[, spec_idx := spec_idx]
    out_all <- rbind(out_all, out)
  }

  # name the control variables being added
  tmp <- data.table(
    spec_idx = 1:(length(controls_list) + 3),
    var_added = c("none_init", "controls_char", "controls_char+controls_liq", controls_list),
    var_type = c(
      rep("", 3), rep("return-predicting chars", length(controls_char)),
      rep("liquidity", length(controls_liq))
    )
  )
  out_all <- merge(out_all, tmp, by = "spec_idx", all.x = T)
  out_all[, type := this_type]

  return(out_all)
}
toc()


# stdev-based specification---
tic("static panel: stdev")
out_stdev <- rbindlist(mclapply(split(data_all, by = "type"), function(x) {
  p.process_one_type(x, reg_spec = "stdev")
}, mc.cores = nc))
toc()

to_dir <- "../tmp/price_impact/regression_contemp/"
dir.create(to_dir, recursive = T, showWarnings = F)
saveRDS(out_stdev, paste0(to_dir, "panel_stdev.RDS"))

# # -- SANITY check: small differences

# new <- readRDS("../tmp/price_impact/regression_contemp/panel_stdev.RDS")
# old <- readRDS("../tmp/price_impact/regression_contemp_todel/panel_stdev.RDS")
# dim(new) == dim(old)

# setkey(new, spec_idx, var, type)
# setkey(old, spec_idx, var, type)

# out <- merge(new[, .(spec_idx, var, type, coef, se)], old[, .(spec_idx, var, type, coef, se)], by = c("spec_idx", "var", "type"))

# out[, cor(coef.x, coef.y, use = "complete.obs"), type] # WAIT BMI changed more?
# out[type == 'BMI', cor(coef.x, coef.y, use = "complete.obs"), spec_idx] # okay almost entirely about the last specification

# out[, cor(se.x, se.y, use = "complete.obs"), type]
# new[, max(obs), type]
# old[, max(obs), type]
