# Augmenting 2a_regression_fm.R to also control for chars X bins
# TODO: merge this into 2a_regression_fm.R at some point
library(this.path)
setwd(this.path::this.dir())
source("../utilities/runmefirst.R")
source("../utilities/regressions.R")

# load regression data
tic("load data")
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
# Fama-MacBeth with nonlinear or stdev-based specification
# ###########################################################

# control variables for different specifications.
# Keep the first few specifications the same as in 2a_regression_fm.R for comparison
controls_char_bin <- paste0(controls_char, " * as.factor(bin)")
controls_liq_bin <- paste0(controls_liq, " * as.factor(bin)")

control_formulas <- c(
  "1",
  paste0(c(controls_char), collapse = "+"),
  paste0(c(controls_char, controls_liq), collapse = "+"),
  paste0(c(controls_char_bin, controls_liq), collapse = "+"),
  paste0(c(controls_char_bin, controls_liq_bin), collapse = "+")
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


  out_all[, type := this_type]

  return(out_all)
}
toc()

# # nonlinear specification---
# tic("static fm: nonlinear")
# out_nonlinear <- rbindlist(mclapply(split(data_all, by = "type"), function(x) {
#   p.process_one_type(x, reg_spec = "nonlinear")
# }, mc.cores = nc))
# toc()

# dir.create("tmp/price_impact/regression_contemp/", recursive = T, showWarnings = F)

# stdev-based specification---
tic("static fm: stdev")
out_stdev <- rbindlist(mclapply(split(data_all, by = "type"), function(x) {
  p.process_one_type(x, reg_spec = "stdev")
}, mc.cores = nc))
toc()

to_dir <- "../tmp/price_impact/regression_contemp/"
dir.create(to_dir, recursive = T, showWarnings = F)
saveRDS(out_stdev, paste0(to_dir, "fm_stdev_substitution.RDS"))

# # === SANITY check: same as earlier?

# new <- readRDS("tmp/price_impact/regression_contemp/fm_stdev_substitution.RDS")[spec_idx %in% 1:3]
# old <- readRDS("tmp/price_impact/regression_contemp/fm_stdev.RDS")

# compare <- merge(new, old, by = c("type", "var", "spec_idx"))
# compare[, cor(coef.x, coef.y)]
# compare[, cor(se.x, se.y)]

# this_type <- "OFI"
# dcast(new[(type == this_type) & (var %in% paste0("ofi_bin", 1:3))], var ~ spec_idx, value.var = "coef")
# dcast(new[(type == this_type) & (var %in% paste0("ofi_bin", 1:3))], var ~ spec_idx, value.var = "se")

# dcast(new[(type == this_type) & grepl("- ofi_bin", var)], var ~ spec_idx, value.var = "coef")
# dcast(new[(type == this_type) & grepl("- ofi_bin", var)], var ~ spec_idx, value.var = "se")
