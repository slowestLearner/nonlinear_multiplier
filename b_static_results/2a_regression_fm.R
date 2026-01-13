# Contemporaneous regression of returns on demand
library(this.path)
setwd(this.path::this.dir())
source("../utilities/runmefirst.R")
source("../utilities/regressions.R")

# load regression data
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


# nonlinear specification---
tic("static fm: nonlinear")
out_nonlinear <- rbindlist(mclapply(split(data_all, by = "type"), function(x) {
  p.process_one_type(x, reg_spec = "nonlinear")
}, mc.cores = nc))
toc()

to_dir <- "../tmp/price_impact/regression_contemp/"
dir.create(to_dir, recursive = T, showWarnings = F)
saveRDS(out_nonlinear, paste0(to_dir, "fm_nonlinear.RDS"))

# stdev-based specification---
tic("static fm: stdev")
out_stdev <- rbindlist(mclapply(split(data_all, by = "type"), function(x) {
  p.process_one_type(x, reg_spec = "stdev")
}, mc.cores = nc))
toc()

# # take a look at results
# tt <- copy(out_stdev[spec_idx %in% 1:3])
# tt <- tt[grepl("ofi", var)]
# tt2 <- unique(tt[, .(var)])[, var_idx := .I]
# tt2[, var_lab := paste0(var_idx, "_", var)]
# tt <- merge(tt, tt2, by = "var")
# rm(tt2)
# tt[, type_lab := paste0(type, "_", spec_idx)]
# tt[, tstat := coef / se]

# options(width = 200)
# dcast(tt[grepl("OFI", type)], var_lab ~ type_lab, value.var = c("coef"))
# dcast(tt[grepl("OFI", type)], var_lab ~ type_lab, value.var = c("tstat"))

to_dir <- "../tmp/price_impact/regression_contemp/"
dir.create(to_dir, recursive = T, showWarnings = F)
saveRDS(out_stdev, paste0(to_dir, "fm_stdev.RDS"))


# # === SANITY check

# # --- nonlinear
# new <- readRDS("tmp/price_impact/regression_contemp/fm_nonlinear.RDS")
# old <- readRDS("../20250117_quarterly/tmp/price_impact/regression_contemp/full_sample.RDS")
# old <- old[type != "OFI"]
# old[type == "OFI_resid", type := "OFI"]
# setnames(old, "idx", "spec_idx")

# # get joint variables
# vv <- intersect(names(new), names(old))
# new <- new[, ..vv]
# old <- old[, ..vv]
# rm(vv)

# # compare
# compare <- merge(new, old, by = c("type", "var", "spec_idx"), all = T)
# rm(new, old)
# compare[, mean(abs(coef.x - coef.y)) / mean(abs(coef.x))]
# compare[, mean(abs(se.x - se.y)) / mean(abs(se.x))]
# compare[, mean(abs(obs.x - obs.y)) / mean(abs(obs.x))]
# compare[, mean(abs(r2.x - r2.y)) / mean(abs(r2.x))]
# compare[, mean(abs(r2_no_ofi.x - r2_no_ofi.y)) / mean(abs(r2_no_ofi.x))]


# --- stdev
# new <- readRDS("tmp/price_impact/regression_contemp/fm_stdev.RDS")
# old <- readRDS("../20250117_quarterly/tmp/price_impact/multiplier_by_shock_size_quarterly/fm_by_stdev_based_bins.RDS")
# old <- old[type != "OFI"]
# old[type == "OFI_resid", type := "OFI"]
# setnames(old, "idx", "spec_idx")

# # get joint variables
# vv <- intersect(names(new), names(old))
# new <- new[, ..vv]
# old <- old[, ..vv]
# rm(vv)

# # rename variables to compare
# new <- new[grepl("ofi_bin", var)]
# new[, var := gsub("ofi_bin", "M", var)]
# new[, var := gsub(" ", "", var)]
# stopifnot(nrow(new) == nrow(old))

# # compare
# compare <- merge(new, old, by = c("type", "var", "spec_idx"), all = T)
# rm(new, old)
# compare[, mean(abs(coef.x - coef.y)) / mean(abs(coef.x))]
# compare[, mean(abs(se.x - se.y)) / mean(abs(se.x))]
# compare[, mean(abs(obs.x - obs.y)) / mean(abs(obs.x))]
# compare[, mean(abs(r2.x - r2.y)) / mean(abs(r2.x))]
# compare[, mean(abs(r2_no_ofi.x - r2_no_ofi.y)) / mean(abs(r2_no_ofi.x))]
