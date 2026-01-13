# Dynamic regression of returns on demand. On a 6-core machine took around 1 min
library(this.path)
setwd(this.path::this.dir())
source("../utilities/runmefirst.R")
source("../utilities/regressions.R")

# load regression data
tic("Loading regression data")
data <- readRDS("../tmp/raw_data/reg_inputs/reg_table_dynamic.RDS")

# separate controls and main variables
vars_id <- c("yyyymm", "permno", "type")
vars_reg <- c("ret", "ofi", "cumofi_1", "cumofi_2", "cumofi_3", "cumofi_4")
vars_controls <- setdiff(names(data), c(vars_id, vars_reg))

# put regression-related data into long format
data_controls <- data[, c(vars_id, vars_controls), with = F]
setnames(data_controls, "type", "demand_type")

data_reg <- data[, c(vars_id, vars_reg), with = F]
data_reg <- melt(data_reg, id.vars = c("yyyymm", "permno", "type", "ret", "ofi"), variable.name = "hor", value.name = "cumofi_lag")
data_reg <- data_reg[0 == rowSums(is.na(data_reg))]
data_reg[, hor := as.integer(gsub("cumofi_", "", hor))]
rm(data, vars_id, vars_reg, vars_controls)

# change the name so old code can easily be applied
data_reg[, demand_type := type]
data_reg[, type := paste0(demand_type, "_", hor, "lag")]
data_reg[, sd_ofi := sd(cumofi_lag), .(yyyymm, type)]
data_reg[, bin := 1]
data_reg[abs(cumofi_lag) > sd_ofi, bin := 2]
data_reg[abs(cumofi_lag) > 2 * sd_ofi, bin := 3]
data_reg[, ofi_bin1 := ifelse(bin == 1, ofi, 0)]
data_reg[, ofi_bin2 := ifelse(bin == 2, ofi, 0)]
data_reg[, ofi_bin3 := ifelse(bin == 3, ofi, 0)]
data_reg[, sd_ofi := NULL]
data_reg[, ofi_absofi := ofi * abs(cumofi_lag)]

# get control variable names
cdata <- readRDS("../tmp/raw_data/controls/controls_classification.RDS")
controls_char <- cdata[control_type == "return-predictor", var]
controls_liq <- cdata[control_type == "liquidity", var]
controls_list <- c(controls_char, controls_liq)
rm(cdata)
toc()

# ###########################################################
# Fama-MacBeth with nonlinear or stdev-based specification
# ###########################################################

# control variables for different specifications
control_formulas <- c(
  "1",
  paste0(c(controls_char), collapse = "+"),
  paste0(c(controls_char, controls_liq), collapse = "+")
)

# -- process one type of data. reg_spec = "nonlinear" or "stdev"
# loops through various specifications. Start with adding direct controls (in steps), and then
# start adding interactios with ofi one step at a time
p.process_one_type <- function(data, reg_spec = "nonlinear") {
  data <- merge(data, data_controls, by = c("yyyymm", "permno", "demand_type"), all.x = T)
  data <- data[0 == rowSums(is.na(data))]
  this_type <- data[1, type]

  out_all <- data.table() # save all results here

  # introduce direct controls
  for (spec_idx in 1:length(control_formulas)) {
    if (reg_spec == "nonlinear") {
      ff <- paste0("ret ~ ", control_formulas[spec_idx], " + ofi + ofi_absofi")
    } else if (reg_spec == "stdev") {
      ff <- paste0("ret ~ ", control_formulas[spec_idx], " + ofi_bin1 + ofi_bin2 + ofi_bin3")
    }
    if (this_type == "BMI") {
      ff <- paste0(ff, "+", paste0(controls_bmi, collapse = "+"))
    }
    out <- p.fama_macbeth(data, ff, compare_coefs = ifelse(reg_spec == "nonlinear", FALSE, TRUE))
    out[, spec_idx := spec_idx]
    out_all <- rbind(out_all, out)
  }

  # then add interactions with ofi
  for (this_v in c(controls_char, controls_liq)) {
    setnames(data, this_v, "xx")
    data[, yy := xx * ofi]
    setnames(data, c("xx", "yy"), c(this_v, paste0("ofi_", this_v)))
    ff <- paste0(ff, " + ofi_", this_v)
    spec_idx <- spec_idx + 1

    out <- p.fama_macbeth(data, ff, compare_coefs = ifelse(reg_spec == "nonlinear", FALSE, TRUE))
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

# nonlinear---
tic("dynamic fm: nonlinear")
out_nonlinear <- rbindlist(mclapply(split(data_reg, by = "type"), function(x) {
  p.process_one_type(x, reg_spec = "nonlinear")
}, mc.cores = nc))

# also get cum_d_sd
sd_data <- data_reg[, .(cum_d_sd = sd(cumofi_lag)), .(type)]
out_nonlinear <- merge(out_nonlinear, sd_data, by = "type")
rm(sd_data)

to_dir <- "../tmp/price_impact/regression_dynamic/"
dir.create(to_dir, recursive = T, showWarnings = F)
saveRDS(out_nonlinear, paste0(to_dir, "fm_nonlinear.RDS"))
toc()

# stdev-based specification---
tic("dynamic fm: stdev")
out_stdev <- rbindlist(mclapply(split(data_reg, by = "type"), function(x) {
  p.process_one_type(x, reg_spec = "stdev")
}, mc.cores = nc))

# also get TS average of XS SD
sd_data <- data_reg[, .(cum_d_sd = sd(cumofi_lag)), .(type, yyyymm)][, .(cum_d_sd = mean(cum_d_sd)), .(type)]
out_stdev <- merge(out_stdev, sd_data, by = "type")
rm(sd_data)

# # --- take a look at results
# out <- out_stdev[spec_idx == 3]
# out <- out[grepl("ofi", var)]
# tt <- unique(out[, .(var)])[, var_idx := .I][, var_lab := paste0(var_idx, "_", var)]
# out <- merge(out, tt, by = "var")
# rm(tt)
# out[, type_lab := paste0(type, "_", spec_idx)]
# out[, tstat := coef / se]
# out <- out[grepl("OFI", type)]
# dcast(out, var_lab ~ type, value.var = c("coef"))
# dcast(out, var_lab ~ type, value.var = c("tstat"))

to_dir <- "../tmp/price_impact/regression_dynamic/"
dir.create(to_dir, recursive = T, showWarnings = F)
saveRDS(out_stdev, paste0(to_dir, "fm_stdev.RDS"))
toc()


# # # === SANITY check (i tink old version had wrong standard errors)

# # --- nonlinear specification
# new <- readRDS("tmp/price_impact/regression_dynamic/fm_nonlinear.RDS")
# old <- readRDS("../20250117_quarterly/tmp/price_impact/regression_with_lagged_interaction/fm_quarterly_without_direction.RDS")
# old <- old[type != "OFI"]
# old[type == "OFI_resid", type := "OFI"]

# # choose specification
# new <- new[spec_idx == 3][, spec_idx := NULL]

# # choose variables
# new <- new[var %in% c("ofi", "ofi_absofi")]
# old[var == "I(ofi * abs(cumofi_prev))", var := "ofi_absofi"]
# common_vars <- intersect(unique(new[, var]), unique(old[, var]))
# new <- new[var %in% common_vars]
# old <- old[var %in% common_vars]
# rm(common_vars)

# # get joint variables
# old[, type := paste0(type, "_", lag, "lag")]
# vv <- intersect(names(new), names(old))
# new <- new[, ..vv]
# old <- old[, ..vv]
# rm(vv)

# # adjust sizes
# new[var == "ofi_absofi", coef := coef / 100]
# new[var == "ofi_absofi", se := se / 100]
# new[, cum_d_sd := cum_d_sd * 100]

# # compare
# compare <- merge(new, old, by = c("type", "var"), all = T)
# stopifnot(nrow(compare) == nrow(new))
# rm(new, old)

# compare[, mean(abs(coef.x - coef.y)) / mean(abs(coef.x))]
# compare[, mean(abs(se.x - se.y)) / mean(abs(se.x))]
# compare[, mean(abs(obs.x - obs.y)) / mean(abs(obs.x))]
# compare[, mean(abs(r2.x - r2.y)) / mean(abs(r2.x))]
# compare[, mean(abs(cum_d_sd.x - cum_d_sd.y)) / mean(abs(cum_d_sd.x))]


# # --- stdev-based specification
# new <- readRDS("tmp/price_impact/regression_dynamic/fm_stdev.RDS")
# old <- readRDS("../20250117_quarterly/tmp/price_impact/regression_with_lagged_interaction/fm_quarterly_by_stdev_based_bins.RDS")
# old <- old[type != "OFI"]
# old[type == "OFI_resid", type := "OFI"]

# # choose specification
# new <- new[spec_idx == 3][, spec_idx := NULL]

# # choose variables
# new[, var := gsub("ofi_bin", "M", var)]
# new[, var := gsub(" ", "", var)]
# common_vars <- intersect(unique(new[, var]), unique(old[, var]))
# new <- new[var %in% common_vars]
# old <- old[var %in% common_vars]
# rm(common_vars)

# # get joint variables
# old[, type := paste0(type, "_", lag, "lag")]
# vv <- intersect(names(new), names(old))
# new <- new[, ..vv]
# old <- old[, ..vv]
# rm(vv)

# # adjust coef size
# new[, cum_d_sd := cum_d_sd * 100]

# # compare
# compare <- merge(new, old, by = c("type", "var"), all = T)
# stopifnot(nrow(compare) == nrow(new))
# rm(new, old)

# compare[, mean(abs(coef.x - coef.y)) / mean(abs(coef.x))]
# compare[, mean(abs(se.x - se.y)) / mean(abs(se.x))]
# compare[, mean(abs(obs.x - obs.y)) / mean(abs(obs.x))]
# compare[, mean(abs(r2.x - r2.y)) / mean(abs(r2.x))]
# compare[, mean(abs(cum_d_sd.x - cum_d_sd.y)) / mean(abs(cum_d_sd.x))]

# # not the same
# compare[type == "OFI_4lag", .(se.x / se.y)]

# compare[, cor(coef.x, coef.y)]
# compare[, cor(se.x, se.y)] # this is somehow not identical, but sufficiently similar?
# compare[, cor(obs.x, obs.y)]
# compare[, cor(r2.x, r2.y)]
