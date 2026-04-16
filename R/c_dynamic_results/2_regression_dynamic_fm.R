# --- Dynamic regression of returns on demand
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
data_reg <- data.table::melt(data_reg, id.vars = c("yyyymm", "permno", "type", "ret", "ofi"), variable.name = "hor", value.name = "cumofi_lag")
data_reg <- data_reg[0 == rowSums(is.na(data_reg))]
data_reg[, hor := as.integer(gsub("cumofi_", "", hor))]
rm(data, vars_id, vars_reg, vars_controls)

# change some variable names
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
# Fama-MacBeth
# ###########################################################

# control variables for different specifications
control_formulas <- c(
  "1",
  paste0(c(controls_char), collapse = "+"),
  paste0(c(controls_char, controls_liq), collapse = "+")
)

# worker function
# currently, on a 6 core machine, this takes 4 mins to run
# TODO: remove nonlinear
# TODO: this code takes quite a while to run, as it really only uses the last spec... consider shortening?
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
    # if (this_type == "BMI") {
    #   ff <- paste0(ff, "+", paste0(controls_bmi, collapse = "+"))
    # }
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


tic("dynamic regression")
out_stdev <- rbindlist(mclapply(split(data_reg, by = "type"), function(x) {
  p.process_one_type(x, reg_spec = "stdev")
}, mc.cores = nc))

# also get time-series average of cross-sectional standard eviations
sd_data <- data_reg[, .(cum_d_sd = sd(cumofi_lag)), .(type, yyyymm)][, .(cum_d_sd = mean(cum_d_sd)), .(type)]
out_stdev <- merge(out_stdev, sd_data, by = "type")
rm(sd_data)

to_dir <- "../tmp/price_impact/regression_dynamic/"
dir.create(to_dir, recursive = T, showWarnings = F)
saveRDS(out_stdev, paste0(to_dir, "fm_stdev.RDS"))
toc()


# # --- SANITY CHECK: close enough

# new <- readRDS("../tmp/price_impact/regression_dynamic/fm_stdev.RDS")
# old <- readRDS("../tmp/price_impact/regression_dynamic_todel/fm_stdev.RDS")
# dim(new) == dim(old)

# out <- merge(new, old, by = c("spec_idx", "var", "type"))
# out[, cor(coef.x, coef.y, use = "complete.obs"), type]
# out[, cor(se.x, se.y, use = "complete.obs"), type]
# out[, max(obs.x) / max(obs.y), type]
