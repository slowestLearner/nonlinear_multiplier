# --- Examine the impact of shock signs as we have increasingly stringent controls
library(this.path)
setwd(this.path::this.dir())
source("../../utilities/runmefirst.R")
source("../../utilities/regressions.R")
options(width = 200)

# load regression data
tic("Loading regression data")
data <- readRDS("../../tmp/raw_data/reg_inputs/reg_table_dynamic.RDS")[type %in% c("FIT", "OFI")]

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
data_reg[, ofi_is_same_sign := ifelse(sign(ofi) == sign(cumofi_lag), ofi, 0)]

# get control variable names
cdata <- readRDS("../../tmp/raw_data/controls/controls_classification.RDS")
controls_char <- cdata[control_type == "return-predictor", var]
controls_liq <- cdata[control_type == "liquidity", var]
controls_list <- c(controls_char, controls_liq)
rm(cdata)

# --- create various types of controls
data_reg[, sd_ofi := sd(ofi), .(type, yyyymm)]
data_reg[, sd_ofi_lag := sd(cumofi_lag), .(type, yyyymm)]

# the 3-bin type controls as in the current paper
data_reg[, bin := ifelse(abs(ofi) > 2 * sd_ofi, 3, ifelse(abs(ofi) > sd_ofi, 2, 1))]
data_reg[, bin_lag := ifelse(abs(cumofi_lag) > 2 * sd_ofi_lag, 3, ifelse(abs(cumofi_lag) > sd_ofi_lag, 2, 1))]

# more refined bins
data_reg[, bin_refined := ntile(abs(ofi), 10), .(yyyymm, type)]
data_reg[, bin_lag_refined := ntile(abs(cumofi_lag), 10), .(yyyymm, type)]

# can also chop up into multiple periods
tmp <- unique(data_reg[, .(yyyymm)])[, yyyy := floor(yyyymm / 100)]
tt <- unique(tmp[, .(yyyy)])[, bin_period := ntile(yyyy, 3)][, lab_period := paste0(min(yyyy), "-", max(yyyy)), bin_period]
tmp <- merge(tmp, tt, by = "yyyy")[, yyyy := NULL]
data_reg <- merge(data_reg, tmp, by = c("yyyymm"))
rm(tmp)

# print(unique(data_reg[, .(bin_period, lab_period)]))
#    bin_period lab_period
#         <int>     <char>
# 1:          1  1993-2002
# 2:          2  2003-2012
# 3:          3  2013-2022

data_reg[, ofi_is_same_sign_period1 := ifelse(bin_period == 1 & sign(ofi) == sign(cumofi_lag), ofi, 0)]
data_reg[, ofi_is_same_sign_period2 := ifelse(bin_period == 2 & sign(ofi) == sign(cumofi_lag), ofi, 0)]
data_reg[, ofi_is_same_sign_period3 := ifelse(bin_period == 3 & sign(ofi) == sign(cumofi_lag), ofi, 0)]
toc()

# ###########################################################
# Fama-MacBeth
# ###########################################################

# worker function
# data <- data_list[[7]]
p.process_one_type <- function(data) {
  data <- merge(data, data_controls, by = c("yyyymm", "permno", "demand_type")) %>% na.omit()
  this_type <- data[1, type]

  controls <- paste0(c(controls_char, controls_liq), collapse = "+")

  # vars_contemp <- paste0(c(paste0("ofi_bin", 1:10)), collapse = " + ")
  # vars_lags <- paste0(c(paste0("ofi_bin_lag", 1:10)), collapse = " + ")

  # specifications
  spec_formulas <- c(
    # "ofi_is_same_sign + ofi",
    "ofi_is_same_sign + ofi * as.factor(bin)",
    "ofi_is_same_sign + ofi * as.factor(bin) + ofi * as.factor(bin_lag)",
    "ofi_is_same_sign + ofi * as.factor(bin_refined) + ofi * as.factor(bin_lag_refined)",
    "ofi_is_same_sign_period1 + ofi_is_same_sign_period2 + ofi_is_same_sign_period3 + ofi * as.factor(bin_refined) + ofi * as.factor(bin_lag_refined)"
  )

  spec_labels <- c(
    "contemp_bin", "contemp_bin + lag_bin",
    rep("contemp_refined_bin + lag_refined_bin", 2)
  )

  # introduce direct controls
  out_all <- data.table() # save all results here
  for (spec_idx in 1:length(spec_formulas)) {
    # print(spec_idx)
    # spec_idx <- 2
    ff <- paste0("ret ~ ", spec_formulas[spec_idx], " + ", controls)
    out <- p.fama_macbeth(data, as.formula(ff), compare_coefs = FALSE)
    out[, spec_idx := spec_idx]
    out_all <- rbind(out_all, out)
  }

  out_all[grepl("ofi_is", var)]

  # name the control specs
  tmp <- data.table(
    spec_idx = 1:length(spec_formulas),
    spec_label = spec_labels
  )
  out_all <- merge(out_all, tmp, by = "spec_idx", all.x = T)
  out_all[, type := this_type]

  return(out_all)
}

tic("regression")
data_list <- split(data_reg, by = "type")
out <- rbindlist(mclapply(data_list, p.process_one_type, mc.cores = nc))
toc()

# save results
to_file <- "tmp/reg_results/sign_with_controls.RDS"
dir.create(dirname(to_file), recursive = T, showWarnings = F)
saveRDS(out, to_file)

# # --- SANITY: take a look at results

# out_bk <- readRDS("tmp/reg_results/sign_with_controls.RDS")
# out <- copy(out_bk)
# tmp <- unique(out[, .(var)])[grepl("ofi_is", var)][, var_idx := .I][, var_lab := paste0(var_idx, "_", var)]
# out <- merge(out, tmp, by = "var")
# out[, coef_round := round(coef, 2)]
# out[, tstat_round := round(coef / se, 2)]
# rm(tmp)

# dcast(out[var == "ofi_is_same_sign"], type ~ spec_idx, value.var = "coef_round")
# dcast(out[var == "ofi_is_same_sign"], type ~ spec_idx, value.var = "tstat_round")

# tt <- out[var != "ofi_is_same_sign"]
# dcast(tt, var_lab ~ type + spec_idx, value.var = "coef_round")
# dcast(tt, var_lab ~ type + spec_idx, value.var = "tstat_round")
