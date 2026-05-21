# --- modified from 2_regression_dynamic_fm
library(this.path)
setwd(this.path::this.dir())
source("../../utilities/runmefirst.R")
source("../../utilities/regressions.R")

# load regression data
tic("Loading regression data")
data <- readRDS("../../tmp/raw_data/reg_inputs/reg_table_dynamic.RDS")[type %in% c("OFI", "FIT")]

# separate controls and main variables
vars_id <- c("yyyymm", "permno", "type")
vars_reg <- c("ret", "ofi", "cumofi_1", "cumofi_2", "cumofi_3", "cumofi_4")
vars_controls <- setdiff(names(data), c(vars_id, vars_reg))

# put regression-related data into long format
data_controls <- data[, c(vars_id, vars_controls), with = F]
setnames(data_controls, "type", "demand_type")

# get various lags
data_reg <- data[, .(ret = last(ret), ofi = last(ofi)), .(yyyymm, permno, type)]
tmp <- readRDS("tmp/raw_files/concentration_of_lagged_demand.RDS")
setnames(tmp, "lag", "hor")
data_reg <- merge(data_reg, tmp, by = c("yyyymm", "permno", "type"), allow.cartesian = T)
rm(tmp)

# change some variable names
data_reg[, demand_type := type]
data_reg[, type := paste0(demand_type, "_", hor, "lag")]

# data_reg[, sd_ofi := sd(cumofi_lag), .(yyyymm, type)]
# data_reg[, bin := 1]
# data_reg[abs(cumofi_lag) > sd_ofi, bin := 2]
# data_reg[abs(cumofi_lag) > 2 * sd_ofi, bin := 3]
# data_reg[, ofi_bin1 := ifelse(bin == 1, ofi, 0)]
# data_reg[, ofi_bin2 := ifelse(bin == 2, ofi, 0)]
# data_reg[, ofi_bin3 := ifelse(bin == 3, ofi, 0)]
# data_reg[, sd_ofi := NULL]
# data_reg[, ofi_absofi := ofi * abs(cumofi_lag)]

# further separate by dispersio into bins
data_reg[, bin_disp := ntile(disp, 3), .(yyyymm, type, hor, disp_type)]

# sort within bin
data_reg[, sd_ofi := sd(cumofi_lag), .(yyyymm, type, disp_type, bin_disp)]
data_reg[, bin := 1]
data_reg[abs(cumofi_lag) > sd_ofi, bin := 2]
data_reg[abs(cumofi_lag) > 2 * sd_ofi, bin := 3]
data_reg[, ofi_bin1 := ifelse(bin == 1, ofi, 0)]
data_reg[, ofi_bin2 := ifelse(bin == 2, ofi, 0)]
data_reg[, ofi_bin3 := ifelse(bin == 3, ofi, 0)]
data_reg[, sd_ofi := NULL]
data_reg[, ofi_absofi := ofi * abs(cumofi_lag)]

# get control variable names
cdata <- readRDS("../../tmp/raw_data/controls/controls_classification.RDS")
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
# data <- data_list[[1]]
p.process_one_type <- function(data) {
  data <- merge(data, data_controls, by = c("yyyymm", "permno", "demand_type"), all.x = T) %>%
    na.omit() %>%
    setDT()
  this_type <- data[1, type]
  this_disp_type <- data[1, disp_type]
  this_bin_disp <- data[1, bin_disp]

  out_all <- data.table() # save all results here

  # introduce direct controls
  for (spec_idx in 1:length(control_formulas)) {
    ff <- paste0("ret ~ ", control_formulas[spec_idx], " + ofi_bin1 + ofi_bin2 + ofi_bin3")
    out <- p.fama_macbeth(data, ff, compare_coefs = TRUE)
    out[, spec_idx := spec_idx]
    out_all <- rbind(out_all, out)
  }

  # # then add interactions with ofi
  # for (this_v in c(controls_char, controls_liq)) {
  #   setnames(data, this_v, "xx")
  #   data[, yy := xx * ofi]
  #   setnames(data, c("xx", "yy"), c(this_v, paste0("ofi_", this_v)))
  #   ff <- paste0(ff, " + ofi_", this_v)
  #   spec_idx <- spec_idx + 1

  #   out <- p.fama_macbeth(data, ff, compare_coefs = ifelse(reg_spec == "nonlinear", FALSE, TRUE))
  #   out[, spec_idx := spec_idx]
  #   out_all <- rbind(out_all, out)
  # }

  # name the control variables being added
  tmp <- data.table(
    spec_idx = 1:3,
    var_added = c("none_init", "controls_char", "controls_char+controls_liq")
  )
  out_all <- merge(out_all, tmp, by = "spec_idx", all.x = T)
  out_all[, type := this_type]
  out_all[, disp_type := this_disp_type]
  out_all[, bin_disp := this_bin_disp]

  return(out_all)
}

# it turns out that we have to resort by subsample
data_reg[, sd_ofi := sd(ofi), .(yyyymm, type, hor)]
data_reg[, bin := 1]
data_reg[abs(ofi) > sd_ofi, bin := 2]
data_reg[abs(ofi) > 2 * sd_ofi, bin := 3]
data_reg[, ofi_bin1 := ifelse(bin == 1, ofi, 0)]
data_reg[, ofi_bin2 := ifelse(bin == 2, ofi, 0)]
data_reg[, ofi_bin3 := ifelse(bin == 3, ofi, 0)]
data_reg[, sd_ofi := NULL]

tic("regression")
data_list <- split(data_reg, by = c("type", "disp_type", "bin_disp"))
out_all <- rbindlist(mclapply(data_list, p.process_one_type, mc.cores = nc))
to_file <- "tmp/reg_by_dispersion/dynamic_coefs_by_3bins.RDS"
dir.create(dirname(to_file), recursive = T, showWarnings = F)
saveRDS(out_all, to_file)
toc()

# # # --- SANITY CHECK: plot results

# out <- readRDS("tmp/reg_by_dispersion/dynamic_coefs_by_3bins.RDS")[spec_idx == 3][, spec_idx := NULL]

# # select vars
# tmp <- unique(out[grepl("ofi_bin", var)][, .(var)])
# tmp <- tmp[c(1:3, 4, 7, 5)][, var_lab := paste0(.I, "_", var)]
# out <- merge(out, tmp, by = "var")
# out_all <- copy(out)

# # only works with ofi_sd
# this_disp_type <- "ofi_sd"
# # this_disp_type <- "hhi_pow2"
# # this_lag <- 8

# # out <- copy(out_all)[disp_type == this_disp_type & nw_lag == this_lag]
# # out[, coef_round := round(coef, 2)]
# # out[, tstat_round := round(coef / se, 2)]

# # dcast(out, var_lab ~ type + bin_disp, value.var = "coef_round")
# # dcast(out, var_lab ~ type + bin_disp, value.var = "tstat_round")
