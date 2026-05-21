# --- Modified from 2a in main code. Separate by vol regime
library(this.path)
setwd(this.path::this.dir())
source("../../utilities/runmefirst.R")
source("../../utilities/regressions.R")
options(width = 120)
nc <- detectCores() - 2

# load regression data
tic("preparing data")
data_all <- readRDS("../../tmp/raw_data/reg_inputs/reg_table_static.RDS")[type %in% c("BMI", "OFI", "FIT")]

# get control variable names
cdata <- readRDS("../../tmp/raw_data/controls/controls_classification.RDS")
controls_char <- cdata[control_type == "return-predictor", var]
controls_liq <- cdata[control_type == "liquidity", var]
controls_list <- c(controls_char, controls_liq)
rm(cdata)

# get names for bmi controls
tmp <- readRDS("../../tmp/raw_data/controls/controls_for_BMI.RDS")
controls_bmi <- setdiff(names(tmp), c("yyyymm", "permno"))
rm(tmp)
toc()

# merge with regimes
n_bins <- 2
data_all[, yyyy := floor(yyyymm / 100)]

tmp <- readRDS("tmp/vols/mkt_and_stock_combined.RDS")
tmp <- melt(tmp, id.vars = "yyyymm", variable.name = "vol_type", value.name = "vol") %>%
  mutate(vol_type = as.character(vol_type)) %>%
  na.omit() %>%
  setDT()
tmp <- merge(tmp, unique(data_all[, .(yyyymm)]), by = "yyyymm")
tmp[, bin_vol := ntile(vol, n_bins), vol_type][, vol := NULL]
data_all <- merge(data_all, tmp, by = "yyyymm", allow.cartesian = TRUE)
rm(tmp)


# ###########################################################
# Fama-MacBeth
# ###########################################################


# worker function to estimate regression with one type of demand variable. reg_spec = "nonlinear" or "stdev"
# data <- data_list[[1]]
p.process_one_type <- function(data) {
  # parse
  this_type <- data[1, type]
  vol_bin <- data[1, bin_vol]
  vol_type <- data[1, vol_type]

  # control variables for different regression specifications
  control_formulas <- c(
    "1",
    paste0(c(controls_char), collapse = "+"),
    paste0(c(controls_char, controls_liq), collapse = "+")
  )

  out_all <- data.table() # save results here

  # introduce direct controls
  for (spec_idx in 1:length(control_formulas)) {
    # spec_idx <- 1
    ff_no_ofi <- paste0("ret ~ ", control_formulas[spec_idx])
    ff <- paste0("ret ~ ", control_formulas[spec_idx], " + ofi_bin1 + ofi_bin2 + ofi_bin3")
    if (this_type == "BMI") {
      ff <- paste0(ff, "+", paste0(controls_bmi, collapse = "+"))
      ff_no_ofi <- paste0(ff_no_ofi, "+", paste0(controls_bmi, collapse = "+"))
    }
    out <- p.fama_macbeth(data, ff, compare_coefs = TRUE)
    out_no_ofi <- p.fama_macbeth(data, ff_no_ofi, compare_coefs = FALSE)

    out[, r2_no_ofi := out_no_ofi[1, r2]]
    out[, spec_idx := spec_idx]
    out_all <- rbind(out_all, out)
  }


  # return
  out_all[, type := this_type]
  out_all[, vol_bin := vol_bin]
  out_all[, vol_type := vol_type]
  return(out_all)
}

tic("regression")
data_list <- split(data_all, by = c("type", "bin_vol", "vol_type"))
out <- rbindlist(mclapply(data_list, p.process_one_type, mc.cores = nc))
toc()

to_file <- "tmp/reg_results/by_vol_regime.RDS"
dir.create(dirname(to_file), recursive = T, showWarnings = F)
saveRDS(out, to_file)

# # --- take a look at results

# out_all <- readRDS("tmp/reg_results/by_vol_regime.RDS")[spec_idx == 3 & grepl("ofi", var)]

# # choose vars
# tmp <- data.table(var = c(paste0("ofi_bin", 1:3), "ofi_bin2 - ofi_bin1", "ofi_bin3 - ofi_bin2", "ofi_bin3 - ofi_bin1")) %>% mutate(var_lab = paste0(row_number(), "_", var))
# out_all <- merge(out_all, tmp, by = "var")
# rm(tmp)

# # rank the vol types
# tmp <- unique(out_all[, .(vol_type)])[order(vol_type)][, idx := c(1, 3, 2)][, vol_type_lab := paste0(idx, "_", vol_type)][, idx := NULL]
# out_all <- merge(out_all, tmp, by = "vol_type")
# rm(tmp)

# # round
# out_all[, coef_round := round(coef, 2)]
# out_all[, tstat_round := round(coef / se, 2)]

# # choose the type
# this_type <- "OFI"
# this_type <- "FIT"
# this_type <- "BMI"

# out <- copy(out_all[type == this_type])
# dcast(out, var_lab ~ vol_type_lab + vol_bin, value.var = "coef_round")
# dcast(out, var_lab ~ vol_type_lab + vol_bin, value.var = "tstat_round")
