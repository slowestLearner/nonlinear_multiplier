# ---- in regression, check if there appears to be differential reversals?
library(this.path)
setwd(this.path::this.dir())
source("../utilities/runmefirst.R")
source("../utilities/regressions.R")

# regression table
data_all <- readRDS("../tmp/raw_data/reg_inputs/reg_table_static.RDS")[, ret := NULL]

# Merge with future returns of different horizons
ret_data <- readRDS("../../../../data/stocks/prices/quarterly_return.RDS")[, .(yyyymm, permno, logret = log(1 + ret))] %>%
  na.omit() %>%
  mutate(idx = frank(yyyymm, ties.method = "dense")) %>%
  setDT()

# get cumulative future log returns
for (i in 1:12) {
  ret_data <- merge(ret_data, ret_data[, .(idx = idx - i, permno, xx = logret)], by = c("idx", "permno"), all.x = T)
  setnames(ret_data, "xx", paste0("ret", i))
}
for (i in 2:12) {
  setnames(ret_data, paste0("ret", c(i - 1, i)), c("xx", "yy"))
  ret_data[, yy := xx + yy]
  setnames(ret_data, c("xx", "yy"), paste0("ret", c(i - 1, i)))
}
ret_data[, c("logret", "idx") := NULL]

# reshape to long format, choose a few horizons
ret_data <- melt(ret_data, id.vars = c("yyyymm", "permno")) %>% setDT()
names(ret_data)[3:4] <- c("ret_hor", "ret")
ret_data[, ret_hor := as.integer(gsub("ret", "", ret_hor))]
ret_data <- ret_data %>%
  filter(ret_hor %in% c(1, 2, 4, 8, 12)) %>%
  na.omit() %>%
  setDT()

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

data_all <- merge(data_all, ret_data, by = c("yyyymm", "permno"), allow.cartesian = T)
rm(ret_data)

# ###########################################################
# Fama-MacBeth
# ###########################################################

tic("regressions")

# control variables for different specifications
control_formulas <- c(
  "1",
  paste0(c(controls_char), collapse = "+"),
  paste0(c(controls_char, controls_liq), collapse = "+")
)

# worker function
# TODO: remove nonlinear
p.process_one_type <- function(data) {
  this_type <- data[1, type]
  this_ret_hor <- data[1, ret_hor]

  out_all <- data.table() # save all results here

  # introduce direct controls
  for (spec_idx in 1:length(control_formulas)) {
    ff_no_ofi <- paste0("ret ~ ", control_formulas[spec_idx])
    ff <- paste0("ret ~ ", control_formulas[spec_idx], " + ofi + ofi_bin2 + ofi_bin3")
    if (this_type == "BMI") {
      ff <- paste0(ff, "+", paste0(controls_bmi, collapse = "+"))
      ff_no_ofi <- paste0(ff_no_ofi, "+", paste0(controls_bmi, collapse = "+"))
    }
    out <- p.fama_macbeth(data, ff, compare_coefs = FALSE)
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

    out <- p.fama_macbeth(data, ff, compare_coefs = FALSE)
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
  out_all[, ret_hor := this_ret_hor]

  return(out_all)
}

# on a computer with 6 cores, took around 4-5 mins
plan(multisession, workers = detectCores() - 2)
out_list <- future_lapply(
  split(data_all, by = c("type", "ret_hor")), p.process_one_type,
  future.packages = "data.table", future.seed = 123
)
out_stdev <- rbindlist(out_list)
plan(sequential)
toc()

to_dir <- "../tmp/price_impact/regression_contemp/reversals/"
dir.create(to_dir, recursive = T, showWarnings = F)
saveRDS(out_stdev, paste0(to_dir, "fm_stdev.RDS"))

# # --- SANITY CHECK: close enough

# new <- readRDS("../tmp/price_impact/regression_contemp/reversals/fm_stdev.RDS")
# old <- readRDS("../tmp/price_impact/regression_contemp_todel/reversals/fm_stdev.RDS")
# dim(new) == dim(old)

# out <- merge(new, old, by = c("spec_idx", "var", "type", "ret_hor"))
# out[, cor(coef.x, coef.y, use = "complete.obs"), type]
# out[, cor(se.x, se.y, use = "complete.obs"), type]
# out[, max(obs.x) / max(obs.y), type]
