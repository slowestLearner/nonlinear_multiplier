# Contemporaneous regression of returns on demand: panel regressions
source("utilities/runmefirst.R")


# load regression data
data_all <- readRDS("tmp/raw_data/reg_inputs/reg_table_static.RDS")

# get control variable names
cdata <- readRDS("tmp/raw_data/controls/controls_classification.RDS")
controls_char <- cdata[control_type == "return-predictor", var]
controls_liq <- cdata[control_type == "liquidity", var]
controls_list <- c(controls_char, controls_liq)
rm(cdata)

# get names for bmi controls
tmp <- readRDS("tmp/raw_data/controls/controls_for_BMI.RDS")
controls_bmi <- setdiff(names(tmp), c("yyyymm", "permno"))
rm(tmp)


# ###########################################################
# Utility Functions ---- TODO: move to separate file
# ###########################################################

# utility function: panel regression
p.panel_regression <- function(data, ff, compare_coefs = FALSE) {
  ols <- feols(as.formula(paste0(ff, " | yyyymm")), data, cluster = c("yyyymm", "permno"))
  out <- data.table(
    var = names(coef(ols)),
    coef = coef(ols),
    se = sqrt(diag(vcov(ols))),
    obs = ols$nobs, r2 = r2(ols)["ar2"]
  )

  # also report coef differences
  if (compare_coefs == TRUE) {
    var_indices <- names(coef(ols)) %in% paste0("ofi_bin", 1:3)

    cc <- matrix(coef(ols)[var_indices])
    C <- vcov(ols)[var_indices, var_indices]

    b_12 <- matrix(c(-1, 1, 0))
    b_23 <- matrix(c(0, -1, 1))
    b_13 <- matrix(c(-1, 0, 1))

    out <- rbind(out, data.table(
      var = c("ofi_bin2 - ofi_bin1", "ofi_bin3 - ofi_bin2", "ofi_bi3 - ofi_bin1"),
      coef = c(
        (t(b_12) %*% cc)[1],
        (t(b_23) %*% cc)[1],
        (t(b_13) %*% cc)[1]
      ),
      se = c(
        sqrt((t(b_12) %*% C %*% b_12)[1]),
        sqrt((t(b_23) %*% C %*% b_23)[1]),
        sqrt((t(b_13) %*% C %*% b_13)[1])
      ),
      obs = ols$nobs, r2 = r2(ols)["ar2"]
    ))
  }

  return(out)
}


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
p.process_one_type <- function(data, reg_spec = "nonlinear") {
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
    out <- p.panel_regression(data, ff, compare_coefs = ifelse(reg_spec == "nonlinear", FALSE, TRUE))
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

    out <- p.panel_regression(data, ff, compare_coefs = ifelse(reg_spec == "nonlinear", FALSE, TRUE))
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

  return(out_all)
}

# nonlinear specification---
tic()
out_nonlinear <- rbindlist(mclapply(split(data_all, by = "type"), function(x) {
  p.process_one_type(x, reg_spec = "nonlinear")
}, mc.cores = nc))
toc()

dir.create("tmp/price_impact/regression_contemp/", recursive = T, showWarnings = F)
saveRDS(out_nonlinear, "tmp/price_impact/regression_contemp/panel_nonlinear.RDS")

# stdev-based specification---
tic()
out_stdev <- rbindlist(mclapply(split(data_all, by = "type"), function(x) {
  p.process_one_type(x, reg_spec = "stdev")
}, mc.cores = nc))
toc()

dir.create("tmp/price_impact/regression_contemp/", recursive = T, showWarnings = F)
saveRDS(out_stdev, "tmp/price_impact/regression_contemp/panel_stdev.RDS")


# # === SANITY check


# # --- nonlinear
# new <- readRDS("tmp/price_impact/regression_contemp/panel_nonlinear.RDS")
# old <- readRDS("../20250117_quarterly/tmp/price_impact/regression_contemp/full_sample_panel.RDS")
# old <- old[type != "OFI"]
# old[type == "OFI_resid", type := "OFI"]

# # get overlapping part
# new <- new[spec_idx %in% 1:3]
# vv <- intersect(names(new), names(old))
# new <- new[, ..vv]
# old <- old[, ..vv]
# rm(vv)
# stopifnot(dim(new) == dim(old))

# # convert coef units
# old[var == "ofi_absofi", coef := coef * 100]
# old[var == "ofi_absofi", se := se * 100]

# # compare
# compare <- merge(new, old, by = c("type", "var", "spec_idx"), all = T)
# stopifnot(nrow(compare) == nrow(new))
# rm(new, old)

# compare[, cor(coef.x, coef.y)]
# compare[, cor(se.x, se.y)]
# compare[, cor(obs.x, obs.y)]
# compare[, cor(r2.x, r2.y)]

# tt = readRDS('../20250117_quarterly/tmp/price_impact/archive/multiplier_by_shock_size_quarterly/panel.RDS')
