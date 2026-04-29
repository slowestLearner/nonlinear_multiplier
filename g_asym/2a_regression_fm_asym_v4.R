# Contemporaneous regression of returns on demand -- rebuilt sign-bin v4
#
# PURPOSE
#   Fama-MacBeth regression with positive- and negative-shock bin slopes
#   estimated in one monthly cross-section:
#     ret ~ ofi_bin1_pos + ofi_bin1_neg + ... + ofi_bin3_neg + controls
#
#   This keeps the architecture and controls of 2a_regression_fm_asym_v2.R:
#   BMI, FIT, and OFI are estimated separately; spec_idx 1-3 are no controls,
#   return-predictor controls, and return-predictor + liquidity controls; the
#   progressive interaction specs are also run; and BMI gets its BMI-specific
#   controls in every specification.
#
# INTERPRETATION
#   v2 estimates positive and negative subsamples separately using bins based
#   on the pooled OFI standard deviation. v4 keeps both signs in the same
#   monthly cross-section and estimates six sign-specific bin slopes with
#   common controls, but rebuilds the bins separately within positive and
#   negative shocks. The positive-side thresholds are based on sd(ofi | ofi>0)
#   within type-month; the negative-side thresholds are based on sd(ofi | ofi<0)
#   within type-month.

library(this.path)
setwd(this.path::this.dir())
source("../R/utilities/runmefirst.R")
library(haven)


# ===========================================================================
# SECTION 1: Control variable names
# ===========================================================================

cdata         <- readRDS("../R/tmp/raw_data/controls/controls_classification.RDS")
controls_char <- cdata[control_type == "return-predictor", var]
controls_liq  <- cdata[control_type == "liquidity", var]
controls_list <- c(controls_char, controls_liq)
rm(cdata)

tmp          <- readRDS("../R/tmp/raw_data/controls/controls_for_BMI.RDS")
controls_bmi <- setdiff(names(tmp), c("yyyymm", "permno"))
rm(tmp)

control_formulas <- c(
  "1",
  paste0(controls_char, collapse = "+"),
  paste0(c(controls_char, controls_liq), collapse = "+")
)

SIGN_BIN_VARS <- c(
  "ofi_bin1_pos", "ofi_bin1_neg",
  "ofi_bin2_pos", "ofi_bin2_neg",
  "ofi_bin3_pos", "ofi_bin3_neg"
)


# ===========================================================================
# SECTION 2: Fama-MacBeth helper for sign-bin variables
# ===========================================================================

p.fama_macbeth_signbins <- function(data, ff, compare_coefs = FALSE) {
  p.get_one_period <- function(this_ym) {
    ols <- lm(ff, data[yyyymm == this_ym])
    dep_var_name <- all.vars(as.formula(ff))[1]
    data.table(
      yyyymm = this_ym,
      var = names(coef(ols)),
      coef = as.numeric(coef(ols)),
      r2 = var(ols$fitted.values) / var(ols$model[[dep_var_name]]),
      obs_month = nobs(ols)
    )
  }

  out <- rbindlist(lapply(sort(unique(data[, yyyymm])), p.get_one_period))
  r2_data <- unique(out[, .(yyyymm, r2)])[, .(r2 = mean(r2))]
  obs_data <- unique(out[, .(yyyymm, obs_month)])[, .(obs = sum(obs_month))]
  out[, c("r2", "obs_month") := NULL]

  if (compare_coefs) {
    diff_specs <- list(
      pos_2_1 = c("ofi_bin2_pos", "ofi_bin1_pos"),
      pos_3_2 = c("ofi_bin3_pos", "ofi_bin2_pos"),
      pos_3_1 = c("ofi_bin3_pos", "ofi_bin1_pos"),
      neg_2_1 = c("ofi_bin2_neg", "ofi_bin1_neg"),
      neg_3_2 = c("ofi_bin3_neg", "ofi_bin2_neg"),
      neg_3_1 = c("ofi_bin3_neg", "ofi_bin1_neg")
    )
    for (nm in names(diff_specs)) {
      lhs <- diff_specs[[nm]][1]
      rhs <- diff_specs[[nm]][2]
      diff_data <- merge(
        out[var == lhs, .(yyyymm, coef_lhs = coef)],
        out[var == rhs, .(yyyymm, coef_rhs = coef)],
        by = "yyyymm",
        all = FALSE
      )
      out <- rbind(out, data.table(
        yyyymm = diff_data[, yyyymm],
        var = paste0(lhs, " - ", rhs),
        coef = diff_data[, coef_lhs - coef_rhs]
      ), fill = TRUE)
    }
  }

  this_hor <- if ("hor" %in% names(data)) data[1, hor] else 1

  coef_data <- data.table()
  for (this_var in unique(out[, var])) {
    this_series <- out[var == this_var & !is.na(coef), .(coef)]
    mm <- lm(coef ~ 1, this_series)
    coef_data <- rbind(coef_data, data.table(
      var = this_var,
      coef = mm$coef[1],
      se = sqrt(vcov(mm)[1, 1]),
      se_nw = sqrt(NeweyWest(mm, this_hor)[1, 1]),
      n_months = nrow(this_series)
    ))
  }
  coef_data[, r2 := r2_data[, r2]]
  coef_data[, obs := obs_data[, obs]]
  coef_data[, type := data[1, type]]
  coef_data[, nw_lag := this_hor]
  coef_data[]
}


# ===========================================================================
# SECTION 3: Process one demand type through v2-style specifications
# ===========================================================================

p.process_one_type <- function(data) {
  this_type <- data[1, type]
  out_all <- data.table()

  for (spec_idx in 1:length(control_formulas)) {
    ff_no_ofi <- paste0("ret ~ ", control_formulas[spec_idx])
    ff <- paste0("ret ~ ", control_formulas[spec_idx], " + ",
                 paste0(SIGN_BIN_VARS, collapse = " + "))

    if (this_type == "BMI") {
      ff <- paste0(ff, "+", paste0(controls_bmi, collapse = "+"))
      ff_no_ofi <- paste0(ff_no_ofi, "+", paste0(controls_bmi, collapse = "+"))
    }

    out <- p.fama_macbeth_signbins(data, ff, compare_coefs = TRUE)
    out_no_ofi <- p.fama_macbeth_signbins(data, ff_no_ofi, compare_coefs = FALSE)
    out[, r2_no_ofi := out_no_ofi[1, r2]]
    out[, spec_idx := spec_idx]
    out_all <- rbind(out_all, out)
  }

  for (this_v in c(controls_char, controls_liq)) {
    setnames(data, this_v, "xx")
    data[, yy := xx * ofi]
    setnames(data, c("xx", "yy"), c(this_v, paste0("ofi_", this_v)))
    ff <- paste0(ff, " + ofi_", this_v)
    ff_no_ofi <- paste0(ff_no_ofi, " + ofi_", this_v)
    spec_idx <- spec_idx + 1

    out <- p.fama_macbeth_signbins(data, ff, compare_coefs = TRUE)
    out_no_ofi <- p.fama_macbeth_signbins(data, ff_no_ofi, compare_coefs = FALSE)
    out[, r2_no_ofi := out_no_ofi[1, r2]]
    out[, spec_idx := spec_idx]
    out_all <- rbind(out_all, out)
  }

  spec_labels <- data.table(
    spec_idx = 1:(length(controls_list) + 3),
    var_added = c("none_init", "controls_char", "controls_char+controls_liq", controls_list),
    var_type = c(rep("", 3),
                 rep("return-predicting chars", length(controls_char)),
                 rep("liquidity", length(controls_liq)))
  )
  out_all <- merge(out_all, spec_labels, by = "spec_idx", all.x = TRUE)
  out_all[, type := this_type]
  out_all[]
}


# ===========================================================================
# SECTION 4: Load data, rebuild sign bins, and run all three demand types
# ===========================================================================

data_all <- readRDS("../R/tmp/raw_data/reg_inputs/reg_table_static.RDS")
data_all <- data_all[type %in% c("BMI", "FIT", "OFI")]

cat("Input rows after table-type filter:", nrow(data_all), "\n")
print(data_all[, .N, by = type][order(type)])

data_all[, ofi_sd_pos := sd(ofi[ofi > 0]), by = .(type, yyyymm)]
data_all[, ofi_sd_neg := sd(ofi[ofi < 0]), by = .(type, yyyymm)]

cat("\nSign-specific SD diagnostics:\n")
print(data_all[, .(
  months = uniqueN(yyyymm),
  missing_pos_sd_months = uniqueN(yyyymm[is.na(ofi_sd_pos)]),
  missing_neg_sd_months = uniqueN(yyyymm[is.na(ofi_sd_neg)]),
  min_pos_sd = min(ofi_sd_pos, na.rm = TRUE),
  min_neg_sd = min(ofi_sd_neg, na.rm = TRUE)
), by = type][order(type)])

data_all[, ofi_bin1_pos := ofi * (ofi > 0 & ofi < ofi_sd_pos)]
data_all[, ofi_bin2_pos := ofi * (ofi >= ofi_sd_pos & ofi < 2 * ofi_sd_pos)]
data_all[, ofi_bin3_pos := ofi * (ofi >= 2 * ofi_sd_pos)]

data_all[, ofi_bin1_neg := ofi * (ofi < 0 & ofi > -ofi_sd_neg)]
data_all[, ofi_bin2_neg := ofi * (ofi <= -ofi_sd_neg & ofi > -2 * ofi_sd_neg)]
data_all[, ofi_bin3_neg := ofi * (ofi <= -2 * ofi_sd_neg)]

cat("\nSign/bin diagnostics:\n")
print(data_all[ofi != 0, .(
  obs = .N,
  mean_ofi = mean(ofi),
  sd_ofi = sd(ofi),
  mean_abs_ofi = mean(abs(ofi))
), by = .(
  type,
  shock = fifelse(ofi > 0, "pos", "neg"),
  bin = fcase(
    ofi_bin1_pos != 0 | ofi_bin1_neg != 0, 1L,
    ofi_bin2_pos != 0 | ofi_bin2_neg != 0, 2L,
    ofi_bin3_pos != 0 | ofi_bin3_neg != 0, 3L,
    default = NA_integer_
  )
)][order(type, shock, bin)])

mc_cores <- max(1L, as.integer(nc), na.rm = TRUE)

tic("static fm rebuilt pooled sign-bin v4")
out_stdev <- rbindlist(mclapply(split(data_all, by = "type"), function(x) {
  p.process_one_type(x)
}, mc.cores = mc_cores))
toc()

saveRDS(out_stdev, "fm_stdev_posneg_rebuilt_v4.RDS")
write_dta(as.data.frame(out_stdev), "fm_stdev_posneg_rebuilt_v4.dta")

monthly_like <- out_stdev[
  spec_idx %in% 1:3 &
    var %in% c(SIGN_BIN_VARS,
               paste0("ofi_bin2_pos - ofi_bin1_pos"),
               paste0("ofi_bin3_pos - ofi_bin2_pos"),
               paste0("ofi_bin3_pos - ofi_bin1_pos"),
               paste0("ofi_bin2_neg - ofi_bin1_neg"),
               paste0("ofi_bin3_neg - ofi_bin2_neg"),
               paste0("ofi_bin3_neg - ofi_bin1_neg"))
]

cat("\nRebuilt pooled sign-bin v4 summary, spec_idx 1-3:\n")
print(monthly_like[order(type, spec_idx, var)])
