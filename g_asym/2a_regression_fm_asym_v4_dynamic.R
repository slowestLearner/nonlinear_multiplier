# Dynamic regression of returns on demand -- rebuilt sign-bin v4
#
# PURPOSE
#   Dynamic Fama-MacBeth regression with positive- and negative-lagged-shock
#   bin slopes estimated in one monthly cross-section:
#     ret ~ ofi_bin1_pos + ofi_bin1_neg + ... + ofi_bin3_neg + controls
#
#   This mirrors code/R/c_dynamic_results/2_regression_dynamic_fm.R:
#   dynamic bins are based on cumofi_1,...,cumofi_4; the regressor inside each
#   dynamic bin is contemporaneous ofi; controls and progressive ofi-control
#   interactions follow the dynamic script.
#
# BIN / SIGN CONVENTION
#   Positive/negative status and bin thresholds are based on lagged shock
#   magnitude, cumofi_lag, not contemporaneous ofi. Positive bins use
#   sd(cumofi_lag | cumofi_lag > 0) within type-horizon-month; negative bins
#   use sd(cumofi_lag | cumofi_lag < 0) within type-horizon-month. The
#   regressor value inside each bin remains contemporaneous ofi.
#
# OUTPUTS
#   Unbalanced:
#     fm_dynamic_stdev_posneg_rebuilt_v4.RDS/.dta
#   Balanced:
#     fm_dynamic_stdev_posneg_rebuilt_v4_balanced.RDS/.dta
#     fm_dynamic_stdev_posneg_rebuilt_v4_balanced_dropped_type_months.RDS/.dta

library(this.path)
setwd(this.path::this.dir())
source("../R/utilities/runmefirst.R")
library(haven)


# ===========================================================================
# SECTION 1: Helpers
# ===========================================================================

SIGN_BIN_VARS <- c(
  "ofi_bin1_pos", "ofi_bin1_neg",
  "ofi_bin2_pos", "ofi_bin2_neg",
  "ofi_bin3_pos", "ofi_bin3_neg"
)

DIFF_SPECS <- list(
  pos_2_1 = c("ofi_bin2_pos", "ofi_bin1_pos"),
  pos_3_2 = c("ofi_bin3_pos", "ofi_bin2_pos"),
  pos_3_1 = c("ofi_bin3_pos", "ofi_bin1_pos"),
  neg_2_1 = c("ofi_bin2_neg", "ofi_bin1_neg"),
  neg_3_2 = c("ofi_bin3_neg", "ofi_bin2_neg"),
  neg_3_1 = c("ofi_bin3_neg", "ofi_bin1_neg")
)

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
    for (nm in names(DIFF_SPECS)) {
      lhs <- DIFF_SPECS[[nm]][1]
      rhs <- DIFF_SPECS[[nm]][2]
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
  coef_data[, demand_type := data[1, demand_type]]
  coef_data[, hor := data[1, hor]]
  coef_data[, nw_lag := this_hor]
  coef_data[]
}


# ===========================================================================
# SECTION 2: Load controls and dynamic regression input
# ===========================================================================

data <- readRDS("../R/tmp/raw_data/reg_inputs/reg_table_dynamic.RDS")

vars_id <- c("yyyymm", "permno", "type")
vars_reg <- c("ret", "ofi", "cumofi_1", "cumofi_2", "cumofi_3", "cumofi_4")
vars_controls <- setdiff(names(data), c(vars_id, vars_reg))

data_controls <- data[, c(vars_id, vars_controls), with = FALSE]
setnames(data_controls, "type", "demand_type")

data_reg <- data[, c(vars_id, vars_reg), with = FALSE]
data_reg <- melt(
  data_reg,
  id.vars = c("yyyymm", "permno", "type", "ret", "ofi"),
  variable.name = "hor",
  value.name = "cumofi_lag"
)
data_reg <- data_reg[0 == rowSums(is.na(data_reg))]
data_reg[, hor := as.integer(gsub("cumofi_", "", hor))]
rm(data, vars_id, vars_reg, vars_controls)

data_reg[, demand_type := type]
data_reg[, type := paste0(demand_type, "_", hor, "lag")]

data_reg[, cumofi_sd_pos := sd(cumofi_lag[cumofi_lag > 0]), by = .(yyyymm, type)]
data_reg[, cumofi_sd_neg := sd(cumofi_lag[cumofi_lag < 0]), by = .(yyyymm, type)]

data_reg[, ofi_bin1_pos := ofi * (cumofi_lag > 0 & cumofi_lag < cumofi_sd_pos)]
data_reg[, ofi_bin2_pos := ofi * (cumofi_lag >= cumofi_sd_pos & cumofi_lag < 2 * cumofi_sd_pos)]
data_reg[, ofi_bin3_pos := ofi * (cumofi_lag >= 2 * cumofi_sd_pos)]

data_reg[, ofi_bin1_neg := ofi * (cumofi_lag < 0 & cumofi_lag > -cumofi_sd_neg)]
data_reg[, ofi_bin2_neg := ofi * (cumofi_lag <= -cumofi_sd_neg & cumofi_lag > -2 * cumofi_sd_neg)]
data_reg[, ofi_bin3_neg := ofi * (cumofi_lag <= -2 * cumofi_sd_neg)]

data_reg <- merge(data_reg, data_controls, by = c("yyyymm", "permno", "demand_type"), all.x = TRUE)
data_reg <- data_reg[0 == rowSums(is.na(data_reg))]
rm(data_controls)

cdata <- readRDS("../R/tmp/raw_data/controls/controls_classification.RDS")
controls_char <- cdata[control_type == "return-predictor", var]
controls_liq <- cdata[control_type == "liquidity", var]
controls_list <- c(controls_char, controls_liq)
rm(cdata)

control_formulas <- c(
  "1",
  paste0(controls_char, collapse = "+"),
  paste0(c(controls_char, controls_liq), collapse = "+")
)


# ===========================================================================
# SECTION 3: Diagnostics
# ===========================================================================

cat("Dynamic input rows after melt/control merge:", nrow(data_reg), "\n")
print(data_reg[, .(rows = .N, months = uniqueN(yyyymm)), by = .(demand_type, hor)][order(demand_type, hor)])

cat("\nDynamic sign-specific SD diagnostics:\n")
print(data_reg[, .(
  months = uniqueN(yyyymm),
  missing_pos_sd_months = uniqueN(yyyymm[is.na(cumofi_sd_pos)]),
  missing_neg_sd_months = uniqueN(yyyymm[is.na(cumofi_sd_neg)]),
  min_pos_sd = min(cumofi_sd_pos, na.rm = TRUE),
  min_neg_sd = min(cumofi_sd_neg, na.rm = TRUE)
), by = .(demand_type, hor)][order(demand_type, hor)])

cat("\nDynamic rebuilt sign/bin diagnostics:\n")
print(data_reg[cumofi_lag != 0, .(
  obs = .N,
  mean_cumofi_lag = mean(cumofi_lag),
  sd_cumofi_lag = sd(cumofi_lag),
  mean_abs_cumofi_lag = mean(abs(cumofi_lag)),
  mean_ofi = mean(ofi)
), by = .(
  demand_type,
  hor,
  shock = fifelse(cumofi_lag > 0, "pos", "neg"),
  bin = fcase(
    ofi_bin1_pos != 0 | ofi_bin1_neg != 0, 1L,
    ofi_bin2_pos != 0 | ofi_bin2_neg != 0, 2L,
    ofi_bin3_pos != 0 | ofi_bin3_neg != 0, 3L,
    default = NA_integer_
  )
)][order(demand_type, hor, shock, bin)])


# ===========================================================================
# SECTION 4: Regression driver
# ===========================================================================

p.process_one_type <- function(data) {
  this_type <- data[1, type]
  out_all <- data.table()

  for (spec_idx in seq_along(control_formulas)) {
    ff_no_ofi <- paste0("ret ~ ", control_formulas[spec_idx])
    ff <- paste0("ret ~ ", control_formulas[spec_idx], " + ", paste0(SIGN_BIN_VARS, collapse = " + "))

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

p.run_variant <- function(data_in, variant_name, balanced = FALSE) {
  data <- copy(data_in)

  if (balanced) {
    bin_pop <- data[, as.list(setNames(
      lapply(SIGN_BIN_VARS, function(v) sum(get(v) != 0, na.rm = TRUE)),
      SIGN_BIN_VARS
    )), by = .(type, demand_type, hor, yyyymm)]
    bin_pop[, all_bins_populated := Reduce(`&`, lapply(.SD, function(x) x > 0)),
            .SDcols = SIGN_BIN_VARS]
    dropped <- bin_pop[all_bins_populated == FALSE]

    cat("\nDynamic balanced-bin type-month filter:\n")
    print(bin_pop[, .(
      months_total = .N,
      months_kept = sum(all_bins_populated),
      months_dropped = sum(!all_bins_populated)
    ), by = .(type, demand_type, hor)][order(type)])

    saveRDS(dropped, paste0("fm_dynamic_stdev_posneg_rebuilt_v4_", variant_name, "_dropped_type_months.RDS"))
    write_dta(as.data.frame(dropped), paste0("fm_dynamic_stdev_posneg_rebuilt_v4_", variant_name, "_dropped_type_months.dta"))

    data <- merge(
      data,
      bin_pop[all_bins_populated == TRUE, .(type, yyyymm)],
      by = c("type", "yyyymm"),
      all = FALSE
    )
  }

  mc_cores <- max(1L, as.integer(nc), na.rm = TRUE)
  tic(paste("dynamic rebuilt pooled sign-bin", variant_name))
  out <- rbindlist(mclapply(split(data, by = "type"), function(x) {
    p.process_one_type(x)
  }, mc.cores = mc_cores))
  toc()

  saveRDS(out, paste0("fm_dynamic_stdev_posneg_rebuilt_v4_", variant_name, ".RDS"))
  write_dta(as.data.frame(out), paste0("fm_dynamic_stdev_posneg_rebuilt_v4_", variant_name, ".dta"))

  cat("\nDynamic rebuilt pooled sign-bin summary:", variant_name, "\n")
  print(out[
    spec_idx == 3 &
      var %in% c(SIGN_BIN_VARS,
                 paste0("ofi_bin2_pos - ofi_bin1_pos"),
                 paste0("ofi_bin3_pos - ofi_bin2_pos"),
                 paste0("ofi_bin3_pos - ofi_bin1_pos"),
                 paste0("ofi_bin2_neg - ofi_bin1_neg"),
                 paste0("ofi_bin3_neg - ofi_bin2_neg"),
                 paste0("ofi_bin3_neg - ofi_bin1_neg"))
  ][order(type, var)])

  invisible(out)
}

p.run_variant(data_reg, "unbalanced", balanced = FALSE)
p.run_variant(data_reg, "balanced", balanced = TRUE)
