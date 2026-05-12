# Formal v4 tests of positive-vs-negative concavity differences.
#
# For each demand type, specification, and (for dynamic results) horizon, this
# script tests whether the pairwise bin difference for positive shocks equals
# the analogous pairwise bin difference for negative shocks:
#
#   H0: Delta_{j-k}(+) = Delta_{j-k}(-)
#
# where Delta_{j-k}(+) = M_j(+) - M_k(+) and Delta_{j-k}(-) = M_j(-) - M_k(-).
# The test statistic is formed from the monthly series
#   [Delta_{j-k}(+) - Delta_{j-k}(-)]_t
# using a standard Fama-MacBeth mean and standard error across months. Omitted
# monthly coefficients remain missing; months enter a test only when all four
# required bin coefficients are available.

library(this.path)
setwd(this.path::this.dir())
source("../R/utilities/runmefirst.R")
source("helpers.R")
library(haven)


# ===========================================================================
# Shared helpers
# ===========================================================================

DELTA_SPECS <- list(
  "2-1" = list(pos_lhs = "ofi_bin2_pos", pos_rhs = "ofi_bin1_pos",
               neg_lhs = "ofi_bin2_neg", neg_rhs = "ofi_bin1_neg"),
  "3-2" = list(pos_lhs = "ofi_bin3_pos", pos_rhs = "ofi_bin2_pos",
               neg_lhs = "ofi_bin3_neg", neg_rhs = "ofi_bin2_neg"),
  "3-1" = list(pos_lhs = "ofi_bin3_pos", pos_rhs = "ofi_bin1_pos",
               neg_lhs = "ofi_bin3_neg", neg_rhs = "ofi_bin1_neg")
)

fm_mean <- function(x) {
  x <- x[!is.na(x)]
  data.table(
    coef = mean(x),
    se = sd(x) / sqrt(length(x)),
    t = mean(x) / (sd(x) / sqrt(length(x))),
    n_months = length(x)
  )
}

get_monthly_coefs <- function(data, ff) {
  rbindlist(lapply(sort(unique(data[, yyyymm])), function(this_ym) {
    ols <- lm(ff, data[yyyymm == this_ym])
    data.table(
      yyyymm = this_ym,
      var = names(coef(ols)),
      coef = as.numeric(coef(ols)),
      obs_month = nobs(ols)
    )
  }))
}

compute_delta_tests <- function(monthly, meta) {
  monthly_wide <- dcast(
    monthly[var %in% SIGN_BIN_VARS],
    yyyymm ~ var,
    value.var = "coef"
  )

  out <- rbindlist(lapply(names(DELTA_SPECS), function(delta_name) {
    spec <- DELTA_SPECS[[delta_name]]
    needed <- unlist(spec, use.names = FALSE)
    test_data <- monthly_wide[complete.cases(monthly_wide[, ..needed])]
    test_data[, delta_pos := get(spec$pos_lhs) - get(spec$pos_rhs)]
    test_data[, delta_neg := get(spec$neg_lhs) - get(spec$neg_rhs)]
    test_data[, delta_pos_minus_neg := delta_pos - delta_neg]

    ans <- fm_mean(test_data[, delta_pos_minus_neg])
    ans[, `:=`(
      delta = delta_name,
      mean_delta_pos = mean(test_data[, delta_pos]),
      mean_delta_neg = mean(test_data[, delta_neg])
    )]
    ans
  }), fill = TRUE)

  for (nm in names(meta)) out[, (nm) := meta[[nm]]]
  out[]
}

add_spec_labels <- function(out, controls_char, controls_liq) {
  controls_list <- c(controls_char, controls_liq)
  spec_labels <- data.table(
    spec_idx = 1:(length(controls_list) + 3),
    var_added = c("none_init", "controls_char", "controls_char+controls_liq", controls_list),
    var_type = c(rep("", 3),
                 rep("return-predicting chars", length(controls_char)),
                 rep("liquidity", length(controls_liq)))
  )
  merge(out, spec_labels, by = "spec_idx", all.x = TRUE)
}


# ===========================================================================
# Static v4 tests
# ===========================================================================

run_static_tests <- function() {
  cdata <- readRDS("../R/tmp/raw_data/controls/controls_classification.RDS")
  controls_char <- cdata[control_type == "return-predictor", var]
  controls_liq <- cdata[control_type == "liquidity", var]
  controls_list <- c(controls_char, controls_liq)
  rm(cdata)

  tmp <- readRDS("../R/tmp/raw_data/controls/controls_for_BMI.RDS")
  controls_bmi <- setdiff(names(tmp), c("yyyymm", "permno"))
  rm(tmp)

  control_formulas <- c(
    "1",
    paste0(controls_char, collapse = "+"),
    paste0(c(controls_char, controls_liq), collapse = "+")
  )

  data_all <- as.data.table(readRDS("../R/tmp/raw_data/reg_inputs/reg_table_static.RDS"))
  data_all <- data_all[type %in% c("BMI", "FIT", "OFI")]

  data_all[, ofi_sd_pos := sd(ofi[ofi > 0]), by = .(type, yyyymm)]
  data_all[, ofi_sd_neg := sd(ofi[ofi < 0]), by = .(type, yyyymm)]

  data_all[, ofi_bin1_pos := ofi * (ofi > 0 & ofi < ofi_sd_pos)]
  data_all[, ofi_bin2_pos := ofi * (ofi >= ofi_sd_pos & ofi < 2 * ofi_sd_pos)]
  data_all[, ofi_bin3_pos := ofi * (ofi >= 2 * ofi_sd_pos)]

  data_all[, ofi_bin1_neg := ofi * (ofi < 0 & ofi > -ofi_sd_neg)]
  data_all[, ofi_bin2_neg := ofi * (ofi <= -ofi_sd_neg & ofi > -2 * ofi_sd_neg)]
  data_all[, ofi_bin3_neg := ofi * (ofi <= -2 * ofi_sd_neg)]

  out <- data.table()
  for (this_type in c("BMI", "FIT", "OFI")) {
    d <- data_all[type == this_type]
    for (spec_idx in seq_along(control_formulas)) {
      ff <- paste0(
        "ret ~ ", control_formulas[spec_idx],
        " + ", paste0(SIGN_BIN_VARS, collapse = " + ")
      )
      if (this_type == "BMI") {
        ff <- paste0(ff, " + ", paste0(controls_bmi, collapse = " + "))
      }

      monthly <- get_monthly_coefs(d, ff)
      tests <- compute_delta_tests(monthly, list(
        sample = "static",
        type = this_type,
        spec_idx = spec_idx
      ))
      out <- rbind(out, tests, fill = TRUE)
    }
  }

  add_spec_labels(out, controls_char, controls_liq)
}


# ===========================================================================
# Dynamic v4 tests
# ===========================================================================

run_dynamic_tests <- function() {
  data <- as.data.table(readRDS("../R/tmp/raw_data/reg_inputs/reg_table_dynamic.RDS"))

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

  out <- data.table()
  for (this_type in sort(unique(data_reg[, type]))) {
    d <- data_reg[type == this_type]
    for (spec_idx in seq_along(control_formulas)) {
      ff <- paste0(
        "ret ~ ", control_formulas[spec_idx],
        " + ", paste0(SIGN_BIN_VARS, collapse = " + ")
      )

      monthly <- get_monthly_coefs(d, ff)
      tests <- compute_delta_tests(monthly, list(
        sample = "dynamic",
        type = this_type,
        demand_type = d[1, demand_type],
        hor = d[1, hor],
        spec_idx = spec_idx
      ))
      out <- rbind(out, tests, fill = TRUE)
    }
  }

  add_spec_labels(out, controls_char, controls_liq)
}


# ===========================================================================
# Run and save
# ===========================================================================

cat("Running static v4 positive-vs-negative delta tests...\n")
static_tests <- run_static_tests()
saveRDS(static_tests, "v4_delta_posneg_tests_static.RDS")
write_dta(as.data.frame(static_tests), "v4_delta_posneg_tests_static.dta")
fwrite(static_tests, "v4_delta_posneg_tests_static.csv")

cat("Running dynamic v4 positive-vs-negative delta tests...\n")
dynamic_tests <- run_dynamic_tests()
saveRDS(dynamic_tests, "v4_delta_posneg_tests_dynamic.RDS")
write_dta(as.data.frame(dynamic_tests), "v4_delta_posneg_tests_dynamic.dta")
fwrite(dynamic_tests, "v4_delta_posneg_tests_dynamic.csv")

cat("\nStatic v4 tests, most controlled specification (spec_idx = 3):\n")
print(static_tests[
  spec_idx == 3,
  .(type, delta, mean_delta_pos, mean_delta_neg, coef, se, t, n_months)
][order(type, delta)])

cat("\nDynamic v4 tests, most controlled specification (spec_idx = 3):\n")
print(dynamic_tests[
  spec_idx == 3,
  .(demand_type, hor, delta, mean_delta_pos, mean_delta_neg, coef, se, t, n_months)
][order(demand_type, hor, delta)])
