# helpers.R — shared constants and FM helper for g_asym sign-bin analysis.
# Sourced by: 2a_regression_fm_asym_v4.R, 2a_regression_fm_asym_v4_dynamic.R,
#             test_v4_delta_posneg.R

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
  coef_data[, nw_lag := this_hor]
  coef_data[]
}
