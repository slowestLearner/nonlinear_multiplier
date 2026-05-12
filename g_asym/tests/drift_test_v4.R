# drift_test_v4.R
#
# PURPOSE
#   Guard four key v4 results against unintended changes during refactoring or
#   future modifications. Loads pre-computed RDS outputs and asserts that
#   protected coefficients and standard errors match expected values within
#   tolerance. Does NOT re-run any regression.
#
# PROTECTED RESULTS (confirmed 2026-05-07)
#   1. Static v4 — fm_stdev_posneg_rebuilt_v4.RDS
#        spec_idx = 3, var_added = "controls_char+controls_liq"
#        OFI bin1_pos / bin3_pos / (bin3_pos - bin1_pos) coef + se
#        FIT bin1_pos / bin3_pos / (bin3_pos - bin1_pos) coef + se
#        OFI bin1_neg / bin3_neg / (bin3_neg - bin1_neg) coef + se
#        FIT bin1_neg / bin3_neg / (bin3_neg - bin1_neg) coef + se
#
#   2. Dynamic v4 — fm_dynamic_stdev_posneg_rebuilt_v4_unbalanced.RDS
#        demand_type in {OFI, FIT}, hor = 4, spec_idx = 3
#        Same 12 coef/se pairs as above
#
# TOLERANCE
#   Coefficients and standard errors: abs(actual - expected) < 1e-4
#   This tolerance is tight relative to reporting precision (4 decimal places)
#   but generous relative to floating-point solver noise (<1e-7). It catches
#   any meaningful result drift while being immune to platform floating-point
#   differences in the last few ULPs.
#
# RUN
#   Rscript tests/drift_test_v4.R   (from the g_asym/ working directory)
#   Rscript drift_test_v4.R         (from the tests/ directory)

library(data.table)
library(this.path)

# Set working directory to the script's own location (tests/).
# RDS files live one level up in g_asym/.
setwd(this.path::this.dir())

# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

n_pass <- 0L

expect_near <- function(actual, expected, tol = 1e-4, label) {
  if (abs(actual - expected) >= tol) {
    stop(sprintf(
      "FAIL: %s\n  actual   = %.10f\n  expected = %.10f\n  diff     = %.2e  (tol = %.2e)",
      label, actual, expected, abs(actual - expected), tol
    ))
  }
  cat(sprintf("PASS: %s\n", label))
  n_pass <<- n_pass + 1L
}

# ---------------------------------------------------------------------------
# 1. Static v4 assertions
#    File : ../fm_stdev_posneg_rebuilt_v4.RDS
#    Filter: spec_idx == 3, var_added == "controls_char+controls_liq"
# ---------------------------------------------------------------------------

static <- readRDS("../fm_stdev_posneg_rebuilt_v4.RDS")

s <- static[spec_idx == 3L & var_added == "controls_char+controls_liq"]

# Helper: pull a single cell value, stop clearly if not found or ambiguous.
get_val <- function(dt, type_name, var_name, col) {
  rows <- dt[type == type_name & var == var_name]
  if (nrow(rows) == 0L)
    stop(sprintf("get_val: no row for type='%s' var='%s'", type_name, var_name))
  if (nrow(rows) > 1L)
    stop(sprintf("get_val: %d rows for type='%s' var='%s'", nrow(rows), type_name, var_name))
  rows[[col]]
}

## OFI — positive side
expect_near(get_val(s, "OFI", "ofi_bin1_pos", "coef"),  4.3361616, label = "Static OFI ofi_bin1_pos coef")
expect_near(get_val(s, "OFI", "ofi_bin1_pos", "se"),    0.1956185, label = "Static OFI ofi_bin1_pos se")
expect_near(get_val(s, "OFI", "ofi_bin3_pos", "coef"),  1.9151927, label = "Static OFI ofi_bin3_pos coef")
expect_near(get_val(s, "OFI", "ofi_bin3_pos", "se"),    0.1310172, label = "Static OFI ofi_bin3_pos se")
expect_near(get_val(s, "OFI", "ofi_bin3_pos - ofi_bin1_pos", "coef"), -2.4209689, label = "Static OFI ofi_bin3_pos - ofi_bin1_pos coef")
expect_near(get_val(s, "OFI", "ofi_bin3_pos - ofi_bin1_pos", "se"),   0.1444924,  label = "Static OFI ofi_bin3_pos - ofi_bin1_pos se")

## FIT — positive side
expect_near(get_val(s, "FIT", "ofi_bin1_pos", "coef"),  7.3154233, label = "Static FIT ofi_bin1_pos coef")
expect_near(get_val(s, "FIT", "ofi_bin1_pos", "se"),    1.0601004, label = "Static FIT ofi_bin1_pos se")
expect_near(get_val(s, "FIT", "ofi_bin3_pos", "coef"),  4.1383540, label = "Static FIT ofi_bin3_pos coef")
expect_near(get_val(s, "FIT", "ofi_bin3_pos", "se"),    0.4950520, label = "Static FIT ofi_bin3_pos se")
expect_near(get_val(s, "FIT", "ofi_bin3_pos - ofi_bin1_pos", "coef"), -3.1770693, label = "Static FIT ofi_bin3_pos - ofi_bin1_pos coef")
expect_near(get_val(s, "FIT", "ofi_bin3_pos - ofi_bin1_pos", "se"),   1.0144457,  label = "Static FIT ofi_bin3_pos - ofi_bin1_pos se")

## OFI — negative side
expect_near(get_val(s, "OFI", "ofi_bin1_neg", "coef"),  2.9279669, label = "Static OFI ofi_bin1_neg coef")
expect_near(get_val(s, "OFI", "ofi_bin1_neg", "se"),    0.1412138, label = "Static OFI ofi_bin1_neg se")
expect_near(get_val(s, "OFI", "ofi_bin3_neg", "coef"),  0.7037404, label = "Static OFI ofi_bin3_neg coef")
expect_near(get_val(s, "OFI", "ofi_bin3_neg", "se"),    0.0622382, label = "Static OFI ofi_bin3_neg se")
expect_near(get_val(s, "OFI", "ofi_bin3_neg - ofi_bin1_neg", "coef"), -2.2242265, label = "Static OFI ofi_bin3_neg - ofi_bin1_neg coef")
expect_near(get_val(s, "OFI", "ofi_bin3_neg - ofi_bin1_neg", "se"),   0.1195421,  label = "Static OFI ofi_bin3_neg - ofi_bin1_neg se")

## FIT — negative side
expect_near(get_val(s, "FIT", "ofi_bin1_neg", "coef"),  3.6023594, label = "Static FIT ofi_bin1_neg coef")
expect_near(get_val(s, "FIT", "ofi_bin1_neg", "se"),    1.3983009, label = "Static FIT ofi_bin1_neg se")
expect_near(get_val(s, "FIT", "ofi_bin3_neg", "coef"),  3.0796500, label = "Static FIT ofi_bin3_neg coef")
expect_near(get_val(s, "FIT", "ofi_bin3_neg", "se"),    0.5416002, label = "Static FIT ofi_bin3_neg se")
expect_near(get_val(s, "FIT", "ofi_bin3_neg - ofi_bin1_neg", "coef"), -0.5227094, label = "Static FIT ofi_bin3_neg - ofi_bin1_neg coef")
expect_near(get_val(s, "FIT", "ofi_bin3_neg - ofi_bin1_neg", "se"),   1.1342760,  label = "Static FIT ofi_bin3_neg - ofi_bin1_neg se")

# ---------------------------------------------------------------------------
# 2. Dynamic v4 assertions
#    File  : ../fm_dynamic_stdev_posneg_rebuilt_v4_unbalanced.RDS
#    Filter: demand_type in {OFI, FIT}, hor == 4, spec_idx == 3
# ---------------------------------------------------------------------------

dynamic <- readRDS("../fm_dynamic_stdev_posneg_rebuilt_v4_unbalanced.RDS")

d <- dynamic[spec_idx == 3L & hor == 4L]

get_dval <- function(dt, dtype, var_name, col) {
  rows <- dt[demand_type == dtype & var == var_name]
  if (nrow(rows) == 0L)
    stop(sprintf("get_dval: no row for demand_type='%s' var='%s'", dtype, var_name))
  if (nrow(rows) > 1L)
    stop(sprintf("get_dval: %d rows for demand_type='%s' var='%s'", nrow(rows), dtype, var_name))
  rows[[col]]
}

## OFI — positive side
expect_near(get_dval(d, "OFI", "ofi_bin1_pos", "coef"),  1.6328028, label = "Dynamic OFI ofi_bin1_pos coef")
expect_near(get_dval(d, "OFI", "ofi_bin1_pos", "se"),    0.0991117, label = "Dynamic OFI ofi_bin1_pos se")
expect_near(get_dval(d, "OFI", "ofi_bin3_pos", "coef"),  0.5685903, label = "Dynamic OFI ofi_bin3_pos coef")
expect_near(get_dval(d, "OFI", "ofi_bin3_pos", "se"),    0.1078185, label = "Dynamic OFI ofi_bin3_pos se")
expect_near(get_dval(d, "OFI", "ofi_bin3_pos - ofi_bin1_pos", "coef"), -1.0642125, label = "Dynamic OFI ofi_bin3_pos - ofi_bin1_pos coef")
expect_near(get_dval(d, "OFI", "ofi_bin3_pos - ofi_bin1_pos", "se"),   0.1130603,  label = "Dynamic OFI ofi_bin3_pos - ofi_bin1_pos se")

## FIT — positive side
expect_near(get_dval(d, "FIT", "ofi_bin1_pos", "coef"),  3.9742536, label = "Dynamic FIT ofi_bin1_pos coef")
expect_near(get_dval(d, "FIT", "ofi_bin1_pos", "se"),    0.3401716, label = "Dynamic FIT ofi_bin1_pos se")
expect_near(get_dval(d, "FIT", "ofi_bin3_pos", "coef"),  2.2932218, label = "Dynamic FIT ofi_bin3_pos coef")
expect_near(get_dval(d, "FIT", "ofi_bin3_pos", "se"),    0.4178467, label = "Dynamic FIT ofi_bin3_pos se")
expect_near(get_dval(d, "FIT", "ofi_bin3_pos - ofi_bin1_pos", "coef"), -1.6810318, label = "Dynamic FIT ofi_bin3_pos - ofi_bin1_pos coef")
expect_near(get_dval(d, "FIT", "ofi_bin3_pos - ofi_bin1_pos", "se"),   0.3885760,  label = "Dynamic FIT ofi_bin3_pos - ofi_bin1_pos se")

## OFI — negative side
expect_near(get_dval(d, "OFI", "ofi_bin1_neg", "coef"),  1.9503377, label = "Dynamic OFI ofi_bin1_neg coef")
expect_near(get_dval(d, "OFI", "ofi_bin1_neg", "se"),    0.0889477, label = "Dynamic OFI ofi_bin1_neg se")
expect_near(get_dval(d, "OFI", "ofi_bin3_neg", "coef"),  1.3463481, label = "Dynamic OFI ofi_bin3_neg coef")
expect_near(get_dval(d, "OFI", "ofi_bin3_neg", "se"),    0.0846016, label = "Dynamic OFI ofi_bin3_neg se")
expect_near(get_dval(d, "OFI", "ofi_bin3_neg - ofi_bin1_neg", "coef"), -0.6039895, label = "Dynamic OFI ofi_bin3_neg - ofi_bin1_neg coef")
expect_near(get_dval(d, "OFI", "ofi_bin3_neg - ofi_bin1_neg", "se"),   0.0828366,  label = "Dynamic OFI ofi_bin3_neg - ofi_bin1_neg se")

## FIT — negative side
expect_near(get_dval(d, "FIT", "ofi_bin1_neg", "coef"),  3.9687287, label = "Dynamic FIT ofi_bin1_neg coef")
expect_near(get_dval(d, "FIT", "ofi_bin1_neg", "se"),    0.4178630, label = "Dynamic FIT ofi_bin1_neg se")
expect_near(get_dval(d, "FIT", "ofi_bin3_neg", "coef"),  1.5402823, label = "Dynamic FIT ofi_bin3_neg coef")
expect_near(get_dval(d, "FIT", "ofi_bin3_neg", "se"),    0.4106457, label = "Dynamic FIT ofi_bin3_neg se")
expect_near(get_dval(d, "FIT", "ofi_bin3_neg - ofi_bin1_neg", "coef"), -2.4284464, label = "Dynamic FIT ofi_bin3_neg - ofi_bin1_neg coef")
expect_near(get_dval(d, "FIT", "ofi_bin3_neg - ofi_bin1_neg", "se"),   0.5311371,  label = "Dynamic FIT ofi_bin3_neg - ofi_bin1_neg se")

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

cat(sprintf("\nAll %d drift tests passed.\n", n_pass))
