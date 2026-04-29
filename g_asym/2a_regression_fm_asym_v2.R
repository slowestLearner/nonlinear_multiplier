# Contemporaneous regression of returns on demand — Asymmetric specification v2
#
# PURPOSE
#   Asymmetric Fama-MacBeth regression: runs separately on positive- and
#   negative-shock subsamples and writes a combined LaTeX table.
#
# CHANGES FROM 2a_regression_fm_asym_AC.R
#   (search "[CHANGE N]" in this file to jump to each change)
#
#   [CHANGE 1] Input data: reads the POOLED reg_table_static.RDS instead of
#              the pre-split reg_table_static_{pos|neg}_ofi.RDS files produced by
#              1_construct_regression_table_static_asym_AC.R.
#              Consequence: the bin thresholds (σ, 2σ) are now based on the full
#              OFI distribution — the same σ used in 2a_regression_fm.R — making
#              the bin assignments directly comparable to the pooled results.
#              The pooled file already contains ofi_bin1/2/3 and ofi_absofi;
#              no recomputation is needed.
#
#   [CHANGE 2] Opposite-sign observations are DROPPED, not zeroed out.
#              Original approach (in 1_construct_regression_table_static_asym_AC.R):
#                data_pos[, ofi := fifelse(ofi > 0, ofi, 0)]
#              This kept ALL rows but set negative OFI to 0, so those observations
#              appeared as "zero demand" in the regression and inflated the
#              reference group. Here rows with the wrong sign are removed entirely.
#
#   [CHANGE 3] Both directions (pos and neg) are run in the same script execution,
#              enabling the table to be written at the end without a second run.
#              Original required changing pos_neg_suffix manually and re-running.
#
#   [CHANGE 4] Output RDS files use a _v2 suffix to avoid overwriting the originals.
#
#   [CHANGE 5] A LaTeX table is generated at the end of this script.
#              Original required running a separate Python script (table_asym_stdev.py).
#              The table format mirrors the reference table: two panels (Positive /
#              Negative shocks), each containing raw multipliers (with control rows
#              and obs / R2) followed by coefficient differences.

library(this.path)
setwd(this.path::this.dir())
source("../R/utilities/runmefirst.R")
source("../R/utilities/regressions.R")
library(scales)   # for comma() in obs formatting


# ===========================================================================
# SECTION 1: Control variable names
# (loaded once here; same for both pos and neg runs)
# ===========================================================================

cdata        <- readRDS("../R/tmp/raw_data/controls/controls_classification.RDS")
controls_char <- cdata[control_type == "return-predictor", var]
controls_liq  <- cdata[control_type == "liquidity", var]
controls_list <- c(controls_char, controls_liq)
rm(cdata)

tmp           <- readRDS("../R/tmp/raw_data/controls/controls_for_BMI.RDS")
controls_bmi  <- setdiff(names(tmp), c("yyyymm", "permno"))
rm(tmp)


# ===========================================================================
# SECTION 2: Fama-MacBeth regression function
# (identical to 2a_regression_fm_asym_AC.R)
# ===========================================================================

# control formulas for the three base specifications
control_formulas <- c(
  "1",
  paste0(controls_char, collapse = "+"),
  paste0(c(controls_char, controls_liq), collapse = "+")
)

# Run one demand type through all specs. reg_spec: "nonlinear" or "stdev"
p.process_one_type <- function(data, reg_spec = "nonlinear") {
  this_type <- data[1, type]
  out_all   <- data.table()

  # --- three base specifications (no interactions) ---
  for (spec_idx in 1:length(control_formulas)) {
    ff_no_ofi <- paste0("ret ~ ", control_formulas[spec_idx])
    if (reg_spec == "nonlinear") {
      ff <- paste0("ret ~ ", control_formulas[spec_idx], " + ofi + ofi_absofi")
    } else {  # "stdev"
      ff <- paste0("ret ~ ", control_formulas[spec_idx], " + ofi_bin1 + ofi_bin2 + ofi_bin3")
    }
    if (this_type == "BMI") {
      ff        <- paste0(ff,        "+", paste0(controls_bmi, collapse = "+"))
      ff_no_ofi <- paste0(ff_no_ofi, "+", paste0(controls_bmi, collapse = "+"))
    }
    out        <- p.fama_macbeth(data, ff,        compare_coefs = (reg_spec == "stdev"))
    out_no_ofi <- p.fama_macbeth(data, ff_no_ofi, compare_coefs = FALSE)
    out[, r2_no_ofi := out_no_ofi[1, r2]]
    out[, spec_idx  := spec_idx]
    out_all <- rbind(out_all, out)
  }

  # --- progressively add characteristic interactions ---
  for (this_v in c(controls_char, controls_liq)) {
    setnames(data, this_v, "xx")
    data[, yy := xx * ofi]
    setnames(data, c("xx", "yy"), c(this_v, paste0("ofi_", this_v)))
    ff        <- paste0(ff,        " + ofi_", this_v)
    ff_no_ofi <- paste0(ff_no_ofi, " + ofi_", this_v)
    spec_idx  <- spec_idx + 1

    out        <- p.fama_macbeth(data, ff,        compare_coefs = (reg_spec == "stdev"))
    out_no_ofi <- p.fama_macbeth(data, ff_no_ofi, compare_coefs = FALSE)
    out[, r2_no_ofi := out_no_ofi[1, r2]]
    out[, spec_idx  := spec_idx]
    out_all <- rbind(out_all, out)
  }

  # label each spec with what was added
  spec_labels <- data.table(
    spec_idx  = 1:(length(controls_list) + 3),
    var_added = c("none_init", "controls_char", "controls_char+controls_liq", controls_list),
    var_type  = c(rep("", 3),
                  rep("return-predicting chars", length(controls_char)),
                  rep("liquidity",               length(controls_liq)))
  )
  out_all <- merge(out_all, spec_labels, by = "spec_idx", all.x = TRUE)
  out_all[, type := this_type]
  return(out_all)
}


# ===========================================================================
# SECTION 3: Run regressions for positive and negative shocks
# [CHANGE 3] Both directions run in one execution
# ===========================================================================

to_dir <- "./"
dir.create(to_dir, recursive = TRUE, showWarnings = FALSE)

results_stdev <- list()  # store stdev results for table generation

for (pos_neg_suffix in c("pos", "neg")) {

  cat("\n=== Running:", pos_neg_suffix, "===\n")

  # -----------------------------------------------------------------------
  # [CHANGE 1] Read the POOLED regression table.
  # Original: readRDS(paste0("...reg_table_static_", pos_neg_suffix, "_ofi.RDS"))
  # New:      readRDS("...reg_table_static.RDS")
  # This file was built by 1_construct_regression_table_static.R using
  # sd(ofi) on the full OFI distribution, so ofi_bin1/2/3 already reflect
  # the pooled thresholds — no recomputation needed.
  # -----------------------------------------------------------------------
  data_all <- readRDS("../R/tmp/raw_data/reg_inputs/reg_table_static.RDS")

  # -----------------------------------------------------------------------
  # [CHANGE 2] Drop opposite-sign rows entirely.
  # Original: set negative OFI to 0 (kept rows, inflated "zero demand" group).
  # New:      remove rows with the wrong sign so the regression is estimated
  #           on the sign-specific subsample only.
  # -----------------------------------------------------------------------
  if (pos_neg_suffix == "pos") {
    data_all <- data_all[ofi > 0]  # [CHANGE 2] drop negative-OFI rows
  } else {
    data_all <- data_all[ofi < 0]  # [CHANGE 2] drop positive-OFI rows
  }
  # The retained rows already have correct ofi_bin1/2/3 values:
  #   pos rows: ofi_bin{k} = ofi > 0 for the relevant bin, 0 elsewhere
  #   neg rows: ofi_bin{k} = ofi < 0 for the relevant bin, 0 elsewhere
  cat("Rows after sign filter:", nrow(data_all), "\n")
  print(data_all[, .(
    n = .N,
    mean_ofi = mean(ofi),
    sd_ofi = sd(ofi),
    p01 = quantile(ofi, 0.01),
    p50 = quantile(ofi, 0.50),
    p99 = quantile(ofi, 0.99)
  ), by = type][order(type)])
  print(data_all[, .(
    bin1_obs = sum(ofi_bin1 != 0),
    bin2_obs = sum(ofi_bin2 != 0),
    bin3_obs = sum(ofi_bin3 != 0),
    mean_abs_bin1 = mean(abs(ofi[ofi_bin1 != 0])),
    mean_abs_bin2 = mean(abs(ofi[ofi_bin2 != 0])),
    mean_abs_bin3 = mean(abs(ofi[ofi_bin3 != 0]))
  ), by = type][order(type)])

  # --- nonlinear specification ---
  tic(paste("static fm nonlinear:", pos_neg_suffix))
  out_nonlinear <- rbindlist(mclapply(split(data_all, by = "type"), function(x) {
    p.process_one_type(x, reg_spec = "nonlinear")
  }, mc.cores = nc))
  toc()
  # [CHANGE 4] _v2 suffix avoids overwriting the original output files
  saveRDS(out_nonlinear, paste0(to_dir, "fm_nonlinear_", pos_neg_suffix, "_v2.RDS"))

  # --- stdev-based specification ---
  tic(paste("static fm stdev:", pos_neg_suffix))
  out_stdev <- rbindlist(mclapply(split(data_all, by = "type"), function(x) {
    p.process_one_type(x, reg_spec = "stdev")
  }, mc.cores = nc))
  toc()
  # [CHANGE 4] _v2 suffix avoids overwriting the original output files
  saveRDS(out_stdev, paste0(to_dir, "fm_stdev_", pos_neg_suffix, "_v2.RDS"))

  results_stdev[[pos_neg_suffix]] <- out_stdev
  rm(data_all, out_nonlinear, out_stdev)
}


# ===========================================================================
# SECTION 4: Generate LaTeX table
# [CHANGE 5] Table is produced here; no separate Python script needed.
#
# Table structure (mirrors the reference table format):
#   Panel A: Positive shocks
#     — raw multipliers: ofi_bin1, ofi_bin2, ofi_bin3 (coef + SE per row)
#     — control indicators, Obs, R2, Marginal R2
#     — coefficient differences: bin2-bin1, bin3-bin2, bin3-bin1
#   Panel B: Negative shocks  (same structure)
# ===========================================================================

# --- 4a: helper functions ---

cuts <- abs(qnorm(c(.01, .05, .1) / 2))  # 2.576, 1.960, 1.645

fmt_coef <- function(coef, se) {
  n_stars  <- sum(abs(coef / se) > cuts)
  stars    <- if (n_stars > 0) paste0("^{", strrep("*", n_stars), "}$") else "$"
  paste0("$", sprintf("%.2f", coef), stars)
}

fmt_se <- function(se) paste0("$(", sprintf("%.2f", se), ")$")

# Return 9 formatted cells (BMI spec1-3 | FIT spec1-3 | OFI spec1-3) for one variable
cells_coef <- function(data, var_name) {
  vapply(list(
    c("BMI",1), c("BMI",2), c("BMI",3),
    c("FIT",1), c("FIT",2), c("FIT",3),
    c("OFI",1), c("OFI",2), c("OFI",3)
  ), function(ts) {
    r <- data[type == ts[1] & spec_idx == as.integer(ts[2]) & var == var_name]
    fmt_coef(r$coef, r$se)
  }, character(1))
}

cells_se <- function(data, var_name) {
  vapply(list(
    c("BMI",1), c("BMI",2), c("BMI",3),
    c("FIT",1), c("FIT",2), c("FIT",3),
    c("OFI",1), c("OFI",2), c("OFI",3)
  ), function(ts) {
    r <- data[type == ts[1] & spec_idx == as.integer(ts[2]) & var == var_name]
    fmt_se(r$se)
  }, character(1))
}

cells_scalar <- function(data, var_name, fmt_fn) {
  vapply(list(
    c("BMI",1), c("BMI",2), c("BMI",3),
    c("FIT",1), c("FIT",2), c("FIT",3),
    c("OFI",1), c("OFI",2), c("OFI",3)
  ), function(ts) {
    r <- data[type == ts[1] & spec_idx == as.integer(ts[2]) & var == var_name]
    fmt_fn(r)
  }, character(1))
}

tex_row <- function(label, cells, vspace = FALSE) {
  v <- if (vspace) " \\vspace{5pt}" else ""
  paste0("  ", label, " & ", paste(cells, collapse = " & "), " \\\\", v, " \n")
}


# --- 4b: variable labels ---

BIN_VARS <- paste0("ofi_bin", 1:3)
BIN_LABS <- c(
  "$M_{\\{|d_{n,t}| < \\sigma\\}}$",
  "$M_{\\{|d_{n,t}| \\in [\\sigma, 2\\sigma]\\}}$",
  "$M_{\\{|d_{n,t}| > 2 \\sigma\\}}$"
)

DIFF_VARS <- c("ofi_bin2 - ofi_bin1", "ofi_bin3 - ofi_bin2", "ofi_bin3 - ofi_bin1")
DIFF_LABS <- c(
  "$M_{\\{|d_{n,t}| \\in [\\sigma, 2\\sigma]\\}} - M_{\\{|d_{n,t}| < \\sigma\\}}$",
  "$M_{\\{|d_{n,t}| > 2 \\sigma\\}} - M_{\\{|d_{n,t}| \\in [\\sigma, 2\\sigma]\\}}$",
  "$M_{\\{|d_{n,t}| > 2 \\sigma\\}} - M_{\\{|d_{n,t}| < \\sigma\\}}$"
)


# --- 4c: build one panel block (raw multipliers + stats + differences) ---

build_panel <- function(data, panel_label) {
  # filter to spec_idx 1:3 and relevant types only
  d <- data[spec_idx %in% 1:3 & type %in% c("BMI", "FIT", "OFI")]

  lines <- character(0)

  # panel header
  lines <- c(lines,
    paste0("  \\hline \\multicolumn{10}{c}{", panel_label, "} \\\\\n"),
    " \\hline\n",
    "  & \\multicolumn{3}{c}{BMI} & \\multicolumn{3}{c}{FIT} & \\multicolumn{3}{c}{OFI} \\\\ \n",
    "  \\cmidrule(l){2-4} \\cmidrule(l){5-7} \\cmidrule(l){8-10} \n",
    "  & (1) & (2) & (3) & (4) & (5) & (6) & (7) & (8) & (9) \\\\ \n"
  )

  # raw multipliers (bin1, bin2, bin3)
  for (i in seq_along(BIN_VARS)) {
    is_last <- (i == length(BIN_VARS))
    lines <- c(lines,
      tex_row(BIN_LABS[i], cells_coef(d, BIN_VARS[i])),
      tex_row("",           cells_se(d,   BIN_VARS[i]), vspace = is_last)
    )
  }

  # control indicator rows (N/Y/Y pattern repeats across all three types)
  pred_ctrls <- rep(c("N", "Y", "Y"), 3)
  liq_ctrls  <- rep(c("N", "N", "Y"), 3)
  lines <- c(lines,
    tex_row("Predictor controls", pred_ctrls),
    tex_row("Liquidity controls", liq_ctrls, vspace = TRUE)
  )

  # obs, R2, marginal R2
  obs_cells <- cells_scalar(d, "ofi_bin1", function(r) comma(r$obs))
  r2_cells  <- cells_scalar(d, "ofi_bin1", function(r) paste0("$", sprintf("%.3f", r$r2), "$"))
  mr2_cells <- cells_scalar(d, "ofi_bin1", function(r) paste0("$", sprintf("%.3f", r$r2 - r$r2_no_ofi), "$"))

  lines <- c(lines,
    tex_row("Obs",               obs_cells),
    tex_row("$R^2$",             r2_cells),
    tex_row("Marginal $R^2(d_{n,t})$", mr2_cells)
  )

  # section divider before differences
  lines <- c(lines, " \\hline\n")

  # coefficient differences
  for (i in seq_along(DIFF_VARS)) {
    is_last <- (i == length(DIFF_VARS))
    lines <- c(lines,
      tex_row(DIFF_LABS[i], cells_coef(d, DIFF_VARS[i])),
      tex_row("",            cells_se(d,   DIFF_VARS[i]), vspace = is_last)
    )
  }

  paste(lines, collapse = "")
}


# --- 4d: assemble and write ---

pos_panel <- build_panel(results_stdev[["pos"]], "Panel A: Positive shocks")
neg_panel <- build_panel(results_stdev[["neg"]], "Panel B: Negative shocks")

tabular <- paste0(
  "\\begin{tabular}{lccccccccc}\n",
  pos_panel,
  neg_panel,
  "   \\hline \n",
  "\\end{tabular}\n"
)

out_dir  <- "./"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_path <- paste0(out_dir, "reg_asym_stdev_v2.tex")   # [CHANGE 4] new filename
writeLines(tabular, out_path)
cat("Table written to:", out_path, "\n")


# ===========================================================================
# SECTION 5: Sanity check (compare v2 vs v1 multipliers)
# Uncomment to run after both v1 and v2 results are available.
# Expect: v2 multipliers smaller than v1 (pooled σ > zeroed-out σ → larger
# bin ranges → smaller per-unit multipliers within each bin).
# ===========================================================================
# v2_pos <- readRDS(paste0(to_dir, "fm_stdev_pos_v2.RDS"))
# v1_pos <- readRDS(paste0(to_dir, "fm_stdev_pos.RDS"))
#
# compare <- merge(
#   v2_pos[spec_idx %in% 1:3 & grepl("ofi_bin[123]$", var) & type %in% c("BMI","FIT","OFI"),
#          .(type, var, spec_idx, coef_v2 = coef, se_v2 = se)],
#   v1_pos[spec_idx %in% 1:3 & grepl("ofi_bin[123]$", var) & type %in% c("BMI","FIT","OFI"),
#          .(type, var, spec_idx, coef_v1 = coef, se_v1 = se)],
#   by = c("type", "var", "spec_idx")
# )
# compare[, ratio := coef_v2 / coef_v1]
# print(compare[order(type, spec_idx, var)])
