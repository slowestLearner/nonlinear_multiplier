# Generate LaTeX tables for the v4 tests of
#   Delta_{j-k}(+) = Delta_{j-k}(-).
#
# Tables use the most-controlled specification, spec_idx = 3. The dynamic table
# reports FIT and OFI only; OFI_pre_whitened is intentionally excluded from the
# note table.

library(data.table)

fmt_num <- function(x) sprintf("%.2f", x)

delta_label <- function(x) {
  fifelse(x == "2-1", "$\\Delta_{2-1}$",
          fifelse(x == "3-2", "$\\Delta_{3-2}$", "$\\Delta_{3-1}$"))
}

row_static <- function(r) {
  paste0(
    r$type, " & ", delta_label(r$delta),
    " & ", fmt_num(r$mean_delta_pos),
    " & ", fmt_num(r$mean_delta_neg),
    " & ", fmt_num(r$coef),
    " & ", fmt_num(r$se),
    " & ", fmt_num(r$t),
    # " & ", r$n_months,
    " \\\\"
  )
}

row_dynamic <- function(r) {
  paste0(
    r$demand_type, " & ", r$hor, " & ", delta_label(r$delta),
    " & ", fmt_num(r$mean_delta_pos),
    " & ", fmt_num(r$mean_delta_neg),
    " & ", fmt_num(r$coef),
    " & ", fmt_num(r$se),
    " & ", fmt_num(r$t),
    # " & ", r$n_months,
    " \\\\"
  )
}

row_dynamic_nolag <- function(r) {
  paste0(
    r$demand_type, " & ", delta_label(r$delta),
    " & ", fmt_num(r$mean_delta_pos),
    " & ", fmt_num(r$mean_delta_neg),
    " & ", fmt_num(r$coef),
    " & ", fmt_num(r$se),
    " & ", fmt_num(r$t),
    # " & ", r$n_months,
    " \\\\"
  )
}

static <- fread("v4_delta_posneg_tests_static.csv")
static <- static[spec_idx == 3]
static[, type := factor(type, levels = c("BMI", "FIT", "OFI"))]
static[, delta := factor(delta, levels = c("2-1", "3-2", "3-1"))]
setorder(static, type, delta)

dynamic <- fread("v4_delta_posneg_tests_dynamic.csv")
dynamic <- dynamic[spec_idx == 3 & demand_type %in% c("FIT", "OFI")]
dynamic[, demand_type := factor(demand_type, levels = c("FIT", "OFI"))]
dynamic[, delta := factor(delta, levels = c("2-1", "3-2", "3-1"))]
setorder(dynamic, demand_type, hor, delta)

static_lines <- c(
  "\\begin{tabular}{llrrrrr}",
  "\\toprule",
  "Type & Contrast & $\\Delta(+)$ & $\\Delta(-)$ & $\\Delta(+) - \\Delta(-)$ & SE & $t$ \\\\",
  "\\midrule"
)

for (i in seq_len(nrow(static))) {
  if (i > 1 && static$type[i] != static$type[i - 1]) {
    static_lines <- c(static_lines, "\\addlinespace")
  }
  static_lines <- c(static_lines, row_static(static[i]))
}

static_lines <- c(static_lines, "\\bottomrule", "\\end{tabular}")
writeLines(static_lines, "test_v4_delta_posneg_static_spec3.tex")

dynamic_lines <- c(
  "\\begin{tabular}{llrrrrrr}",
  "\\toprule",
  "Type & Lag & Contrast & $\\Delta(+)$ & $\\Delta(-)$ & $\\Delta(+) - \\Delta(-)$ & SE & $t$ \\\\",
  "\\midrule"
)

for (i in seq_len(nrow(dynamic))) {
  if (i > 1 && dynamic$demand_type[i] != dynamic$demand_type[i - 1]) {
    dynamic_lines <- c(dynamic_lines, "\\addlinespace")
  }
  dynamic_lines <- c(dynamic_lines, row_dynamic(dynamic[i]))
}

dynamic_lines <- c(dynamic_lines, "\\bottomrule", "\\end{tabular}")
writeLines(dynamic_lines, "test_v4_delta_posneg_dynamic_spec3.tex")

dynamic_lag4 <- dynamic[hor == 4]
dynamic_lag4_lines <- c(
  "\\begin{tabular}{lrrrrrr}",
  "\\toprule",
  "Type & Contrast & $\\Delta(+)$ & $\\Delta(-)$ & $\\Delta(+) - \\Delta(-)$ & SE & $t$ \\\\",
  "\\midrule"
)

for (i in seq_len(nrow(dynamic_lag4))) {
  if (i > 1 && dynamic_lag4$demand_type[i] != dynamic_lag4$demand_type[i - 1]) {
    dynamic_lag4_lines <- c(dynamic_lag4_lines, "\\addlinespace")
  }
  dynamic_lag4_lines <- c(dynamic_lag4_lines, row_dynamic_nolag(dynamic_lag4[i]))
}

dynamic_lag4_lines <- c(dynamic_lag4_lines, "\\bottomrule", "\\end{tabular}")
writeLines(dynamic_lag4_lines, "test_v4_delta_posneg_dynamic_spec3_lag4.tex")
