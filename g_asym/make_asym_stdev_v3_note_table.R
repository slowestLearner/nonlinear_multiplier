# Generate the LaTeX tabular used by asym_stdev_v3_note.tex.
#
# The table is generated from the v3 pooled six-variable output:
#   fm_stdev_posneg_v3.RDS

library(data.table)
library(scales)

out <- readRDS("fm_stdev_posneg_v3.RDS")

cuts <- abs(qnorm(c(.01, .05, .1) / 2))

fmt_coef <- function(coef, se) {
  n_stars <- sum(abs(coef / se) > cuts)
  stars <- if (n_stars > 0) paste0("^{", strrep("*", n_stars), "}$") else "$"
  paste0("$", sprintf("%.2f", coef), stars)
}

fmt_se <- function(se) paste0("$(", sprintf("%.2f", se), ")$")

cols <- list(
  c("BMI", 1), c("BMI", 2), c("BMI", 3),
  c("FIT", 1), c("FIT", 2), c("FIT", 3),
  c("OFI", 1), c("OFI", 2), c("OFI", 3)
)

get_result <- function(var_name, what) {
  vapply(cols, function(ts) {
    r <- out[type == ts[1] & spec_idx == as.integer(ts[2]) & var == var_name]
    if (nrow(r) != 1) stop("Missing/non-unique result for ", var_name)
    if (what == "coef") fmt_coef(r$coef, r$se) else fmt_se(r$se)
  }, character(1))
}

row_line <- function(label, cells) {
  paste0("  ", label, " & ", paste(cells, collapse = " & "), " \\\\\n")
}

panel_lines <- function(panel_name, rows) {
  lines <- c(
    paste0("  \\hline \\multicolumn{10}{c}{", panel_name, "} \\\\\n"),
    "  \\hline\n",
    "  & \\multicolumn{3}{c}{BMI} & \\multicolumn{3}{c}{FIT} & \\multicolumn{3}{c}{OFI} \\\\\n",
    "  \\cmidrule(lr){2-4} \\cmidrule(lr){5-7} \\cmidrule(lr){8-10}\n",
    "  & (1) & (2) & (3) & (1) & (2) & (3) & (1) & (2) & (3) \\\\\n"
  )
  for (rw in rows) {
    lines <- c(
      lines,
      row_line(rw[1], get_result(rw[2], "coef")),
      row_line("", get_result(rw[2], "se"))
    )
  }
  lines
}

pos_rows <- list(
  c("Bin 1", "ofi_bin1_pos"),
  c("Bin 2", "ofi_bin2_pos"),
  c("Bin 3", "ofi_bin3_pos"),
  c("Bin 2 - Bin 1", "ofi_bin2_pos - ofi_bin1_pos"),
  c("Bin 3 - Bin 2", "ofi_bin3_pos - ofi_bin2_pos"),
  c("Bin 3 - Bin 1", "ofi_bin3_pos - ofi_bin1_pos")
)

neg_rows <- list(
  c("Bin 1", "ofi_bin1_neg"),
  c("Bin 2", "ofi_bin2_neg"),
  c("Bin 3", "ofi_bin3_neg"),
  c("Bin 2 - Bin 1", "ofi_bin2_neg - ofi_bin1_neg"),
  c("Bin 3 - Bin 2", "ofi_bin3_neg - ofi_bin2_neg"),
  c("Bin 3 - Bin 1", "ofi_bin3_neg - ofi_bin1_neg")
)

get_scalar <- function(type, spec_idx, field) {
  out[type == type & spec_idx == spec_idx][1][[field]]
}

scalar_cells <- function(field, transform = identity) {
  vapply(cols, function(ts) {
    r <- out[type == ts[1] & spec_idx == as.integer(ts[2])][1]
    transform(r[[field]])
  }, character(1))
}

r2_cells <- function(field) {
  vapply(cols, function(ts) {
    r <- out[type == ts[1] & spec_idx == as.integer(ts[2])][1]
    sprintf("%.3f", r[[field]])
  }, character(1))
}

marginal_r2_cells <- vapply(cols, function(ts) {
  r <- out[type == ts[1] & spec_idx == as.integer(ts[2])][1]
  sprintf("%.3f", r$r2 - r$r2_no_ofi)
}, character(1))

lines <- c(
  "\\begin{tabular}{lccccccccc}\n",
  panel_lines("Panel A: Positive shocks", pos_rows),
  panel_lines("Panel B: Negative shocks", neg_rows),
  "  \\hline\n",
  row_line("Predictor controls", c("No", "Yes", "Yes", "No", "Yes", "Yes", "No", "Yes", "Yes")),
  row_line("Liquidity controls", c("No", "No", "Yes", "No", "No", "Yes", "No", "No", "Yes")),
  row_line("Observations", scalar_cells("obs", comma)),
  row_line("$R^2$", r2_cells("r2")),
  row_line("Marginal $R^2$", marginal_r2_cells),
  "  \\hline\n",
  "\\end{tabular}\n"
)

writeLines(lines, "reg_asym_stdev_posneg_v3.tex")
