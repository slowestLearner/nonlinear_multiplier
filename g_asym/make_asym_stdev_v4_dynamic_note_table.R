# Generate the dynamic LaTeX tabular used by asym_stdev_v4_note.tex.
#
# The table is generated from the v4 dynamic pooled six-variable output with
# rebuilt positive- and negative-side bins. The note reports the 4-lag horizon;
# balanced and unbalanced v4 dynamic outputs are identical because the balanced
# filter drops zero type-horizon-month cells.

library(data.table)
library(scales)

out <- as.data.table(readRDS("fm_dynamic_stdev_posneg_rebuilt_v4_unbalanced.RDS"))

cuts <- abs(qnorm(c(.01, .05, .1) / 2))

fmt_coef <- function(coef, se) {
  n_stars <- sum(abs(coef / se) > cuts)
  stars <- if (n_stars > 0) paste0("^{", strrep("*", n_stars), "}$") else "$"
  paste0("$", sprintf("%.2f", coef), stars)
}

fmt_se <- function(se) paste0("$(", sprintf("%.2f", se), ")$")

cols <- list(
  c("FIT", 1), c("FIT", 2), c("FIT", 3),
  c("OFI", 1), c("OFI", 2), c("OFI", 3)
)

get_result <- function(var_name, what) {
  vapply(cols, function(ts) {
    r <- out[demand_type == ts[1] & hor == 4L & spec_idx == as.integer(ts[2]) &
               var == var_name]
    if (nrow(r) != 1) stop("Missing/non-unique result for ", var_name)
    if (what == "coef") fmt_coef(r$coef, r$se) else fmt_se(r$se)
  }, character(1))
}

row_line <- function(label, cells) {
  paste0("  ", label, " & ", paste(cells, collapse = " & "), " \\\\")
}

panel_lines <- function(panel_name, main_rows, diff_rows) {
  lines <- c(
    paste0("  \\hline \\multicolumn{7}{c}{", panel_name, "} \\\\"),
    "  \\hline",
    "  & \\multicolumn{3}{c}{FIT} & \\multicolumn{3}{c}{OFI} \\\\",
    "  \\cmidrule(l){2-4} \\cmidrule(l){5-7}",
    "  & (1) & (2) & (3) & (4) & (5) & (6) \\\\"
  )
  for (rw in main_rows) {
    lines <- c(
      lines,
      row_line(rw[1], get_result(rw[2], "coef")),
      row_line("", get_result(rw[2], "se"))
    )
  }
  lines <- c(
    lines,
    "  \\vspace{5pt}",
    row_line("Predictor controls", c("N", "Y", "Y", "N", "Y", "Y")),
    row_line("Liquidity controls", c("N", "N", "Y", "N", "N", "Y")),
    "  \\vspace{5pt}",
    row_line("Obs", scalar_cells("obs", comma)),
    row_line("$R^2$", r2_cells("r2")),
    row_line("Marginal $R^2(d_{n,t})$", marginal_r2_cells),
    "  \\hline"
  )
  for (rw in diff_rows) {
    lines <- c(
      lines,
      row_line(rw[1], get_result(rw[2], "coef")),
      row_line("", get_result(rw[2], "se"))
    )
  }
  c(lines, "  \\vspace{5pt}")
}

pos_main_rows <- list(
  c("$M_{\\{0 < c^{(4)}_{n,t} < \\sigma^+\\}}$", "ofi_bin1_pos"),
  c("$M_{\\{c^{(4)}_{n,t} \\in [\\sigma^+, 2\\sigma^+)\\}}$", "ofi_bin2_pos"),
  c("$M_{\\{c^{(4)}_{n,t} \\ge 2\\sigma^+\\}}$", "ofi_bin3_pos")
)

pos_diff_rows <- list(
  c("$M_{\\{c^{(4)}_{n,t} \\in [\\sigma^+, 2\\sigma^+)\\}} - M_{\\{0 < c^{(4)}_{n,t} < \\sigma^+\\}}$", "ofi_bin2_pos - ofi_bin1_pos"),
  c("$M_{\\{c^{(4)}_{n,t} \\ge 2\\sigma^+\\}} - M_{\\{c^{(4)}_{n,t} \\in [\\sigma^+, 2\\sigma^+)\\}}$", "ofi_bin3_pos - ofi_bin2_pos"),
  c("$M_{\\{c^{(4)}_{n,t} \\ge 2\\sigma^+\\}} - M_{\\{0 < c^{(4)}_{n,t} < \\sigma^+\\}}$", "ofi_bin3_pos - ofi_bin1_pos")
)

neg_main_rows <- list(
  c("$M_{\\{-\\sigma^- < c^{(4)}_{n,t} < 0\\}}$", "ofi_bin1_neg"),
  c("$M_{\\{c^{(4)}_{n,t} \\in (-2\\sigma^-, -\\sigma^-]\\}}$", "ofi_bin2_neg"),
  c("$M_{\\{c^{(4)}_{n,t} \\le -2\\sigma^-\\}}$", "ofi_bin3_neg")
)

neg_diff_rows <- list(
  c("$M_{\\{c^{(4)}_{n,t} \\in (-2\\sigma^-, -\\sigma^-]\\}} - M_{\\{-\\sigma^- < c^{(4)}_{n,t} < 0\\}}$", "ofi_bin2_neg - ofi_bin1_neg"),
  c("$M_{\\{c^{(4)}_{n,t} \\le -2\\sigma^-\\}} - M_{\\{c^{(4)}_{n,t} \\in (-2\\sigma^-, -\\sigma^-]\\}}$", "ofi_bin3_neg - ofi_bin2_neg"),
  c("$M_{\\{c^{(4)}_{n,t} \\le -2\\sigma^-\\}} - M_{\\{-\\sigma^- < c^{(4)}_{n,t} < 0\\}}$", "ofi_bin3_neg - ofi_bin1_neg")
)

scalar_cells <- function(field, transform = identity) {
  vapply(cols, function(ts) {
    r <- out[demand_type == ts[1] & hor == 4L & spec_idx == as.integer(ts[2])][1]
    transform(r[[field]])
  }, character(1))
}

r2_cells <- function(field) {
  vapply(cols, function(ts) {
    r <- out[demand_type == ts[1] & hor == 4L & spec_idx == as.integer(ts[2])][1]
    sprintf("%.3f", r[[field]])
  }, character(1))
}

marginal_r2_cells <- vapply(cols, function(ts) {
  r <- out[demand_type == ts[1] & hor == 4L & spec_idx == as.integer(ts[2])][1]
  sprintf("%.3f", r$r2 - r$r2_no_ofi)
}, character(1))

lines <- c(
  "\\begin{tabular}{lcccccc}",
  panel_lines("Panel A: Positive lagged cumulative shocks", pos_main_rows, pos_diff_rows),
  panel_lines("Panel B: Negative lagged cumulative shocks", neg_main_rows, neg_diff_rows),
  "  \\hline",
  "\\end{tabular}"
)

writeLines(lines, "reg_dynamic_asym_stdev_posneg_rebuilt_v4_4lag.tex")
