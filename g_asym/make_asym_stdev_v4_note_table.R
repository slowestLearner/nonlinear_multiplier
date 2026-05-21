# Generate the LaTeX tabular used by asym_stdev_v4_note.tex.
#
# The table is generated from the v4 pooled six-variable output with rebuilt
# positive- and negative-side bins:
#   fm_stdev_posneg_rebuilt_v4.RDS

library(data.table)
library(scales)

out <- readRDS("fm_stdev_posneg_rebuilt_v4.RDS")

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
  paste0("  ", label, " & ", paste(cells, collapse = " & "), " \\\\")
}


pos_main_rows <- list(
  c("$M_{\\{0 < d_{n,t} < \\sigma^+\\}}$", "ofi_bin1_pos"),
  c("$M_{\\{d_{n,t} \\in [\\sigma^+, 2\\sigma^+)\\}}$", "ofi_bin2_pos"),
  c("$M_{\\{d_{n,t} \\ge 2\\sigma^+\\}}$", "ofi_bin3_pos")
)

pos_diff_rows <- list(
  c("$M_{\\{d_{n,t} \\in [\\sigma^+, 2\\sigma^+)\\}} - M_{\\{0 < d_{n,t} < \\sigma^+\\}}$", "ofi_bin2_pos - ofi_bin1_pos"),
  c("$M_{\\{d_{n,t} \\ge 2\\sigma^+\\}} - M_{\\{d_{n,t} \\in [\\sigma^+, 2\\sigma^+)\\}}$", "ofi_bin3_pos - ofi_bin2_pos"),
  c("$M_{\\{d_{n,t} \\ge 2\\sigma^+\\}} - M_{\\{0 < d_{n,t} < \\sigma^+\\}}$", "ofi_bin3_pos - ofi_bin1_pos")
)

neg_main_rows <- list(
  c("$M_{\\{-\\sigma^- < d_{n,t} < 0\\}}$", "ofi_bin1_neg"),
  c("$M_{\\{d_{n,t} \\in (-2\\sigma^-, -\\sigma^-]\\}}$", "ofi_bin2_neg"),
  c("$M_{\\{d_{n,t} \\le -2\\sigma^-\\}}$", "ofi_bin3_neg")
)

neg_diff_rows <- list(
  c("$M_{\\{d_{n,t} \\in (-2\\sigma^-, -\\sigma^-]\\}} - M_{\\{-\\sigma^- < d_{n,t} < 0\\}}$", "ofi_bin2_neg - ofi_bin1_neg"),
  c("$M_{\\{d_{n,t} \\le -2\\sigma^-\\}} - M_{\\{d_{n,t} \\in (-2\\sigma^-, -\\sigma^-]\\}}$", "ofi_bin3_neg - ofi_bin2_neg"),
  c("$M_{\\{d_{n,t} \\le -2\\sigma^-\\}} - M_{\\{-\\sigma^- < d_{n,t} < 0\\}}$", "ofi_bin3_neg - ofi_bin1_neg")
)

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

all_main_rows <- c(pos_main_rows, neg_main_rows)
all_diff_rows <- c(pos_diff_rows, neg_diff_rows)

main_lines <- c()
for (i in seq_along(all_main_rows)) {
  rw <- all_main_rows[[i]]
  main_lines <- c(main_lines, row_line(rw[1], get_result(rw[2], "coef")))
  se_label <- if (i == length(all_main_rows)) "\\vspace{5pt} " else ""
  main_lines <- c(main_lines, row_line(se_label, get_result(rw[2], "se")))
}

diff_lines <- c()
for (rw in all_diff_rows) {
  diff_lines <- c(
    diff_lines,
    row_line(rw[1], get_result(rw[2], "coef")),
    row_line("", get_result(rw[2], "se"))
  )
}

lines <- c(
  "\\begin{tabular}{lccccccccc}",
  "  \\hline \\multicolumn{10}{c}{Panel A: Regression coefficients} \\\\",
  " \\hline",
  "                                    & \\multicolumn{9}{c}{Dependent variable: stock return $r_{n,t}$} \\\\",
  "",
  "                                    \\cmidrule(l){2-10} & \\multicolumn{3}{c}{BMI} & \\multicolumn{3}{c}{FIT} & \\multicolumn{3}{c}{OFI} \\\\",
  " \\cmidrule(l){2-4} \\cmidrule(l){5-7} \\cmidrule(l){8-10} ",
  "  & (1) & (2) & (3) & (4) & (5) & (6) & (7) & (8) & (9) \\\\ ",
  main_lines,
  row_line("Predictor controls", c("N", "Y", "Y", "N", "Y", "Y", "N", "Y", "Y")),
  row_line("\\vspace{5pt}Liquidity controls", c("N", "N", "Y", "N", "N", "Y", "N", "N", "Y")),
  row_line("Obs", scalar_cells("obs", comma)),
  row_line("$R^2$", r2_cells("r2")),
  row_line("Marginal $R^2(d_{n,t})$", marginal_r2_cells),
  "   \\hline \\multicolumn{10}{c}{Panel B: Coefficient differences} \\\\",
  paste0(" \\hline", diff_lines[1]),
  diff_lines[-1],
  "   \\hline ",
  "\\end{tabular}"
)

writeLines(lines, "reg_asym_stdev_posneg_rebuilt_v4.tex")
