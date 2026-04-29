# Generate the dynamic bin-summary LaTeX table used by asym_stdev_v4_note.tex.
#
# The dynamic v4 regression bins observations by lagged cumulative demand
# shocks, while the regressor inside each bin is contemporaneous demand. This
# table reports both variables for the 4-lag FIT/OFI specification used in the
# dynamic note table.

library(data.table)
library(scales)

data <- as.data.table(readRDS("../R/tmp/raw_data/reg_inputs/reg_table_dynamic.RDS"))

vars_id <- c("yyyymm", "permno", "type")
vars_reg <- c("ret", "ofi", "cumofi_1", "cumofi_2", "cumofi_3", "cumofi_4")

data_reg <- data[, c(vars_id, vars_reg), with = FALSE]
data_reg <- melt(
  data_reg,
  id.vars = c("yyyymm", "permno", "type", "ret", "ofi"),
  variable.name = "hor",
  value.name = "cumofi_lag"
)
data_reg <- data_reg[0 == rowSums(is.na(data_reg))]
data_reg[, hor := as.integer(gsub("cumofi_", "", hor))]
data_reg <- data_reg[type %in% c("FIT", "OFI") & hor == 4L]

data_reg[, cumofi_sd_pos := sd(cumofi_lag[cumofi_lag > 0]), by = .(yyyymm, type, hor)]
data_reg[, cumofi_sd_neg := sd(cumofi_lag[cumofi_lag < 0]), by = .(yyyymm, type, hor)]

summary_data <- data_reg[, .(
  shock = fcase(
    cumofi_lag > 0, "Positive",
    cumofi_lag < 0, "Negative",
    default = NA_character_
  ),
  bin = fcase(
    cumofi_lag > 0 & cumofi_lag < cumofi_sd_pos, 1L,
    cumofi_lag >= cumofi_sd_pos & cumofi_lag < 2 * cumofi_sd_pos, 2L,
    cumofi_lag >= 2 * cumofi_sd_pos, 3L,
    cumofi_lag < 0 & cumofi_lag > -cumofi_sd_neg, 1L,
    cumofi_lag <= -cumofi_sd_neg & cumofi_lag > -2 * cumofi_sd_neg, 2L,
    cumofi_lag <= -2 * cumofi_sd_neg, 3L,
    default = NA_integer_
  )
)]
summary_data[, `:=`(
  type = data_reg$type,
  cumofi_lag = data_reg$cumofi_lag,
  ofi = data_reg$ofi
)]

summary_data <- summary_data[!is.na(shock) & !is.na(bin), .(
  obs = .N,
  mean_cumofi_lag = mean(cumofi_lag),
  sd_cumofi_lag = sd(cumofi_lag),
  mean_ofi = mean(ofi),
  sd_ofi = sd(ofi)
), by = .(shock, type, bin)]

summary_data[, shock_order := fifelse(shock == "Negative", 1L, 2L)]
setorder(summary_data, shock_order, type, bin)

fmt_num <- function(x) sprintf("%.6f", x)

row_line <- function(r) {
  paste0(
    r$shock, " & ", r$type, " & ", r$bin, " & ", comma(r$obs),
    " & ", fmt_num(r$mean_cumofi_lag),
    " & ", fmt_num(r$sd_cumofi_lag),
    " & ", fmt_num(r$mean_ofi),
    " & ", fmt_num(r$sd_ofi),
    " \\\\"
  )
}

lines <- c(
  "\\begin{tabular}{lllrrrrr}",
  "\\toprule",
  "Shock & Type & Bin & Obs & Mean $c^{(4)}_{n,t}$ & SD $c^{(4)}_{n,t}$ & Mean $d_{n,t}$ & SD $d_{n,t}$ \\\\",
  "\\midrule"
)

for (i in seq_len(nrow(summary_data))) {
  if (i > 1 && summary_data$shock[i] != summary_data$shock[i - 1]) {
    lines <- c(lines, "\\addlinespace")
  }
  lines <- c(lines, row_line(summary_data[i]))
}

lines <- c(lines, "\\bottomrule", "\\end{tabular}")

writeLines(lines, "summary_dynamic_asym_stdev_posneg_rebuilt_v4_4lag.tex")
