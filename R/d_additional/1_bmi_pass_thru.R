# --- figure out whether BMI's pass-through rate depend has variability
# code modified from 26_bmi_pass_through/2_bmi_regressions.R
library(this.path)
setwd(this.path::this.dir())
source("../utilities/runmefirst.R")

# merge bmi with io changes
data <- readRDS("../tmp/raw_data/reg_inputs/reg_table_static.RDS")[type == "BMI"][, type := NULL]
tmp <- readRDS("../../../data/institutional/s34_io_changes.RDS")
data <- merge(data, tmp, by = c("yyyymm", "permno"))
rm(tmp)

# # winsorize weird stuff
# data[, dio := Winsorize(dio, quantile(dio, probs = c(0.01, 0.99)))]

# get control variable names
cdata <- readRDS("../tmp/raw_data/controls/controls_classification.RDS")
controls_char <- cdata[control_type == "return-predictor", var]
controls_liq <- cdata[control_type == "liquidity", var]
controls_list <- c(controls_char, controls_liq)
rm(cdata)

# get the bmi controls
tmp <- readRDS("..//tmp/raw_data/controls/controls_for_BMI.RDS")
controls_bmi <- setdiff(names(tmp), c("yyyymm", "permno"))
rm(tmp)

# regression specifications
formula_data <- data.table(idx = 0, formula = "dio ~ ofi_bin1 + ofi_bin2 + ofi_bin3 | yyyymm")
formula_data <- rbind(formula_data, data.table(idx = 1, formula = paste0("dio ~ ofi_bin1 + ofi_bin2 + ofi_bin3 + ", paste(controls_bmi, collapse = " + "), " | yyyymm")))
formula_data <- rbind(formula_data, data.table(idx = 2, formula = paste0("dio ~ ofi_bin1 + ofi_bin2 + ofi_bin3 + ", paste(c(controls_bmi, controls_char), collapse = " + "), " | yyyymm")))
formula_data <- rbind(formula_data, data.table(idx = 3, formula = paste0("dio ~ ofi_bin1 + ofi_bin2 + ofi_bin3 + ", paste(c(controls_bmi, controls_char, controls_liq), collapse = " + "), " | yyyymm")))


# helper function
p.get_one_regression <- function(this_data) {
  # this_data <- data_list[[1]]
  reg_out <- data.table()
  for (this_idx in formula_data[, idx]) {
    # this_idx <- formula_data[, idx]
    this_formula <- formula_data[idx == this_idx, formula]
    this_reg <- feols(as.formula(this_formula), data = this_data)
    reg_out <- rbind(reg_out, data.table(
      idx = this_idx, var = names(this_reg$coefficients)[1:3], coef = this_reg$coefficients[1:3],
      obs = nrow(this_data), r2 = r2(this_reg)["r2"]
    ))
  }
  reg_out[, yyyymm := this_data[1, yyyymm]]
  return(reg_out)
}

data_list <- split(data, by = "yyyymm")
out <- rbindlist(lapply(data_list, p.get_one_regression))

# Let's summarize including the pairwise differences
out_list <- split(out, by = "idx")

p.get_one_summary <- function(this_out) {
  # this_out <- out_list[[1]]


  X <- dcast(this_out, yyyymm ~ var, value.var = "coef")
  mu <- matrix(colMeans(X[, -1]))
  C <- cor(X[, -1]) / sqrt(nrow(X))

  out <- data.table(
    var = paste0("ofi_bin", 1:3),
    coef = colMeans(X[, -1]),
    se = sqrt(diag(C))
  )

  b <- matrix(c(1, -1, 0))
  out <- rbind(out, data.table(var = "ofi_bin1 - ofi_bin2", coef = (t(b) %*% mu)[1], se = sqrt((t(b) %*% C %*% b)[1])))
  b <- matrix(c(0, 1, -1))
  out <- rbind(out, data.table(var = "ofi_bin2 - ofi_bin3", coef = (t(b) %*% mu)[1], se = sqrt((t(b) %*% C %*% b)[1])))
  b <- matrix(c(1, 0, -1))
  out <- rbind(out, data.table(var = "ofi_bin1 - ofi_bin3", coef = (t(b) %*% mu)[1], se = sqrt((t(b) %*% C %*% b)[1])))

  out[, obs := unique(this_out[, list(yyyymm, obs)])[, sum(obs)]]
  out[, r2 := mean(this_out[, r2])]
  out[, var_idx := .I]
  out[, idx := this_out[1, idx]]
}

out <- rbindlist(lapply(out_list, p.get_one_summary))
out[idx == 0, spec_lab := "none"]
out[idx == 1, spec_lab := "bmi"]
out[idx == 2, spec_lab := "bmi + char"]
out[idx == 3, spec_lab := "bmi + char + liq"]
out[, spec_lab := paste0(idx, "_", spec_lab)]
out[, tstat := coef / se]
out[, var_lab := paste0(var_idx, "_", var)]

to_file <- "../tmp/additional/bmi_pass_thru.RDS"
dir.create(dirname(to_file), recursive = T, showWarnings = F)
saveRDS(out, to_file)

# options(width = 200)
# dcast(out, var_lab ~ spec_lab, value.var = "coef")
# dcast(out, var_lab ~ spec_lab, value.var = "tstat")
