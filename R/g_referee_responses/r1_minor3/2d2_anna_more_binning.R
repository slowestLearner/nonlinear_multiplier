# --- Can add bin-specific controls
library(this.path)
setwd(this.path::this.dir())
source("../../utilities/runmefirst.R")
source("../../utilities/regressions.R")

# load regression data
tic("preparing data")
data_all <- readRDS("../../tmp/raw_data/reg_inputs/reg_table_static.RDS")[type %in% c("BMI", "FIT", "OFI")]

# get control variable names
cdata <- readRDS("../../tmp/raw_data/controls/controls_classification.RDS")
controls_char <- cdata[control_type == "return-predictor", var]
controls_liq <- cdata[control_type == "liquidity", var]
controls_list <- c(controls_char, controls_liq)
rm(cdata)

# get names for bmi controls
tmp <- readRDS("../../tmp/raw_data/controls/controls_for_BMI.RDS")
controls_bmi <- setdiff(names(tmp), c("yyyymm", "permno"))
rm(tmp)
toc()

# --- TODO: put this into utility function later
# ff <- paste0("ret ~ ", control_formulas[spec_idx], " + ofi_bin1 + ofi_bin2 + ofi_bin3")
# compare_coefs <- TRUE
# output_cov <- TRUE
p.fama_macbeth_with_cov <- function(data, ff, compare_coefs = FALSE, output_cov = FALSE) {
  # regression for one period
  p.get_one_period <- function(this_ym) {
    ols <- lm(ff, data[yyyymm == this_ym])
    # Get the dependent variable name dynamically
    dep_var_name <- all.vars(as.formula(ff))[1]
    return(data.table(
      yyyymm = this_ym, var = names(coef(ols)), coef = ols$coef,
      r2 = var(ols$fitted.values) / var(ols$model[[dep_var_name]])
    ))
  }

  # regression by period
  out <- rbindlist(lapply(unique(data[, yyyymm]), p.get_one_period))

  # extract R2
  r2_data <- unique(out[, .(yyyymm, r2)])[, .(r2 = mean(r2))]
  out[, r2 := NULL]

  # let's do Newey-West
  yms <- sort(unique(out[, yyyymm]))

  # # compare coefficient differences
  # if (compare_coefs) {
  #   out <- rbind(out, data.table(
  #     yyyymm = yms, var = "ofi_bin2 - ofi_bin1",
  #     coef = out[var == "ofi_bin2", coef] - out[var == "ofi_bin1", coef]
  #   ))
  #   out <- rbind(out, data.table(
  #     yyyymm = yms, var = "ofi_bin3 - ofi_bin1",
  #     coef = out[var == "ofi_bin3", coef] - out[var == "ofi_bin1", coef]
  #   ))
  #   out <- rbind(out, data.table(
  #     yyyymm = yms, var = "ofi_bin3 - ofi_bin2",
  #     coef = out[var == "ofi_bin3", coef] - out[var == "ofi_bin2", coef]
  #   ))
  # }

  if (compare_coefs) {
    # 1. Identify all 'ofi_bin' variables and extract their numeric suffixes
    all_vars <- unique(out$var)
    bin_vars <- grep("^ofi_bin[0-9]+$", all_vars, value = TRUE)
    bin_nums <- sort(as.numeric(gsub("ofi_bin", "", bin_vars)))

    # 2. Generate all pairs (i, j) where j > i
    # combn(bin_nums, 2) creates a matrix where each column is a pair [i, j]
    pairs <- combn(bin_nums, 2)

    # 3. Use lapply to create the new data tables for each comparison
    diff_list <- apply(pairs, 2, function(p) {
      i <- p[1]
      j <- p[2]

      var_i <- paste0("ofi_bin", i)
      var_j <- paste0("ofi_bin", j)

      # Calculate the difference for all yyyymm
      # We assume 'out' is sorted or keyed so that yyyymm matches correctly
      res <- out[var == var_j, .(yyyymm, coef)]
      res[, coef := coef - out[var == var_i, coef]]
      res[, var := paste0(var_j, " - ", var_i)]

      return(res)
    })

    # 4. Bind the new comparisons back to the main table
    out <- rbindlist(c(list(out), diff_list), use.names = TRUE)
  }

  # newey-west lag
  if ("hor" %in% names(data)) {
    this_hor <- data[1, hor]
  } else {
    this_hor <- 1
  }

  coef_data <- data.table()
  for (this_var in unique(out[, var])) {
    mm <- lm(coef ~ 1, out[var == this_var])
    coef_data <- rbind(coef_data, data.table(
      var = this_var, coef = mm$coef[1],
      se = sqrt(vcov(mm)[1, 1]),
      se_nw = sqrt(NeweyWest(mm, this_hor)[1, 1])
    ))
  }
  coef_data[, r2 := r2_data[, r2]]
  coef_data[, obs := nrow(data)]
  coef_data[, type := data[1, type]]
  coef_data[, nw_lag := this_hor]

  # let's also get covariance matrix
  if (output_cov) {
    tt <- dcast(out,
      yyyymm ~ var,
      value.var = "coef"
    )[, yyyymm := NULL]
    vv <- names(tt)
    mat <- as.matrix(tt)
    fit <- lm(mat ~ 1)
    cov_matrix <- NeweyWest(fit, lag = this_hor, prewhite = FALSE)
    rownames(cov_matrix) <- vv
    colnames(cov_matrix) <- vv
  } else {
    cov_matrix <- NULL
  }

  return(list(coef_data = coef_data, cov_matrix = cov_matrix))
}

# ###########################################################
# Fama-MacBeth
# ###########################################################


# worker function to estimate regression with one type of demand variable. reg_spec = "nonlinear" or "stdev"
# TODO: take out nonlinear
# data <- data_list[[1]]
p.process_one_type <- function(data) {
  this_type <- data[1, type] # parse

  # control variables for different regression specifications
  control_formulas <- c(
    "1",
    paste0(c(controls_char), collapse = "+"),
    paste0(c(controls_char, controls_liq), collapse = "+"),
    paste0(paste0(c(controls_char, controls_liq), " * as.factor(bin)"), collapse = "+")
  )

  out_all <- data.table() # save results here

  max_spec_idx <- length(control_formulas)
  cov_all <- vector("list", max_spec_idx)

  # introduce direct controls
  for (spec_idx in 1:max_spec_idx) {
    ff_no_ofi <- paste0("ret ~ ", control_formulas[spec_idx])
    ff <- paste0("ret ~ ", control_formulas[spec_idx], " + ofi_bin1 + ofi_bin2 + ofi_bin3")
    if (this_type == "BMI") {
      ff <- paste0(ff, "+", paste0(controls_bmi, collapse = "+"))
      ff_no_ofi <- paste0(ff_no_ofi, "+", paste0(controls_bmi, collapse = "+"))
    }

    tt <- p.fama_macbeth_with_cov(data, ff, compare_coefs = TRUE, output_cov = TRUE)
    out <- tt$coef_data
    vars_select <- grepl("ofi_bin", colnames(tt$cov_matrix))
    cov_matrix <- tt$cov_matrix[vars_select, vars_select]

    out_no_ofi <- p.fama_macbeth(data, ff_no_ofi, compare_coefs = FALSE)
    out[, r2_no_ofi := out_no_ofi[1, r2]]
    out[, spec_idx := spec_idx]

    # keep track
    out_all <- rbind(out_all, out)
    cov_all[[spec_idx]] <- cov_matrix
  }

  # #  then add interactions with demand
  # for (this_v in c(controls_char, controls_liq)) {
  #   setnames(data, this_v, "xx")
  #   data[, yy := xx * ofi]
  #   setnames(data, c("xx", "yy"), c(this_v, paste0("ofi_", this_v)))
  #   ff <- paste0(ff, " + ofi_", this_v)
  #   ff_no_ofi <- paste0(ff_no_ofi, " + ofi_", this_v)
  #   spec_idx <- spec_idx + 1

  #   tt <- p.fama_macbeth_with_cov(data, ff, compare_coefs = ifelse(reg_spec == "nonlinear", FALSE, TRUE), output_cov = TRUE)
  #   out <- tt$coef_data
  #   vars <- paste0("ofi_bin", 1:3)
  #   cov_matrix <- tt$cov_matrix[vars, vars]
  #   out_no_ofi <- p.fama_macbeth(data, ff_no_ofi, compare_coefs = FALSE)

  #   out[, r2_no_ofi := out_no_ofi[1, r2]]
  #   out[, spec_idx := spec_idx]

  #   # keep track
  #   out_all <- rbind(out_all, out)
  #   cov_all[[spec_idx]] <- cov_matrix
  # }

  # # name the control variables being added
  # tmp <- data.table(
  #   spec_idx = 1:(length(controls_list) + 3),
  #   var_added = c("none_init", "controls_char", "controls_char+controls_liq", controls_list),
  #   var_type = c(
  #     rep("", 3), rep("return-predicting chars", length(controls_char)),
  #     rep("liquidity", length(controls_liq))
  #   )
  # )
  # out_all <- merge(out_all, tmp, by = "spec_idx", all.x = T)

  out_all[, type := this_type]

  # how to mark this_type in cov_all?
  names(cov_all) <- paste0(this_type, "_", names(cov_all))

  # return
  return(list(out_all = out_all, cov_all = cov_all))
}

tic("regressions")

data_list <- split(data_all, by = "type")
raw_results <- mclapply(data_list, p.process_one_type, mc.cores = nc)

# 2. Extract and stack all the data.tables (estimates)
out_coef <- rbindlist(lapply(raw_results, `[[`, "out_all"))

# 3. Extract all covariance lists and name them by 'type'
out_cov <- lapply(raw_results, `[[`, "cov_all")

toc()

to_dir <- "tmp/anna/reg_with_flexible_binning/"
dir.create(to_dir, recursive = T, showWarnings = F)
saveRDS(out_coef, paste0(to_dir, "fm_coefs.RDS"))
saveRDS(out_cov, paste0(to_dir, "fm_cov.RDS"))

# ----- just check results
options(width = 150)

out <- readRDS("tmp/anna/reg_with_flexible_binning/fm_coefs.RDS")[grepl("ofi_bin", var) & spec_idx %in% c(3, 4)]
out[, var_type := ifelse(var %in% paste0("ofi_bin", 1:3), "coef", "diff")]
out[, coef_round := round(coef, 2)]
out[, se_round := round(se, 2)]
out[, tstat_round := round(coef / se, 2)]

# rank coefs
tmp <- unique(out[, .(var)]) %>% mutate(var_lab = paste0(row_number(), "_", var))
out <- merge(out, tmp, by = "var")
rm(tmp)

dcast(out, var_lab ~ type + spec_idx, value.var = "coef_round")
dcast(out, var_lab ~ type + spec_idx, value.var = "se_round")
dcast(out, var_lab ~ type + spec_idx, value.var = "tstat_round")


this_type <- "coef"
this_type <- "diff"

dcast(out[var_type == this_type], var ~ type + spec_idx, value.var = "coef_round")
dcast(out[var_type == this_type], var ~ type + spec_idx, value.var = "se_round")
dcast(out[var_type == this_type], var ~ type + spec_idx, value.var = "tstat_round")

# # === SANITY check from earlier

# old <- readRDS("../../tmp/price_impact/regression_contemp/fm_stdev.RDS")[spec_idx %in% 1:3 & grepl("ofi_bin", var) & type != 'OFI_pre_whitened'][, c('var_type','var_added') := NULL]
# new <- readRDS("tmp/anna/reg_with_flexible_binning/fm_coefs.RDS")[spec_idx %in% 1:3 & grepl("ofi_bin", var)]
# dim(new) == dim(old)

# out <- merge(new[, .(spec_idx, var, type, coef, se, se_nw)], old[, .(spec_idx, var, type, coef, se, se_nw)], by = c("spec_idx", "var", "type"))
# out[, max(abs(coef.x - coef.y)), type]
# out[, cor(coef.x, coef.y, use = "complete.obs"), type]
# out[, cor(se.x, se.y, use = "complete.obs"), type]
