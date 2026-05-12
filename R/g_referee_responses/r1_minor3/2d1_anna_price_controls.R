# --- Answers to Anna's question. Let's add the controls one at a time
# code directly taken from the main static specification
library(this.path)
setwd(this.path::this.dir())
source("../../utilities/runmefirst.R")
source("../../utilities/regressions.R")
options(width = 120)
library(patchwork)

# load regression data
tic("preparing data")
data_all <- readRDS("../../tmp/raw_data/reg_inputs/reg_table_static.RDS")[type != "OFI_pre_whitened"]

# get control variable names
# cdata <- readRDS("../../tmp/raw_data/controls/controls_classification.RDS")
cdata <- fread("tmp/anna/control_classification.csv")
cdata[is.na(cdata)] <- 0

# Let's rank the order in which they are actually added
cdata <- cdata[order(-None, -Indirect, -Direct)]
# controls_char <- cdata[control_type == "return-predictor", var]
# controls_liq <- cdata[control_type == "liquidity", var]
# controls_list <- c(controls_char, controls_liq)

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

  # compare coefficient differences
  if (compare_coefs) {
    out <- rbind(out, data.table(
      yyyymm = yms, var = "ofi_bin2 - ofi_bin1",
      coef = out[var == "ofi_bin2", coef] - out[var == "ofi_bin1", coef]
    ))
    out <- rbind(out, data.table(
      yyyymm = yms, var = "ofi_bin3 - ofi_bin1",
      coef = out[var == "ofi_bin3", coef] - out[var == "ofi_bin1", coef]
    ))
    out <- rbind(out, data.table(
      yyyymm = yms, var = "ofi_bin3 - ofi_bin2",
      coef = out[var == "ofi_bin3", coef] - out[var == "ofi_bin2", coef]
    ))
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
    tt <- dcast(out
      # [var %in% c("ofi_bin1", "ofi_bin2", "ofi_bin3")]
      , yyyymm ~ var,
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

# control variables for different regression specifications
spec_data <- data.table(var_added = c("1", cdata[, var]))
spec_data[, spec_idx := 1:nrow(spec_data)]
spec_data[1, ff := "1"]
for (i in 2:nrow(spec_data)) {
  spec_data[i, ff := paste0(spec_data[2:i, var_added], collapse = " + ")]
}

# worker function to estimate regression with one type of demand variable. reg_spec = "nonlinear" or "stdev"
# TODO: take out nonlinear
# data <- data_list[[1]]
p.process_one_type <- function(data) {
  this_type <- data[1, type] # parse

  out_all <- data.table() # save results here

  max_spec_idx <- nrow(spec_data)
  cov_all <- vector("list", max_spec_idx)
  names(cov_all) <- paste0("spec_", 1:max_spec_idx)

  # introduce direct controls
  for (spec_idx in spec_data[, spec_idx]) {
    # spec_idx <- 1
    ff_no_ofi <- paste0("ret ~ ", spec_data[spec_idx, ff])
    ff <- paste0("ret ~ ", spec_data[spec_idx, ff], " + ofi_bin1 + ofi_bin2 + ofi_bin3")
    if (this_type == "BMI") {
      ff <- paste0(ff, "+", paste0(controls_bmi, collapse = "+"))
      ff_no_ofi <- paste0(ff_no_ofi, "+", paste0(controls_bmi, collapse = "+"))
    }

    tt <- p.fama_macbeth_with_cov(data, ff, compare_coefs = TRUE, output_cov = TRUE)
    out <- tt$coef_data
    vars <- paste0("ofi_bin", 1:3)
    cov_matrix <- tt$cov_matrix[vars, vars]
    out_no_ofi <- p.fama_macbeth(data, ff_no_ofi, compare_coefs = FALSE)

    out[, r2_no_ofi := out_no_ofi[1, r2]]
    out[, spec_idx := spec_idx]

    # keep track
    out_all <- rbind(out_all, out)
    cov_all[[spec_idx]] <- cov_matrix
  }

  out_all <- merge(out_all, spec_data[, .(spec_idx, var_added)], by = "spec_idx", all.x = T)
  out_all[, type := this_type]

  return(list(out_all = out_all, cov_all = cov_all))
}

# stdev-based specification---
tic("anna regression")

data_list <- split(data_all, by = "type")
raw_results <- mclapply(data_list, p.process_one_type, mc.cores = nc)

# 2. Extract and stack all the data.tables (estimates)
# This works because we are specifically grabbing the 'out_all' element from each list
out_stdev <- rbindlist(lapply(raw_results, `[[`, "out_all"))

# 3. Extract all covariance lists and name them by 'type'
# This creates a nested list: cov_stdev$BMI$spec_1, cov_stdev$OtherType$spec_1
cov_stdev <- lapply(raw_results, `[[`, "cov_all")
toc()

to_dir <- "tmp/anna/regression_results/"
dir.create(to_dir, recursive = T, showWarnings = F)
saveRDS(out_stdev, paste0(to_dir, "fm_stdev.RDS"))
saveRDS(cov_stdev, paste0(to_dir, "fm_stdev_cov.RDS"))

# --- let's just plot

out_all <- readRDS("tmp/anna/regression_results/fm_stdev.RDS")
out_all[, var_type := ifelse(var %in% paste0("ofi_bin", 1:3), "coef", "diff")]
out_all[, spec_lab := var_added]
out_all[var_add == "1", spec_lab := ""]

spec_cuts <- c(8.5, 12.5)

this_type <- "FIT"
this_type <- "OFI"
this_type <- "BMI"

out <- copy(out_all[var_type == "coef" & type == this_type]) %>% setorder(var, spec_idx)
p1 <- ggplot(
  out,
  aes(x = spec_idx, y = coef, color = var)
) +
  geom_line(lwd = 2) +
  geom_point(cex = 5) +
  geom_ribbon(aes(ymin = coef - se, ymax = coef + se, fill = var), alpha = 0.2, color = NA) +
  theme_classic() +
  scale_x_continuous(breaks = out[, spec_idx], labels = out[, spec_lab]) +
  labs(x = element_blank(), y = "Coefficient") +
  geom_hline(yintercept = 0, lty = 3) +
  theme(text = element_text(size = 25), legend.position = "bottom", legend.title = element_blank(), axis.text.x = element_text(angle = 45, vjust = .6)) +
  ggtitle(paste0(this_type, " coefs")) +
  geom_vline(xintercept = spec_cuts, lty = 3)

# Plot 2: Differences
out <- copy(out_all[var_type == "diff" & type == this_type]) %>% setorder(var, spec_idx)
p2 <- ggplot(
  out,
  aes(x = spec_idx, y = coef, color = var)
) +
  geom_line(lwd = 2) +
  geom_point(cex = 5) +
  geom_ribbon(aes(ymin = coef - se, ymax = coef + se, fill = var), alpha = 0.2, color = NA) +
  theme_classic() +
  scale_x_continuous(breaks = out[, spec_idx], labels = out[, spec_lab]) +
  labs(x = element_blank(), y = "Coefficient") +
  geom_hline(yintercept = 0, lty = 3) +
  theme(text = element_text(size = 25), legend.position = "bottom", legend.title = element_blank(), axis.text.x = element_text(angle = 45, vjust = .6)) +
  ggtitle(paste0(this_type, " diffs")) +
  geom_vline(xintercept = spec_cuts, lty = 3)

# plot
p1 / p2

# # === SANITY check. quite similar

# # nonlinear
# new <- readRDS("../tmp/price_impact/regression_contemp/fm_stdev.RDS")
# old <- readRDS("../tmp/price_impact/regression_contemp_todel/fm_stdev.RDS")
# dim(new) == dim(old)

# out <- merge(new[, .(spec_idx, var, type, coef, se, se_nw)], old[, .(spec_idx, var, type, coef, se, se_nw)], by = c("spec_idx", "var", "type"))
# out[, cor(coef.x, coef.y, use = "complete.obs"), type]
# out[, cor(se.x, se.y, use = "complete.obs"), type]
