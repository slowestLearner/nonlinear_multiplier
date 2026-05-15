# ---- estimate multipliers with more bins. TODO: make sure coefs are the same for lower bins
# code taken from b_static_results/2a_regression_fm.R
library(this.path)
setwd(this.path::this.dir())
options(width = 200)
source("../../utilities/runmefirst.R")
source("../../utilities/regressions.R")

# load regression data
tic("preparing data")
data_all <- readRDS("../../tmp/raw_data/reg_inputs/reg_table_static.RDS")[type %in% c("FIT", "OFI", "BMI")]

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

# get 4 bins
data_all[, sd_ofi := sd(ofi), by = .(yyyymm, type)]
data_all[, bin := ifelse(abs(ofi) > 3 * sd_ofi, 4, ifelse(abs(ofi) > 2 * sd_ofi, 3, ifelse(abs(ofi) > sd_ofi, 2, 1)))]
data_all[, ofi_bin3 := ifelse(bin == 3, ofi, 0)]
data_all[, ofi_bin4 := ifelse(bin == 4, ofi, 0)]
toc()

# this is updated to 4 bins
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
        out <- rbind(out, data.table(
            yyyymm = yms, var = "ofi_bin4 - ofi_bin1",
            coef = out[var == "ofi_bin4", coef] - out[var == "ofi_bin1", coef]
        ))
        out <- rbind(out, data.table(
            yyyymm = yms, var = "ofi_bin4 - ofi_bin2",
            coef = out[var == "ofi_bin4", coef] - out[var == "ofi_bin2", coef]
        ))
        out <- rbind(out, data.table(
            yyyymm = yms, var = "ofi_bin4 - ofi_bin3",
            coef = out[var == "ofi_bin4", coef] - out[var == "ofi_bin3", coef]
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
control_formulas <- c(
    "1",
    paste0(c(controls_char), collapse = "+"),
    paste0(c(controls_char, controls_liq), collapse = "+")
)

# worker function to estimate regression with one type of demand variable. reg_spec = "nonlinear" or "stdev"
# just need the first 3 specifications
p.process_one_type <- function(data, reg_spec = "nonlinear") {
    this_type <- data[1, type] # parse

    out_all <- data.table() # save results here

    # max_spec_idx <- length(c(control_formulas, controls_char, controls_liq))
    max_spec_idx <- length(c(control_formulas))
    cov_all <- vector("list", max_spec_idx)
    names(cov_all) <- paste0("spec_", 1:max_spec_idx)

    # introduce direct controls
    for (spec_idx in 1:length(control_formulas)) {
        ff_no_ofi <- paste0("ret ~ ", control_formulas[spec_idx])
        if (reg_spec == "nonlinear") {
            ff <- paste0("ret ~ ", control_formulas[spec_idx], " + ofi + ofi_absofi")
        } else if (reg_spec == "stdev") {
            ff <- paste0("ret ~ ", control_formulas[spec_idx], " + ofi_bin1 + ofi_bin2 + ofi_bin3 + ofi_bin4")
        }
        if (this_type == "BMI") {
            ff <- paste0(ff, "+", paste0(controls_bmi, collapse = "+"))
            ff_no_ofi <- paste0(ff_no_ofi, "+", paste0(controls_bmi, collapse = "+"))
        }

        tt <- p.fama_macbeth_with_cov(data, ff, compare_coefs = ifelse(reg_spec == "nonlinear", FALSE, TRUE), output_cov = TRUE)
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

    # #  then add interactions with demand
    # for (this_v in c(controls_char, controls_liq)) {
    #     setnames(data, this_v, "xx")
    #     data[, yy := xx * ofi]
    #     setnames(data, c("xx", "yy"), c(this_v, paste0("ofi_", this_v)))
    #     ff <- paste0(ff, " + ofi_", this_v)
    #     ff_no_ofi <- paste0(ff_no_ofi, " + ofi_", this_v)
    #     spec_idx <- spec_idx + 1

    #     tt <- p.fama_macbeth_with_cov(data, ff, compare_coefs = ifelse(reg_spec == "nonlinear", FALSE, TRUE), output_cov = TRUE)
    #     out <- tt$coef_data
    #     vars <- paste0("ofi_bin", 1:3)
    #     cov_matrix <- tt$cov_matrix[vars, vars]
    #     out_no_ofi <- p.fama_macbeth(data, ff_no_ofi, compare_coefs = FALSE)

    #     out[, r2_no_ofi := out_no_ofi[1, r2]]
    #     out[, spec_idx := spec_idx]

    #     # keep track
    #     out_all <- rbind(out_all, out)
    #     cov_all[[spec_idx]] <- cov_matrix
    # }

    # # name the control variables being added
    # tmp <- data.table(
    #     spec_idx = 1:(length(controls_list) + 3),
    #     var_added = c("none_init", "controls_char", "controls_char+controls_liq", controls_list),
    #     var_type = c(
    #         rep("", 3), rep("return-predicting chars", length(controls_char)),
    #         rep("liquidity", length(controls_liq))
    #     )
    # )
    # out_all <- merge(out_all, tmp, by = "spec_idx", all.x = T)
    out_all[, type := this_type]

    # return
    return(list(out_all = out_all, cov_all = cov_all))
}

# stdev-based specification---
tic("regression")

# data_list <- split(data_all, by = "type")
raw_results <- mclapply(split(data_all, by = "type"), function(x) {
    p.process_one_type(x, reg_spec = "stdev")
}, mc.cores = nc)

out_stdev <- rbindlist(lapply(raw_results, `[[`, "out_all"))
toc()

# write to file
to_dir <- "../../tmp/referee/r1_minor2/regression_4bins/"
dir.create(to_dir, recursive = T, showWarnings = F)
saveRDS(out_stdev, paste0(to_dir, "fm_stdev.RDS"))
saveRDS(cov_stdev, paste0(to_dir, "fm_stdev_cov.RDS"))

# -- also save magnitude of shocks
data <- data_all[, .(d = mean(abs(ofi))), by = .(type, bin)][order(type, bin)]
saveRDS(data, paste0(to_dir, "abs_ofi_magnitude.RDS"))


# # --- sanity check: lower order coefs fully agree? yep, mostly

# old <- readRDS("../../tmp/price_impact/regression_contemp/fm_stdev.RDS")[type != "OFI_pre_whitened"]
# old <- old[spec_idx %in% 1:3 & var %in% paste0("ofi_bin", 1:3)]
# old <- old[, .(spec_idx, var, type, coef, se)]

# new <- readRDS("../../tmp/referee/r1_minor2/regression_4bins/fm_stdev.RDS")
# new <- new[spec_idx %in% 1:3 & var %in% paste0("ofi_bin", 1:3)]
# new <- new[, .(spec_idx, var, type, coef, se)]

# compare <- merge(old, new, by = c("spec_idx", "var", "type"))
# compare[, mean(abs(coef.x - coef.y)), var]
