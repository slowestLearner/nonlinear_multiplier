# Check if there is anticipatory behavior.
# Code adapted from 23_are_big_shocks_more_anticipated/3_existing_results.Rmd
# scripts Takes 2 mins to run with 6 cores
# NOTE: the controls are not lagged. I guess the correct thing to do is to do "specification-specific lags"
source("utilities/runmefirst.R")
source("utilities/regressions.R")


tic("loading data")

# All demand and returns
data <- readRDS("tmp/raw_data/reg_inputs/all_ofi_and_ret.RDS")
data <- data[yyyymm >= 199306] # liq-chars are not available for 199303
data[, ret_freq := ifelse(type == "BMI", "monthly", "quarterly")]
data[, ret := NULL]

# get returns with more lags
tmp <- readRDS("../../data/stocks/prices/quarterly_return.RDS")[, .(yyyymm, permno, ret_freq = "quarterly", ret)]
tmp <- rbind(tmp, readRDS("../../data/stocks/prices/monthly_return.RDS")[, .(yyyymm, permno, ret_freq = "monthly", ret)])
tmp[, idx := frank(yyyymm, ties.method = "dense"), ret_freq]
for (i in 1:4) {
    tmp <- merge(tmp, tmp[, .(idx = idx + i, ret_freq, permno, xx = ret)], by = c("idx", "ret_freq", "permno"), all.x = T)
    setnames(tmp, "xx", paste0("ret_", i))
}
tmp[, idx := NULL]
out <- tmp[, .(yyyymm, permno, ret_freq, lag = 0, ret)]
out <- rbind(out, tmp[, .(yyyymm, permno, ret_freq, lag = 1, ret = ret + ret_1)])
out <- rbind(out, tmp[, .(yyyymm, permno, ret_freq, lag = 2, ret = ret + ret_1 + ret_2)])
out <- rbind(out, tmp[, .(yyyymm, permno, ret_freq, lag = 3, ret = ret + ret_1 + ret_2 + ret_3)])
out <- rbind(out, tmp[, .(yyyymm, permno, ret_freq, lag = 4, ret = ret + ret_1 + ret_2 + ret_3 + ret_4)])
rm(tmp)
data <- merge(data, out, by = c("yyyymm", "permno", "ret_freq"), all.x = T, allow.cartesian = T) %>% na.omit()
rm(out)

# add controls that are specific to BMI
tmp <- readRDS("tmp/raw_data/controls/controls_for_BMI.RDS")
controls_bmi <- setdiff(names(tmp), c("yyyymm", "permno"))
data <- merge(data, tmp, by = c("yyyymm", "permno"), all.x = T)
rm(tmp)

# get control variable names
cdata <- readRDS("../../formal_tests/20250117_quarterly/tmp/raw_data/controls/controls_classification.RDS")
controls_liq <- cdata[control_type == "liquidity", var]
controls_char <- cdata[control_type == "return-predictor", var]
rm(cdata)

# add the other controls
tmp <- readRDS("../../formal_tests/20250117_quarterly/tmp/raw_data/controls/quarterly_controls_lagged.RDS")
tmp <- tmp[, c("yyyymm", "permno", controls_liq, controls_char), with = F]
data <- merge(data, tmp, by = c("yyyymm", "permno"), all.x = T)
rm(tmp)

# redefine type as type-lag
data[, type := paste0(type, "_", lag)]
data[, lag := NULL]

# Re-center characteristics within each time period and demand type
for (this_v in c(controls_liq, controls_char)) {
    # print(this_v)
    setnames(data, this_v, "xx")
    data <- merge(data, data[, .(m = mean(xx, na.rm = T)), .(yyyymm, type)], by = c("yyyymm", "type"))
    data[, xx := xx - m]
    data[, m := NULL]
    setnames(data, "xx", this_v)
}
rm(this_v)
data[is.na(data)] <- 0 # fill zeros if needed (for these characteristics)

# compute the RHS variables
data[, ofi_absofi := ofi * abs(ofi)]
data[, sd_ofi := sd(ofi), .(yyyymm, type)]
data[, bin := 1]
data[abs(ofi) > sd_ofi, bin := 2]
data[abs(ofi) > 2 * sd_ofi, bin := 3]
data[, ofi_bin1 := ifelse(bin == 1, ofi, 0)]
data[, ofi_bin2 := ifelse(bin == 2, ofi, 0)]
data[, ofi_bin3 := ifelse(bin == 3, ofi, 0)]
data[, sd_ofi := NULL]

controls_list <- c(controls_char, controls_liq)

# also specify lags for nw
data[, hor := as.integer(substr(type, 5, 5)) + 1] # this is used as NW lag
data_all <- copy(data)
rm(data)
toc()

# ###########################################################
# Fama-MacBeth
# ###########################################################

# control variables for different specifications
control_formulas <- c(
    "1",
    paste0(c(controls_char), collapse = "+"),
    paste0(c(controls_char, controls_liq), collapse = "+")
)

# process one type of data. reg_spec = "nonlinear" or "stdev"
# it is a bit wasteful as it does a lot, but can fix later
p.process_one_type <- function(data, reg_spec = "nonlinear") {
    this_type <- data[1, type]

    out_all <- data.table() # save all results here

    # introduce direct controls
    for (spec_idx in 1:length(control_formulas)) {
        ff_no_ofi <- paste0("ret ~ ", control_formulas[spec_idx])
        if (reg_spec == "nonlinear") {
            ff <- paste0("ret ~ ", control_formulas[spec_idx], " + ofi + ofi_absofi")
        } else if (reg_spec == "stdev") {
            ff <- paste0("ret ~ ", control_formulas[spec_idx], " + ofi_bin1 + ofi_bin2 + ofi_bin3")
        }
        if (this_type == "BMI") {
            ff <- paste0(ff, "+", paste0(controls_bmi, collapse = "+"))
            ff_no_ofi <- paste0(ff_no_ofi, "+", paste0(controls_bmi, collapse = "+"))
        }
        out <- p.fama_macbeth(data, ff, compare_coefs = ifelse(reg_spec == "nonlinear", FALSE, TRUE))
        out_no_ofi <- p.fama_macbeth(data, ff_no_ofi, compare_coefs = FALSE)

        out[, r2_no_ofi := out_no_ofi[1, r2]]
        out[, spec_idx := spec_idx]
        out_all <- rbind(out_all, out)
    }

    # then add interactions with ofi
    for (this_v in c(controls_char, controls_liq)) {
        setnames(data, this_v, "xx")
        data[, yy := xx * ofi]
        setnames(data, c("xx", "yy"), c(this_v, paste0("ofi_", this_v)))
        ff <- paste0(ff, " + ofi_", this_v)
        ff_no_ofi <- paste0(ff_no_ofi, " + ofi_", this_v)
        spec_idx <- spec_idx + 1

        out <- p.fama_macbeth(data, ff, compare_coefs = ifelse(reg_spec == "nonlinear", FALSE, TRUE))
        out_no_ofi <- p.fama_macbeth(data, ff_no_ofi, compare_coefs = FALSE)
        out[, r2_no_ofi := out_no_ofi[1, r2]]
        out[, spec_idx := spec_idx]
        out_all <- rbind(out_all, out)
    }

    # name the control variables being added
    tmp <- data.table(
        spec_idx = 1:(length(controls_list) + 3),
        var_added = c("none_init", "controls_char", "controls_char+controls_liq", controls_list),
        var_type = c(
            rep("", 3), rep("return-predicting chars", length(controls_char)),
            rep("liquidity", length(controls_liq))
        )
    )
    out_all <- merge(out_all, tmp, by = "spec_idx", all.x = T)
    out_all[, type := this_type]

    return(out_all)
}

# stdev-based specification---
tic("estimating stdev-based specification")
out_stdev <- rbindlist(mclapply(split(data_all, by = "type"), function(x) {
    p.process_one_type(x, reg_spec = "stdev")
}, mc.cores = nc))
toc()

dir.create("tmp/price_impact/contemp/", recursive = T, showWarnings = F)
saveRDS(out_stdev, "tmp/price_impact/regression_contemp/fm_stdev_anticipation_same_controls.RDS")
