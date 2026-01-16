# Compare rolling standard deviations of FIT and OFI
library(this.path)
setwd(this.path::this.dir())
source("../utilities/runmefirst.R")
source("../utilities/regressions.R")

# load regression data
data <- readRDS("../tmp/raw_data/reg_inputs/reg_table_static.RDS") %>% filter(type != "BMI")

# get control variable names
cdata <- readRDS("../tmp/raw_data/controls/controls_classification.RDS")
controls_char <- cdata[control_type == "return-predictor", var]
controls_liq <- cdata[control_type == "liquidity", var]
controls_list <- c(controls_char, controls_liq)
rm(cdata)

# merge with lagged stdev
demand_vol <- readRDS("../tmp/raw_data/reg_inputs/quarterly_fit_and_ofi_lagged_rolling_stdev.RDS")
demand_vol <- demand_vol[0 == rowSums(is.na(demand_vol))]
demand_vol[, std_ofi := Winsorize(std_ofi, val = quantile(std_ofi, probs = c(.01, .99))), .(type, spec)]

# translate lags
demand_vol[, spec := gsub("sd_ofi_", "", spec)]
demand_vol[, spec := as.integer(gsub("q", "", spec))]
setnames(demand_vol, "spec", "stdev_lag")
data <- merge(data, demand_vol, by = c("yyyymm", "permno", "type"), allow.cartesian = T)
rm(demand_vol)

# standardize demand
setnames(data, "ofi", "ofi_raw")
data[, ofi := ofi_raw / std_ofi]
data[, ofi := Winsorize(ofi, val = quantile(ofi, probs = c(.001, .999))), .(type, stdev_lag)]
data[, scaling := sd(ofi_raw) / sd(ofi), .(type, stdev_lag)]
data[, ofi := ofi * scaling]
data[, c("scaling", "ofi_raw", "std_ofi") := NULL]

# mark specifications
data[, type := paste0(type, "-", stdev_lag)]

# we need to redo the transformation to get nonlinear x variables AFTER dividing by vol
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
data_all <- copy(data)
rm(data)
gc()

# ###########################################################
# Fama-MacBeth with nonlinear or stdev-based specification
# ###########################################################

# control variables for different specifications
control_formulas <- c(
    "1",
    paste0(c(controls_char), collapse = "+"),
    paste0(c(controls_char, controls_liq), collapse = "+")
)

# process one type of data. reg_spec = "nonlinear" or "stdev"
p.process_one_type <- function(data, reg_spec = "nonlinear") {
    this_type <- data[1, type]

    out_all <- data.table() # save all results here

    # introduce direct controls
    for (spec_idx in 1:length(control_formulas)) {
        if (reg_spec == "nonlinear") {
            ff <- paste0("ret ~ ", control_formulas[spec_idx], " + ofi + ofi_absofi")
        } else if (reg_spec == "stdev") {
            ff <- paste0("ret ~ ", control_formulas[spec_idx], " + ofi_bin1 + ofi_bin2 + ofi_bin3")
        }
        if (this_type == "BMI") {
            ff <- paste0(ff, "+", paste0(controls_bmi, collapse = "+"))
        }
        out <- p.fama_macbeth(data, ff, compare_coefs = ifelse(reg_spec == "nonlinear", FALSE, TRUE))
        out[, spec_idx := spec_idx]
        out_all <- rbind(out_all, out)
    }

    # then add interactions with ofi
    for (this_v in c(controls_char, controls_liq)) {
        setnames(data, this_v, "xx")
        data[, yy := xx * ofi]
        setnames(data, c("xx", "yy"), c(this_v, paste0("ofi_", this_v)))
        ff <- paste0(ff, " + ofi_", this_v)
        spec_idx <- spec_idx + 1

        out <- p.fama_macbeth(data, ff, compare_coefs = ifelse(reg_spec == "nonlinear", FALSE, TRUE))
        out[, spec_idx := spec_idx]
        out_all <- rbind(out_all, out)
    }

    # name the control variables being added
    tmp <- data.table(
        spec_idx = 1:(length(controls_list) + 4),
        var_added = c("none_init", "controls_char", "controls_char+controls_liq", "none", controls_list),
        var_type = c(
            rep("", 4), rep("return-predicting chars", length(controls_char)),
            rep("liquidity", length(controls_liq))
        )
    )
    out_all <- merge(out_all, tmp, by = "spec_idx", all.x = T)
    out_all[, type := this_type]

    return(out_all)
}


# in parallel, takes around X mins
tic()
out_stdev <- rbindlist(mclapply(split(data_all, by = "type"), function(x) {
    p.process_one_type(x, reg_spec = "stdev")
}, mc.cores = nc))
gc()
toc()

to_dir <- "../tmp/price_impact/regression_contemp/"
dir.create(to_dir, recursive = T, showWarnings = F)
saveRDS(out_stdev, paste0(to_dir, "standardized_fm_stdev.RDS"))


# # --- SANITY check


# # --- stdev
# new <- readRDS("tmp/price_impact/regression_contemp/standardized_fm_stdev.RDS")
# new[, stdev_lag := as.integer(gsub(".*-", "", type))]
# new[, type := substr(type, 1, 3)]
# old <- readRDS("../20250117_quarterly/tmp/price_impact/regression_contemp/full_sample_standardized_d_by_stdev_bin.RDS")
# old <- old[type != "OFI"]
# old[type == "OFI_resid", type := "OFI"]

# # get joint variables
# new <- new[spec_idx == 3][, spec_idx := NULL]
# old[, var := gsub("M", "ofi_bin", var)]
# new[, var := gsub(" ", "", var)]
# vv <- intersect(names(new), names(old))
# new <- new[, ..vv]
# old <- old[, ..vv]
# rm(vv)

# common_vars <- intersect(unique(new[, var]), unique(old[, var]))
# new <- new[var %in% common_vars]
# old <- old[var %in% common_vars]
# rm(common_vars)
# stopifnot(dim(new) == dim(old))

# # compare - not identical, but similar?
# compare <- merge(new, old, by = c("type", "var", "stdev_lag"), all = T)
# rm(new, old)

# compare[, mean(abs(coef.x - coef.y)) / mean(abs(coef.x))]
# compare[, mean(abs(se.x - se.y)) / mean(abs(se.x))]
# compare[, mean(abs(obs.x - obs.y)) / mean(abs(obs.x))]
# compare[, mean(abs(r2.x - r2.y)) / mean(abs(r2.x))]

# compare[, cor(coef.x, coef.y)]
# compare[, cor(se.x, se.y)]
# compare[, cor(obs.x, obs.y)]
# compare[, cor(r2.x, r2.y)]
