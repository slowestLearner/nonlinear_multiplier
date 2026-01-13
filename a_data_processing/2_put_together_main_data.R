# === Put together the monthly/quarterly data for main regressions

# load libraries
library(this.path)
setwd(this.path::this.dir())
source("../utilities/runmefirst.R")

# == 1) Put together quarterly OFI/FIT and return data.

# quarterly stock returns
stock_data <- readRDS("../../../data/stocks/prices/quarterly_return.RDS")
stock_data <- stock_data[, list(yyyymm, permno, ret)]

# FIT. Pick a cleaned version
data <- readRDS("../../../data/demand_shocks/j_fit/quarterly_updated_20250313.RDS")[, list(yyyymm, permno, type = "FIT", ofi = fit_adj_cut01)]

# pre-whitened OFI
tmp <- readRDS("../tmp/raw_data/cleaning/ofi/residual_daily_ofi.RDS")
tmp[, ofi_resid := DescTools::Winsorize(ofi_resid, val = quantile(ofi_resid, probs = c(.005, .995)))]
tmp[, yyyymm := 100 * year(date) + 3 * quarter(date)]
tmp <- tmp[, .(obs = length(date), type = "OFI_pre_whitened", ofi = sum(ofi_resid)), .(yyyymm, permno)]
tmp[, max_obs := max(obs), yyyymm]
tmp <- tmp[obs >= (max_obs - 5)] # need to have data for most days in the quarter
tmp[, c("obs", "max_obs") := NULL]
data <- rbind(data, tmp)
rm(tmp)

# non-pre-whitened OFI
tmp <- readRDS("../../../data/demand_shocks/ofi/quarterly.RDS")[, type := "OFI"]
data <- rbind(data, tmp)
rm(tmp)

# merge with stock data
data <- merge(stock_data, data, by = c("yyyymm", "permno"), allow.cartesian = T)
rm(stock_data)
data_fit_and_ofi <- copy(data)
rm(data)


# === 2) add BMI to the quarterly OFI/FIT tables.

# TODO: move it to upstrea
# data <- data.table(haven::read_stata("../../../tests/9_russell_rdd/russell_rdd_stock_level_panel.dta"))
data <- readRDS("../../../data/demand_shocks/bmi/russell_rdd_stock_level_panel.RDS")

# Calculate the log of tot_mktcap_r3
data$log_tot_mktcap_r3 <- log(data$tot_mktcap_r3)

# Initialize 'band' to 0
data$band <- 0

# Set 'band' to 1 where conditions are met
data$band[data$may_mkt_cap_ranks > data$upper_cutoff_rank &
    data$may_mkt_cap_ranks < data$lower_cutoff_rank] <- 1

# Initialize 'band_d2' to 0 and then compute its value
data$band_d2 <- 0
data$band_d2 <- data$band * data$was_in_2000_in_may

names(data)[names(data) == "BA_Spread_Percentage_1_Year_MA"] <- "ba_spread"
PS_controls <- c("log_tot_mktcap_r3", "band", "was_in_2000_in_may", "band_d2", "ba_spread")

# set up the band width filter
bandwidth <- 150
data[, filter := 0]
data[((abs(dist_up) <= bandwidth) | (abs(dist_down) <= bandwidth)), filter := 1]

# compute changes in BMI. bmi is in May, bmi_june is in June. Reccale
data[, d_bmi := bmi_june - bmi]
data[, d_bmi := .2 * d_bmi]
data <- data[!is.na(d_bmi) & !is.na(full_ret)]

# just need these
data <- data[filter == 1]
data[, filter := NULL]

# change into the same format
setnames(data, c("full_ret", "d_bmi"), c("ret", "ofi"))
data[, year := 100 * year + 6]
setnames(data, "year", "yyyymm")

# combine with FIT and OFI data
reg_data <- rbind(data[, .(yyyymm, permno, type = "BMI", ret, ofi)], data_fit_and_ofi)
to_dir <- "../tmp/raw_data/reg_inputs/"
dir.create(to_dir, recursive = T, showWarnings = F)
saveRDS(reg_data, paste0(to_dir, "/all_ofi_and_ret.RDS"))
rm(reg_data)

# also save BMI-specific controls
data <- data[, c("yyyymm", "permno", PS_controls), with = F]
to_dir <- "../tmp/raw_data/controls/"
dir.create(to_dir, recursive = T, showWarnings = F)
saveRDS(data, paste0(to_dir, "/controls_for_BMI.RDS"))

# # === SANITY: check with earlier data

# # BMI controls
# old = readRDS('../20250117_quarterly/tmp/raw_data/controls/controls_for_BMI.RDS')
# new = readRDS('tmp/raw_data/controls/controls_for_bmi.RDS')
# mean(old == new, na.rm = T)

# # demand data
# old = readRDS('../20250117_quarterly/tmp/raw_data/reg_inputs/all_ofi_and_ret.RDS')
# old = old[type != 'OFI']
# old[type == 'OFI_resid', type := 'OFI']
# new = readRDS('tmp/raw_data/reg_inputs/all_ofi_and_ret.RDS')
# dim(new) == dim(old)
# compare = merge(old, new, by = c('yyyymm','permno','type')); rm(old,new)
# compare[, cor(ofi.x, ofi.y), type]
