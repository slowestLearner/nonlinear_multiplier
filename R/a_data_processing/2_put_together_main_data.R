# --- put together quarterly data for regressions (1993 to 2022)
library(this.path)
setwd(this.path::this.dir())
source("../utilities/runmefirst.R")

# --- 1) combine OFI, FIT, and returns

# stock returns
stock_data <- readRDS("../../../../data/stocks/prices/quarterly_return.RDS")[, .(yyyymm, permno, ret)]

# FIT
data <- readRDS("../../../../data/demand_shocks/j_fit/quarterly.RDS")[, .(yyyymm, permno, type = "FIT", ofi = fit_adj_cut01)]

# pre-whitened daily OFI, aggregate to quarterly
tmp <- readRDS("../tmp/raw_data/cleaning/ofi/residual_daily_ofi.RDS")
# tmp[, ofi_resid := DescTools::Winsorize(ofi_resid, val = quantile(ofi_resid, probs = c(.005, .995)))]
tmp[, yyyymm := 100 * year(date) + 3 * quarter(date)]
tmp <- tmp[, .(obs = length(date), type = "OFI_pre_whitened", ofi = sum(ofi_resid)), .(yyyymm, permno)]
tmp[, max_obs := max(obs), yyyymm]
tmp <- tmp[obs >= (max_obs - 5)][, c("obs", "max_obs") := NULL] # need to have data for most days in the quarter
data <- rbind(data, tmp)
rm(tmp)

# original OFI
tmp <- readRDS("../../../../data/demand_shocks/ofi/quarterly.RDS")[, type := "OFI"]
data <- rbind(data, tmp)
rm(tmp)

# merge with stock data
data_fit_and_ofi <- merge(stock_data, data, by = c("yyyymm", "permno"), allow.cartesian = T)
rm(stock_data, data)

# -- 2) add BMI to the quarterly OFI/FIT tables. keep track of associated controls as well

# data <- data.table(haven::read_stata("../../../tests/9_russell_rdd/russell_rdd_stock_level_panel.dta"))
data <- readRDS("../../../../data/demand_shocks/bmi/russell_rdd_stock_level_panel.RDS")

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

# compute changes in BMI. bmi is in May, bmi_june is in June. Rescale, following the original BMI paper
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

# take out non-unique entries
data[, idx := 1:.N, .(yyyymm, permno)]
data <- data[idx == 1][, idx := NULL]

# combine with FIT and OFI data
reg_data <- rbind(data[, .(yyyymm, permno, type = "BMI", ret, ofi)], data_fit_and_ofi)
reg_data <- reg_data[yyyymm >= 199306 & yyyymm <= 202212] # the main sample period for the paper
to_file <- "../tmp/raw_data/reg_inputs/all_ofi_and_ret.RDS"
dir.create(dirname(to_file), showWarnings = FALSE, recursive = TRUE)
saveRDS(reg_data, to_file)

# also save BMI-specific controls
data <- data[, c("yyyymm", "permno", PS_controls), with = F]

to_file <- "../tmp/raw_data/controls/controls_for_BMI.RDS"
dir.create(dirname(to_file), showWarnings = FALSE, recursive = TRUE)
saveRDS(data, to_file)

# # --- SANITY: quite similar to earlier data

# # BMI controls
# old <- readRDS("../../../20250117_quarterly/tmp/raw_data/controls/controls_for_BMI.RDS")
# new <- readRDS("../tmp/raw_data/controls/controls_for_bmi.RDS")

# # keep unique ones
# old[, obs := .N, .(yyyymm, permno)]
# old <- old[obs == 1][, obs := NULL]
# new <- merge(new, old[, .(yyyymm, permno)], by = c("yyyymm", "permno"))
# setkey(old, yyyymm, permno)
# setkey(new, yyyymm, permno)
# mean(old == new, na.rm = T)

# # demand data
# old <- readRDS("../tmp/raw_data/reg_inputs_todel/all_ofi_and_ret.RDS")
# new <- readRDS("../tmp/raw_data/reg_inputs/all_ofi_and_ret.RDS")
# out <- merge(old, new, by = c("yyyymm", "permno", "type"), all = T)

# out <- out %>% na.omit() # only old, OFI-pre-whitened may be missing
# out[, cor(ofi.x, ofi.y), type] # only OFI-pre-whitened somewhat changed, and it does not matter
