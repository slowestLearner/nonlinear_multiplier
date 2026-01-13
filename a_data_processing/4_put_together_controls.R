# === script to put together controls for regressions
library(this.path)
setwd(this.path::this.dir())
source("utilities/runmefirst.R")


# liquidity chars
data <- readRDS("../../../data/stocks/controls/monthly_liquidity_measures_not_lagged.RDS")
setnames(data, "me", "size") # to not overlap with size in characteristics
vv_liq <- setdiff(names(data), c("yyyymm", "permno")) # keep track of liquidity variable names

# regular characteristics
tmp <- readRDS("../../../data/stocks/controls/monthly_characteristics_not_lagged.RDS")
tmp[, c("size", "realized_vol") := NULL] # already part of liquidity chars
vv_char <- setdiff(names(tmp), c("yyyymm", "permno"))
data <- merge(data, tmp, by = c("yyyymm", "permno"), all.x = T)
rm(tmp)

# industry classifications. TODO: delete these later
tmp <- readRDS("../../../data/stocks/controls/ff12_industries_zero_mean.RDS")
vv_ind <- setdiff(names(tmp), c("yyyymm", "permno"))
data <- merge(data, tmp, by = c("yyyymm", "permno"), all.x = T)
rm(tmp)
data_all <- copy(data)

# === monthly version, lagged
data <- copy(data_all)
time_data <- unique(data[, list(yyyymm)])
time_data[, yyyymm_next := data.table::shift(yyyymm, -1)]
data <- merge(data, time_data, by = "yyyymm")
rm(time_data)
data[, yyyymm := yyyymm_next][, yyyymm_next := NULL]
dir.create("../tmp/raw_data/controls/", recursive = T, showWarnings = F)
saveRDS(data, "../tmp/raw_data/controls/monthly_controls_lagged.RDS")

# === quarterly version, lagged
data <- copy(data_all)
time_data <- unique(data[, .(yyyymm)])
time_data[, mm := yyyymm - 100 * floor(yyyymm / 100)]
time_data <- time_data[mm %in% c(3, 6, 9, 12)]
time_data[, yyyymm_next := data.table::shift(yyyymm, -1)]
time_data[, mm := NULL]
data <- merge(data, time_data, by = "yyyymm")
rm(time_data)
data[, yyyymm := yyyymm_next][, yyyymm_next := NULL]
saveRDS(data, "../tmp/raw_data/controls/quarterly_controls_lagged.RDS")

# === save the list of control names
tmp <- data.table(
    control_type = c(
        rep("liquidity", length(vv_liq)),
        rep("industry", length(vv_ind)),
        rep("return-predictor", length(vv_char))
    ),
    var = c(vv_liq, vv_ind, vv_char)
)

# give human-readable names
tmp[control_type == "liquidity", var_lab := c(
    "Effective Spread", "Quoted Spread", "Realized Vol",
    "Size", "Turnover", "Volume"
)]
tmp[control_type == "industry", var_lab := c(
    "Manuf", "Chemical", "Telecom", "Utilities", "Bus Equip", "Nondurable",
    "Shops", "Health", "Money", "Durable", "Energy"
)]
tmp[control_type == "return-predictor", var_lab := c(
    "Accruals", "Asset gr", "Beta", "B/M", "Profitability", "Ind Mom", "Intermediate Mom", "1y issuance", "5y issuance", "Mom", "Seasonal Mom",
    "NOA", "Reversal"
)]
saveRDS(tmp, "../tmp/raw_data/controls/controls_classification.RDS")


# # === SANITY CHECK: check with earlier data

# # monthly
# old = readRDS('../20250117_quarterly/tmp/raw_data/controls/monthly_controls_lagged.RDS')
# new = readRDS('tmp/raw_data/controls/monthly_controls_lagged.RDS')
# mean(old == new, na.rm = T)
# mean(is.na(old) == is.na(new))

# # quarterly
# old = readRDS('../20250117_quarterly/tmp/raw_data/controls/quarterly_controls_lagged.RDS')
# new = readRDS('tmp/raw_data/controls/quarterly_controls_lagged.RDS')
# mean(old == new, na.rm = T)
# mean(is.na(old) == is.na(new))

# # compare the list of control names
# old = readRDS('../20250117_quarterly/tmp/raw_data/controls/controls_classification.RDS')
# new = readRDS('tmp/raw_data/controls/controls_classification.RDS')
# mean(old == new)
