# --- combine various controls that are later used in regressions
library(this.path)
setwd(this.path::this.dir())
source("utilities/runmefirst.R")

# --- 1) get data, save quarterly lagged controls

# liquidity char
data <- readRDS("../../../../data/stocks/controls/monthly_liquidity_measures_not_lagged.RDS")
setnames(data, "me", "size") # to not overlap with size in characteristics
vv_liq <- setdiff(names(data), c("yyyymm", "permno")) # keep track of liquidity variable names

# regular characteristics
tmp <- readRDS("../../../../data/stocks/controls/monthly_characteristics_not_lagged.RDS")
tmp[, c("size", "realized_vol") := NULL] # already part of liquidity chars
vv_char <- setdiff(names(tmp), c("yyyymm", "permno"))
data <- merge(data, tmp, by = c("yyyymm", "permno"), all.x = T)
rm(tmp)

# industry classifications
tmp <- readRDS("../../../../data/stocks/controls/ff12_industries_zero_mean.RDS")
vv_ind <- setdiff(names(tmp), c("yyyymm", "permno"))
data <- merge(data, tmp, by = c("yyyymm", "permno"), all.x = T)
rm(tmp)
# data_all <- copy(data)

# # === monthly version, lagged
# data <- copy(data_all)
# time_data <- unique(data[, list(yyyymm)])
# time_data[, yyyymm_next := data.table::shift(yyyymm, -1)]
# data <- merge(data, time_data, by = "yyyymm")
# rm(time_data)
# data[, yyyymm := yyyymm_next][, yyyymm_next := NULL]
# dir.create("../tmp/raw_data/controls/", recursive = T, showWarnings = F)
# saveRDS(data, "../tmp/raw_data/controls/monthly_controls_lagged.RDS")

# put into lagged quarterly controls
time_data <- unique(data[, .(yyyymm)])[order(yyyymm)]
time_data[, mm := yyyymm - 100 * floor(yyyymm / 100)]
time_data <- time_data[mm %in% c(3, 6, 9, 12)]
time_data[, yyyymm_next := data.table::shift(yyyymm, -1)]
time_data[, mm := NULL]
data <- merge(data, time_data, by = "yyyymm")
rm(time_data)
data[, yyyymm := yyyymm_next][, yyyymm_next := NULL]
data <- data[!is.na(yyyymm)]
to_file <- "../tmp/raw_data/controls/quarterly_controls_lagged.RDS"
dir.create(dirname(to_file), recursive = T, showWarnings = F)
saveRDS(data, to_file)

# # --- sanity check, good, cannot recheck

# tmp <- readRDS("../tmp/raw_data/controls/quarterly_controls_lagged.RDS")[!is.na(yyyymm)]
# tmp <- tmp[yyyymm %in% data[, unique(yyyymm)]]
# tmp <- tmp[, names(data), with = F]

# setkey(data, yyyymm, permno)
# setkey(tmp, yyyymm, permno)
# cor(data[, effective_spread], tmp[, effective_spread], use = "complete.obs")
# cor(data[, size], tmp[, size], use = "complete.obs")
# cor(data[, turnover], tmp[, turnover], use = "complete.obs")


# --- 2) save the list of control names

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

# # -- sanity: identical
# tt <- readRDS("../tmp/raw_data/controls/controls_classification.RDS")
# tt <- tt[control_type != 'industry']
# mean(tmp == tt)

to_file <- "../tmp/raw_data/controls/controls_classification.RDS"
dir.create(dirname(to_file), recursive = T, showWarnings = F)
saveRDS(tmp, to_file)
