# ----- put together liquidity control variables for regressions
# code adapted from /Users/slowlearner/Dropbox/SpeculativeIdeas/Nonlinear_multipliers/tests/13_main_result_robustness/6_liquidity_interaction.Rmd
library(this.path)
setwd(this.path::this.dir())
source("../utilities/runmefirst.R")

# --- 1) get liquidity variables

# read spreads from TAQ summaries
data <- readRDS("../../../../data/stocks/liquidity/monthly_file.RDS")
tmp <- readRDS("../../../../data/stocks/liquidity/daily_file.RDS")
data <- rbind(data, tmp)
rm(tmp)

# slightly winsorize, turn to monthly
data[, effective_spread := Winsorize(effective_spread, val = quantile(effective_spread, probs = c(0, .995), na.rm = T))]
data[, quoted_spread := Winsorize(quoted_spread, val = quantile(quoted_spread, probs = c(0, .995), na.rm = T))]

data[, yyyymm := 100 * year(date) + month(date)]
data <- data[, .(
    effective_spread = mean(effective_spread, na.rm = T),
    quoted_spread = mean(quoted_spread, na.rm = T)
), .(yyyymm, permno)]
gc()

# merge with volatilities
tmp <- readRDS("../../../../data/stocks/controls/monthly_characteristics_not_lagged.RDS")[, .(yyyymm, permno, realized_vol)]
data <- merge(data, tmp, by = c("yyyymm", "permno"), all.x = T)
rm(tmp)

# also get turnover and dollar volume
tmp <- readRDS("../../../../data/stocks/prices/monthly_return.RDS")
tmp <- tmp[, .(yyyymm, permno, me = prc * shrout / 1e3, turnover = (vol / 10) / shrout)] %>%
    filter(me > 0) %>%
    na.omit() %>%
    setDT()
tmp[, turnover := Winsorize(turnover, val = quantile(turnover, probs = c(0, .995)))]
tmp[, dollar_volume := turnover * me]
data <- merge(data, tmp, by = c("yyyymm", "permno"), all.x = T)
# data_bk <- copy(data)
rm(tmp)


# --- 2) turn into uniform distributions

# flip directions, so larger values means more liquid (does not impact regression results)
data[, effective_spread := -effective_spread]
data[, quoted_spread := -quoted_spread]
data <- data.table(melt(data, id.vars = c("yyyymm", "permno"))) %>%
    na.omit() %>%
    setDT()
names(data)[3:4] <- c("var", "char")

# rank by (yyyymm, var)
data <- data[order(yyyymm, var, char)]
data[, char_unif := (seq_len(.N) - 1) / (.N - 1) - 0.5, by = .(yyyymm, var)]
data

# make wide and save
data[, char := NULL]
data <- dcast(data, yyyymm + permno ~ var, value.var = "char_unif")
data <- data[yyyymm <= 202212]
to_file <- "../../../../data/stocks/controls/monthly_liquidity_measures_not_lagged.RDS"
dir.create(dirname(to_file), recursive = T, showWarnings = F)
saveRDS(data, to_file)

# # # --- sanity check with earlier. Sufficiently similar

# data <- readRDS("../../../../data/stocks/liquidity/liquidity_measures_put_together.RDS")
# tmp <- readRDS("../../../../data/stocks/controls/monthly_liquidity_measures_not_lagged.RDS")
# tmp <- tmp[yyyymm <= 202212]

# dim(data)
# dim(tmp)

# setkey(data, yyyymm, permno)
# setkey(tmp, yyyymm, permno)

# cor(data[, turnover], tmp[, turnover], use = "complete.obs")
# cor(data[, dollar_volume], tmp[, dollar_volume], use = "complete.obs")
# cor(data[, effective_spread], tmp[, effective_spread], use = "complete.obs")
