# --- Just combine them
library(zoo)
library(this.path)
setwd(this.path::this.dir())
nc <- detectCores() - 2
source("../../utilities/runmefirst.R")

# market part
data <- readRDS("tmp/vols/market_vol.RDS")[, .(date, vol_mkt = rolling_vol_3y)] %>%
    arrange(date) %>%
    mutate(yyyymm = 100 * year(date) + 3 * quarter(date))
data <- data[, .(vol_mkt = last(vol_mkt)), .(yyyymm)] %>%
    na.omit() %>%
    setDT()

# stock vol, can be EW or VW
tmp <- readRDS("tmp/vols/stock_vol_rolling_12q.RDS")[, .(yyyymm, permno, stock_vol = vol_idio)]
tt <- readRDS("../../../../../data/stocks/prices/quarterly_return.RDS")[, .(yyyymm, permno, me_1)]
tmp <- merge(tmp, tt, by = c("yyyymm", "permno")) %>% na.omit()
tmp <- tmp[, .(vol_stock_ew = mean(stock_vol), vol_stock_vw = weighted.mean(stock_vol, me_1)), .(yyyymm)]
data <- merge(data, tmp, by = "yyyymm")
rm(tt, tmp)

saveRDS(data, "tmp/vols/mkt_and_stock_combined.RDS")
