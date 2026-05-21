# --- Compute rolling market volatility
library(zoo)
library(this.path)
setwd(this.path::this.dir())
source("../../utilities/runmefirst.R")

# daily market factor returns
data <- readRDS("../../../../../data/factors/daily_ff.RDS")[, .(date, ret = mktrf)] %>%
  arrange(date) %>%
  setDT()

# compute rolling vol, annualized
data[, rolling_vol_1y := frollapply(ret, 252 * 1, sd) * sqrt(252)]
data[, rolling_vol_3y := frollapply(ret, 252 * 3, sd) * sqrt(252)]
data[, rolling_vol_5y := frollapply(ret, 252 * 5, sd) * sqrt(252)]
data[, ret := NULL]

to_file <- "tmp/vols/market_vol.RDS"
dir.create(dirname(to_file), showWarnings = FALSE, recursive = TRUE)
saveRDS(data, to_file)

# # --- sanity: compare with vix. okay, so vix is much more fast moving, and is forward looking

# # realized vol
# data <- readRDS("tmp/vols/market_vol.RDS") %>%
#   melt(., id.vars = "date", variable.name = "type", value.name = "vol") %>%
#   mutate(type = as.character(type)) %>%
#   na.omit() %>%
#   setDT()

# # vix
# tmp <- readRDS("../../../../../data/factors/daily_vix.RDS")

# common_dates <- intersect(data[, date], tmp[, date])
# data <- data[date %in% common_dates]
# tmp <- tmp[date %in% common_dates]
# tmp <- tmp[, .(date, type = 'vix', vol = vix/100)]
# data <- rbind(data, tmp)
# rm(tmp)

# ggplot(data, aes(x = date, y = vol, color = type)) +
#   geom_line() +
#   theme_classic() +
#   theme(text = element_text(size = 30), legend.position = "bottom", legend.title = element_blank()) +
#   labs(x = element_blank(), y = element_blank()) +
#   ggtitle("Market Volatility") +
#   geom_hline(yintercept = 0, lty = 2)
