# --- compute multipliers using the delta method
library(this.path)
setwd(this.path::this.dir())
source("~/.runmefirst")
options(width = 200)

# get regression results
data <- readRDS("tmp/regression_contemp/4_bins/fm_stdev.RDS")[spec_idx == 3]
data <- data[var %in% paste0("ofi_bin", 1:4)]
data <- data[, .(bin = as.integer(sub("ofi_bin", "", var)), type, M = coef, se)]

# get distance in |d|
tmp <- readRDS("tmp/regression_contemp/4_bins/abs_ofi_magnitude.RDS") %>% dplyr::rename(d = abs_ofi)
data <- merge(data, tmp, by = c("type", "bin"))[order(type, bin)]
rm(tmp)

data[, cum_price_impact := M * d]
data[, se_cum_price_impact := se * d]

# let's plot?
ggplot(data, aes(x = d, y = cum_price_impact, fill = type)) +
    geom_point(cex = 5, aes(color = type)) +
    geom_line(lwd = 3, aes(color = type)) +
    geom_ribbon(aes(ymin = cum_price_impact - 1.96 * se_cum_price_impact, ymax = cum_price_impact + 1.96 * se_cum_price_impact), alpha = 0.2) +
    labs(x = "|d|", y = "Cumulative price impact") +
    theme_classic() +
    theme(text = element_text(size = 35), legend.position = c(.8, .2), legend.title = element_blank()) +
    geom_hline(yintercept = 0, lty = 3)
