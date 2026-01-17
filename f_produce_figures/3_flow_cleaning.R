# ------ plot figures: the process of flow cleaning
source("utilities/runmefirst.R")
library(latex2exp)
library(styler)

# --- create directories
dir.create("output/figs/data/cleaning/fit/", showWarnings = F, recursive = T)

# --- 1) plot flow by bin
# TODO: move this to the code later
# data <- readRDS("../20250117_quarterly/tmp/raw_data/cleaning/fit/trade_and_flow_regression_table.RDS")
data <- readRDS("../20250117_quarterly/tmp/raw_data/cleaning/fit/heterogeneous_trade_to_flow.RDS")

tmp <- data[, .(flow = mean(flow)), bin_flow]

pp <- ggplot(tmp, aes(x = bin_flow, y = flow)) +
  geom_point() +
  geom_line(lwd = 1) +
  theme_classic() +
  labs(x = "Bin", y = "Fund flow") +
  scale_y_continuous(labels = scales::percent_format(1)) +
  theme_classic() +
  geom_vline(xintercept = 0, lty = 3) +
  geom_hline(yintercept = 0, lty = 3) +
  theme(text = element_text(size = 12))
dir.create("output/figs/data/cleaning/fit/", recursive = T)
ggsave("output/figs/data/cleaning/fit/flow_by_bin.png", pp, "png", w = 5, h = 4.5)

# --- 2) plot trade response to flow
rm(list = ls())

# TODO: move this to the code later
data <- readRDS("../20250117_quarterly/tmp/raw_data/cleaning/fit/heterogeneous_trade_to_flow.RDS")
data <- data[bin %in% c(2, 11, 20)]
data[, lab := NULL]

tmp <- unique(data[, .(bin)])[order(bin)]
tmp[, lab := c("Existing weight - bottom", "Existing weight - median", "Existing weight - top")]
data <- merge(data, tmp, by = "bin")
rm(tmp)

xx <- range(c(data[, flow], data[, coef]))
pp <- ggplot(data, aes(x = flow, y = coef, fill = reorder(lab, bin))) +
  geom_line(aes(color = reorder(lab, bin)), lwd = 1) +
  geom_abline(slope = 1, intercept = 0, lty = 3) +
  theme_classic() +
  theme(legend.position = c(.7, .15), legend.title = element_blank()) +
  labs(x = "Fund flow", y = "Trading") +
  coord_cartesian(xlim = xx, ylim = xx) +
  geom_vline(xintercept = 0, lty = 3) +
  geom_hline(yintercept = 0, lty = 3) +
  scale_x_continuous(labels = scales::percent_format(1)) +
  scale_y_continuous(labels = scales::percent_format(1)) +
  theme(text = element_text(size = 12))

ggsave("output/figs/data/cleaning/fit/trade_to_flow.png", pp, "png", w = 5, h = 4.5)


# --- 3) plot flow winsorization
rm(list = ls())

# TODO: move this to the code later
data <- readRDS("../20250117_quarterly/tmp/raw_data/cleaning/fit/trade_and_flow_regression_table.RDS")
data <- data[, .(flow = last(flow)), .(yyyymm, wficn)]

cuts_1pct <- quantile(data[, flow], c(.005, .995))
cuts_5pct <- quantile(data[, flow], c(.025, .975))
cuts_10pct <- quantile(data[, flow], c(.05, .95))

pp <- ggplot(data[abs(flow) < 1], aes(x = flow)) +
  geom_density(alpha = .2, lwd = 1) +
  theme_classic() +
  geom_hline(yintercept = 0, lty = 3) +
  geom_vline(xintercept = cuts_1pct, lty = 2, col = 2) +
  geom_vline(xintercept = cuts_5pct, lty = 3, col = 3) +
  # geom_vline(xintercept = cuts_10pct, lty = 4, col = 4) +
  annotate("text", x = .8, y = 2, label = "1% cutoff", col = 2, size = 4) +
  annotate("text", x = .53, y = 3, label = "5% cutoff", col = 3, size = 4) +
  # annotate("text", x = .4, y = 4, label = "10% cutoff", col = 4, size = 4) +
  labs(x = "Fund flow", y = "Density") +
  scale_x_continuous(labels = scales::percent_format(1)) +
  theme(text = element_text(size = 12))
ggsave("output/figs/data/cleaning/fit/flow_winsorization.png", pp, "png", w = 5, h = 4.5)


# --- 4) compare FIT measures
rm(list = ls())

# TODO: move this to the code later
data <- readRDS("../../data/demand_shocks/j_fit/quarterly_updated_20250313.RDS")
out <- data[, .(yyyymm, permno, fit)]
out[, bin := ntile(fit, 100)]

data[, fit := NULL]
names(data)[3:6] <- paste0("fit_adj_v", 1:4)
data[, fit_adj_v4 := NULL] # too extreme, remove
data <- data.table(melt(data, id.vars = c("yyyymm", "permno")))
names(data)[3:4] <- c("idx", "fit_adj")
data[, idx := as.character(idx)]
data[, idx := as.integer(gsub("fit_adj_v", "", idx))]
data <- merge(out, data, by = c("yyyymm", "permno"))
rm(out)
gc()

out <- data[, .(fit = mean(fit), fit_adj = mean(fit_adj)), .(bin, idx)]
tt <- data.table(idx = 1:3, lab = c(
  "no winsorization",
  "1% winsorization",
  "5% winsorization"
))
out <- merge(out, tt, by = "idx")
rm(tt)

xx <- range(c(out[, fit], out[, fit_adj]))

pp <- ggplot(out, aes(x = fit, y = fit_adj, fill = reorder(lab, idx))) +
  geom_line(aes(color = reorder(lab, idx)), lwd = 1) +
  theme_classic() +
  theme(legend.position = c(.8, .2), legend.title = element_blank()) +
  labs(x = "Raw FIT", y = "Cleaned FIT") +
  geom_abline(slope = 1, intercept = 0, lty = 3) +
  geom_vline(xintercept = 0, lty = 3) +
  geom_hline(yintercept = 0, lty = 3) +
  coord_cartesian(xlim = xx, ylim = xx) +
  scale_x_continuous(labels = scales::percent_format(1)) +
  scale_y_continuous(labels = scales::percent_format(1)) +
  theme(text = element_text(size = 12))
ggsave("output/figs/data/cleaning/fit/fit_before_vs_after.png", pp, "png", w = 5, h = 4.5)
