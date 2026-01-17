# ------ plot figures: multipliers
library(data.table)
library(ggplot2)
library(tictoc)
library(latex2exp)

dir.create("output/figs/price_impact/static/", showWarnings = F, recursive = T)
dir.create("output/figs/price_impact/dynamic/", showWarnings = F, recursive = T)

# --- static

tic("Plotting static multipliers")
data <- readRDS("tmp/price_impact/regression_contemp/fm_stdev.RDS")

# choose specifications
data <- data[spec_idx == 3]
data <- data[var %in% paste0("ofi_bin", 1:3)]
data[, var := gsub("ofi_bin", "M", var)]
data[, type_idx := ifelse(type == "BMI", 1, ifelse(type == "FIT", 2, 3))]
data <- data[order(var, type_idx)]

pd <- position_dodge(width = 0.9) # use same dodge for bar + errorbar

pp <- ggplot(data, aes(x = type_idx, y = coef, fill = var)) +
  geom_bar(stat = "identity", position = pd) +
  geom_errorbar(aes(ymin = coef - 1.96 * se, ymax = coef + 1.96 * se), alpha = .5, position = pd, width = .5) +
  theme_classic() +
  labs(x = "Demand measure", y = "Price multiplier") +
  scale_fill_discrete(
    name = NULL,
    labels = c(
      "M1"  = TeX("$M_{|d| < \\sigma}$"),
      "M2"  = TeX("$M_{|d| \\in [\\sigma, 2\\sigma]}$"),
      "M3"  = TeX("$M_{|d| > 2\\sigma}$")
    )
  ) +
  theme(legend.position = c(.18, .8), legend.text = element_text(size = 12), text = element_text(size = 12)) +
  scale_x_discrete(limits = factor(1:3), labels = data[1:3, type])

ggsave("output/figs/price_impact/static/multiplier_by_shock_size_all_measures.png", pp, "png", w = 6, h = 4, dpi = 300, units = "in")
toc()

# --- dynamic. Plot separately for FIT and OFI

tic("Plotting dynamic multipliers")
data <- readRDS("tmp/price_impact/regression_dynamic/fm_stdev.RDS")

# choose specifications
data <- data[spec_idx == 3]
data <- data[var %in% paste0("ofi_bin", 1:3)]
data[, var := gsub("ofi_bin", "M", var)]

# parse
data[, lag := as.integer(substr(type, 5, 5))]
data[, type := substr(type, 1, 3)]
data_all <- copy(data)
rm(data)

for (this_type in unique(data_all[, type])) {
  data <- copy(data_all[type == this_type])

  pd <- position_dodge(width = 0.9) # use same dodge for bar + errorbar

  pp <- ggplot(data, aes(x = lag, y = coef, fill = var)) +
    geom_bar(stat = "identity", position = pd) +
    geom_errorbar(aes(ymin = coef - 1.96 * se, ymax = coef + 1.96 * se), alpha = .5, position = pd, width = .5) +
    theme_classic() +
    labs(x = "Lookback horizon (h quarters)", y = "Price multiplier") +
    scale_fill_discrete(
      name = NULL,
      labels = c(
        "M1"  = TeX("$M_{|\\sum_{l=1}^h d_{t-l}| < \\sigma}$"),
        "M2"  = TeX("$M_{|\\sum_{l=1}^h d_{t-l}| \\in [\\sigma, 2\\sigma]}$"),
        "M3"  = TeX("$M_{|\\sum_{l=1}^h d_{t-l}| > 2\\sigma}$")
      )
    ) +
    theme(legend.position = c(.85, .88), legend.text = element_text(size = 12), text = element_text(size = 12), axis.text = element_text(size = 12)) +
    scale_x_discrete(limits = factor(1:4)) +
    coord_cartesian(ylim = c(0, 5.5))


  dir.create("output/figs/dynamic/", showWarnings = F, recursive = T)
  ggsave(paste0("output/figs/price_impact/dynamic/multiplier_by_lag_and_stdev_bin_", this_type, ".png"),
    pp, "png",
    w = 5, # Double the width
    h = 4.5, # Double the height (maintains 10:9 ratio)
    dpi = 300, # Higher DPI for better resolution
    units = "in"
  )
}
toc()
