# ------ plot figures: show that progressively adding controls do not impact inference
library(this.path)
setwd(this.path::this.dir())
source("../utilities/runmefirst.R")
library(latex2exp)
library(styler)

# --- create directories
tic("loading data")
to_dir <- "../output/figs/price_impact/static/more_controls/"
dir.create(to_dir, recursive = T, showWarnings = F)

# regression results
data <- readRDS("../tmp/price_impact/regression_contemp/fm_stdev.RDS")[type != "OFI_pre_whitened"]

# start from the spedification before adding demand-based interactions
data <- data[spec_idx >= 3]
data[, spec_idx := spec_idx - 2]

# choose variables
data <- data[grepl("ofi_bin", var)]
data <- data[grepl(" - ", var)]
data[, var := gsub("ofi_bin", "M", var)]

# parse variable names
data <- data[, .(type, spec_idx, coef, se, var, var_added)] # var_added denote the newest control introduced (interacted with demand)

tmp <- readRDS("../tmp/raw_data/controls/controls_classification.RDS")
tmp <- rbind(tmp, data.table(var = "none", control_type = "", var_lab = "None"))
setnames(tmp, "var", "var_added")
data <- merge(data, tmp, by = "var_added", all.x = T)
rm(tmp)

# parse variable names
data[var == "M2 - M1", var_tex := TeX("$M_{|d| \\in [\\sigma, 2\\sigma]} - M_{|d| < \\sigma}$$")]
data[var == "M3 - M1", var_tex := TeX("$M_{|d| > 2\\sigma} - M_{|d| < \\sigma}$")]
data[var == "M3 - M2", var_tex := TeX("$M_{|d| > 2\\sigma} - M_{|d| \\in [\\sigma, 2\\sigma]}$")]

data[spec_idx == 1, var_lab := "(None)"] # no controls added
data_all <- copy(data)
rm(data)
toc()

# Plot
for (this_type in unique(data_all[, type])) {
  # this_type <- 'FIT'
  tic(paste0("Plotting ", this_type))
  data <- copy(data_all[type == this_type])
  data <- data[order(spec_idx)]

  pp <- ggplot(data, aes(x = spec_idx, y = coef, fill = var, color = var)) +
    geom_line(lwd = 1) +
    geom_point() +
    geom_ribbon(aes(ymin = coef - se, ymax = coef + se), color = NA, alpha = .2) +
    theme_classic() +
    geom_hline(yintercept = 0, lty = 3) +
    theme(
      legend.position = c(.2, .17), legend.title = element_blank(),
      axis.text.x = element_text(angle = 45, vjust = .6),
      text = element_text(size = 12),
      legend.text = element_text(size = 10)
    ) +
    scale_x_discrete(
      limits = data[var == first(var), factor(spec_idx)],
      labels = data[var == first(var), var_lab]
    ) +
    labs(x = "Interaction controls added", y = "Regression coefficient") +
    scale_fill_discrete(
      name = NULL, # No legend title
      labels = setNames(data$var_tex, data$var)
    ) +
    scale_color_discrete(
      name = NULL, # No legend title
      labels = setNames(data$var_tex, data$var)
    ) +
    coord_cartesian(ylim = c(-3.5, NA))

  ggsave(paste0(to_dir, this_type, ".png"), pp, "png", w = 5, h = 4, dpi = 300, units = "in")
  toc()
}
