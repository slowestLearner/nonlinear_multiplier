# --- Estimate rolling multipliers over 3y windows
library(this.path)
setwd(this.path::this.dir())
source("../../utilities/runmefirst.R")
source("../../utilities/regressions.R")
options(width = 120)
library(latex2exp)

out_all <- readRDS("tmp/reg_results/rolling_multipliers.RDS")[spec_idx == 3 & var %in% paste0("ofi_bin", 1:3)]

# plot labels
full_label_map <- TeX(c(
  "ofi_bin1" = r"($M_{|d| < \sigma}$)",
  "ofi_bin2" = r"($M_{|d| \in [\sigma, 2\sigma]}$)",
  "ofi_bin3" = r"($M_{|d| > 2\sigma}$)"
))

rm(tmp)

to_dir <- "../../output/figs/referee_responses/r5_minor1/"
dir.create(to_dir, recursive = T, showWarnings = F)

for (this_type in unique(out_all[, type])) {
  # this_type <- 'OFI'
  out <- out_all[type == this_type]
  pp <- ggplot(out, aes(x = date, y = coef, fill = factor(var))) +
    geom_line(aes(color = factor(var)), lwd = 1) +
    geom_ribbon(aes(ymin = coef - se, ymax = coef + se), alpha = 0.2) +
    theme_classic() +
    scale_color_discrete(labels = full_label_map) +
    scale_fill_discrete(labels = full_label_map) +
    labs(x = element_blank(), y = "Multiplier") +
    geom_hline(yintercept = 0, lty = 3, lwd = .75) +
    theme(text = element_text(size = 10), legend.position = c(.8, .8), legend.title = element_blank())

  ggsave(paste0(to_dir, "multiplier_over_time_", this_type, ".png"), pp, "png", w = 4, h = 3.5, dpi = 300, units = "in")
}
