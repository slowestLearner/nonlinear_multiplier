# --- Let's plot
library(this.path)
setwd(this.path::this.dir())
source("../../utilities/runmefirst.R")
source("../../utilities/regressions.R")
library(latex2exp)
options(width = 120)

out_all <- readRDS("tmp/reg_results/reg_results_fm_12ind_separate.RDS")
out_all[, var_type := ifelse(var %in% paste0("ofi_bin", 1:3), "coef", "diff")]
out_all[, se := se_nw][, se_nw := NULL] # use the right standard errors

# mark jkp
out_all[, var_added := gsub("jkp_", "JKP-part ", var_added)]

# start from adding stuff
out_all <- out_all[spec_idx >= 3]
out_all[, spec_idx := spec_idx - 2]
out_all[spec_idx == 1, var_added := ""]

# mark variable names
tmp <- unique(out_all[, .(var_type, var)])[order(var_type, var)]
tmp[, var_idx := .I]
out_all <- merge(out_all, tmp, by = c("var", "var_type"))

full_label_map <- TeX(c(
  "1" = r"($M_{|d| < \sigma}$)",
  "2" = r"($M_{|d| \in [\sigma, 2\sigma]}$)",
  "3" = r"($M_{|d| > 2\sigma}$)",
  "4" = r"($M_{|d| \in [\sigma, 2\sigma]} - M_{|d| < \sigma}$)",
  "5" = r"($M_{|d| > 2\sigma} - M_{|d| < \sigma}$)",
  "6" = r"($M_{|d| > 2\sigma} - M_{|d| \in [\sigma, 2\sigma]}$)"
))

rm(tmp)

# cut_location <- 18
cut_location <- 12.5

to_dir <- "../../output/figs/referee_responses/r1_minor3/"
dir.create(to_dir, recursive = T, showWarnings = F)

for (this_type in unique(out_all[, type])) {
  # this_type <- 'FIT'
  out <- out_all[var_type == "coef" & type == this_type][order(var, spec_idx)]
  pp <- ggplot(out, aes(x = spec_idx, y = coef, fill = factor(var_idx))) +
    geom_line(aes(color = factor(var_idx)), lwd = 1) +
    geom_point(aes(color = factor(var_idx)), cex = 1) +
    geom_ribbon(aes(ymin = coef - se, ymax = coef + se), alpha = 0.2) +
    scale_x_discrete(limits = factor(out[var == first(var), spec_idx]), labels = out[var == first(var), var_added]) +
    theme_classic() +
    scale_color_discrete(labels = full_label_map) +
    scale_fill_discrete(labels = full_label_map) +
    labs(x = "Control added", y = "Coefficient") +
    geom_hline(yintercept = 0, lty = 3, lwd = .75) +
    # geom_vline(xintercept = cut_location, lty = 3, lwd = .75) +
    coord_cartesian(ylim = c(NA, 5.5))

  if (this_type == "FIT") {
    pp <- pp + theme(text = element_text(size = 30), legend.position = c(.2, .22), legend.title = element_blank(), axis.text.x = element_text(angle = 60, vjust = .6))
  } else {
    pp <- pp + theme(text = element_text(size = 30), legend.position = c(.8, 1.8), legend.title = element_blank(), axis.text.x = element_text(angle = 60, vjust = .6))
  }

  ggsave(paste0(to_dir, this_type, "_coef.png"), pp, "png", w = 4, h = 3.5, dpi = 300, units = "in")

  out <- out_all[var_type == "diff" & type == this_type][order(var, spec_idx)]
  pp <- ggplot(out, aes(x = spec_idx, y = coef, fill = factor(var_idx))) +
    geom_line(aes(color = factor(var_idx)), lwd = 1) +
    geom_point(aes(color = factor(var_idx)), cex = 1) +
    geom_ribbon(aes(ymin = coef - se, ymax = coef + se), alpha = 0.2) +
    scale_x_discrete(limits = factor(out[var == first(var), spec_idx]), labels = out[var == first(var), var_added]) +
    theme_classic() +
    scale_color_discrete(labels = full_label_map) +
    scale_fill_discrete(labels = full_label_map) +
    labs(x = "Control added", y = "Coefficient") +
    geom_hline(yintercept = 0, lty = 3, lwd = .75) +
    # geom_vline(xintercept = cut_location, lty = 3, lwd = .75) +
    coord_cartesian(ylim = c(-2.3, NA))

  if (this_type == "FIT") {
    pp <- pp + theme(text = element_text(size = 30), legend.position = c(.82, .2), legend.title = element_blank(), axis.text.x = element_text(angle = 60, vjust = .6))
  } else {
    pp <- pp + theme(text = element_text(size = 30), legend.position = c(.82, 1.2), legend.title = element_blank(), axis.text.x = element_text(angle = 60, vjust = .6))
  }

  ggsave(paste0(to_dir, this_type, "_diff.png"), pp, "png", w = 4, h = 3.5, dpi = 300, units = "in")
}
