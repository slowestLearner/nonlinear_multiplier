# ------ plot figures: show what happens to FIT when we do PCA on flow
library(this.path)
setwd(this.path::this.dir())
source("../utilities/runmefirst.R")
library(latex2exp)
library(styler)

# --- create directories
# to_dir <- "../output/figs/data/cleaning/fit/"
# dir.create(to_dir, showWarnings = F, recursive = T)

# data
data <- readRDS("../tmp/additional/fit_alternative_construction.RDS")[origin == "flow" & !(type %in% c("flow_origin", "flow_resid"))]

# parse lables
tmp <- unique(data[, .(type)])
tmp[, num_pcs := as.integer(gsub("flow_took_out_|_pcs", "", type))]
tmp[, type_lab := paste0("Removed ", num_pcs, " PCs")]
data <- merge(data, tmp, by = "type")
rm(tmp)

this_var <- "FIT"
pp <- ggplot(data[var == this_var], aes(x = flow_raw, y = flow_mean, color = reorder(type_lab, num_pcs))) +
    geom_line(lwd = 1) +
    geom_point(cex = 1) +
    theme_classic() +
    theme(legend.position = c(.2, .8)) +
    labs(x = "Original FIT", y = "PC-removed FIT") +
    scale_x_continuous(labels = scales::percent_format(0.1)) +
    scale_y_continuous(labels = scales::percent_format(0.1)) +
    theme(text = element_text(size = 25), legend.title = element_blank(), legend.position = c(.2, .8)) +
    #   geom_hline(yintercept = 0, lty = 2) + geom_vline(xintercept = 0, lty = 2) +
    geom_abline(slope = 1, intercept = 0, lty = 2)

to_dir <- "../output/figs/data/cleaning/fit_pca/"
dir.create(to_dir, showWarnings = F, recursive = T)
ggsave(paste0(to_dir, "fit_removing_pcs.png"), pp, "png", w = 5, h = 4.5)
