# ------ plot figures: show what happens to FIT when we remove PCs from flow
library(this.path)
setwd(this.path::this.dir())
source("../utilities/runmefirst.R")
library(latex2exp)


# data - PAUSE
data <- readRDS("../tmp/additional/clean_fit/flow_residuals/2_fit.RDS")
data <- merge(data[flow_type == "flow", .(yyyymm, permno, flow_raw = fit2shrout)],
    data[flow_type != "flow", .(yyyymm, permno, type = flow_type, flow = fit2shrout)],
    by = c("yyyymm", "permno"), allow.cartesian = T
)

# data <- readRDS("../tmp/additional/fit_alternative_construction.RDS")[origin == "flow" & !(type %in% c("flow_origin", "flow_resid"))]

# parse labels
tmp <- unique(data[, .(type)])
tmp[, num_pcs := as.integer(gsub("flow_took_out_|_pcs", "", type))]
tmp[, type_lab := paste0("Removed ", num_pcs, " PCs")]
data <- merge(data, tmp, by = "type")[, type := NULL]
rm(tmp)

# summarize by bin
# data_bk <- copy(data)
data[, bin := ntile(flow_raw, 20), .(yyyymm, type_lab)]
out <- data[, .(flow_raw = mean(flow_raw), flow_mean = mean(flow)), .(bin, type_lab, num_pcs)]

# recenter, as we focus on cross-sectional variation
out[, flow_raw := flow_raw - mean(flow_raw), .(type_lab)]
out[, flow_mean := flow_mean - mean(flow_mean), .(type_lab)]

pp <- ggplot(out, aes(x = flow_raw, y = flow_mean, color = reorder(type_lab, num_pcs))) +
    geom_line(lwd = 1) +
    geom_point(cex = 1) +
    theme_classic() +
    theme(legend.position = c(.25, .8)) +
    labs(x = "Original FIT", y = "PC-removed FIT") +
    scale_x_continuous(labels = scales::percent_format(0.1)) +
    scale_y_continuous(labels = scales::percent_format(0.1)) +
    theme(text = element_text(size = 12), legend.title = element_blank(), legend.position = c(.2, .8)) +
    geom_abline(slope = 1, intercept = 0, lty = 2)

to_dir <- "../output/figs/data/cleaning/fit_pca/"
dir.create(to_dir, showWarnings = F, recursive = T)
ggsave(paste0(to_dir, "fit_removing_pcs.png"), pp, "png", w = 5, h = 4.5)
