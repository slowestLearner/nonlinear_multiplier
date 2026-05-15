# --- choose a few examples to plot and visualize
library(this.path)
setwd(this.path::this.dir())
source("../../utilities/runmefirst.R")
options(width = 200)

# focus on extreme bins
n_extreme_bin <- 20

# daily OFI
tic("loading data")
data <- readRDS("../../../../../data/demand_shocks/ofi/daily.RDS")

# map to time within a quarter
tmp <- unique(data[, .(date)])[order(date)][, yyyymm := 100 * year(date) + 3 * quarter(date)]
tmp[, first_date := first(date), yyyymm]
tmp[, last_date := last(date), yyyymm]
tmp[, time_within_quarter := as.numeric(date - first_date) / as.numeric(last_date - first_date)]
tmp[, c("first_date", "last_date") := NULL]
data <- merge(data, tmp, by = "date")
rm(tmp)

# choose some to plot
tmp <- readRDS("tmp/raw_files/ofi_dispersion.RDS")[, .(yyyymm, permno, ofi = ofi_sum, abs_ofi_sum, hhi = hhi_pow2, ofi_sd)]

# by type
tmp <- melt(tmp, id.vars = c("yyyymm", "permno", "ofi", "abs_ofi_sum"), variable.name = "disp_type", value.name = "disp") %>%
    mutate(disp_type = as.character(disp_type)) %>%
    filter(abs_ofi_sum < .4) %>%
    setDT()

# directions, etc
tmp[, dir := sign(ofi)]
tmp[, bin_abs_ofi := ntile(abs(ofi), 10), .(disp_type, dir)]
tmp[, bin_disp := ntile(disp, n_extreme_bin), .(disp_type, dir, bin_abs_ofi)]
tmp[, ofi := NULL]
data <- merge(data, tmp, by = c("yyyymm", "permno"), allow.cartesian = TRUE)
rm(tmp)

# map (yyyymm, permno) to idx
tt <- unique(data[, .(yyyymm, permno)])[, idx := .I]
data <- merge(data, tt, by = c("yyyymm", "permno"))
data_bk <- copy(data)
rm(tt)
toc()

tic("plotting examples")

# condition on ofi being large, and...
data <- copy(data_bk)
tmp <- copy(data)[bin_abs_ofi == 10 & bin_disp %in% c(1, n_extreme_bin), .(abs_ofi_sum = last(abs_ofi_sum), disp = last(disp)), .(idx, dir, disp_type, bin_disp)][order(disp_type, dir, bin_disp, disp)]

# fix a randomization scheme
set.seed(123)
tmp[, random_num := runif(nrow(tmp))]
tmp <- tmp[order(dir, disp_type, bin_disp, random_num)]

cut <- 3
tmp <- tmp[, plot_idx := .I][, plot_idx := frank(plot_idx, ties.method = "first"), .(disp_type, dir, bin_disp)][plot_idx <= cut]
tmp <- tmp[, .(idx, plot_idx)]
data <- merge(data, tmp, by = "idx")

# plot to visualize
for (this_disp_type in c("hhi", "ofi_sd")) {
    to_dir <- paste0("../../output/figs/referee_responses/r5_major2/visualize/", this_disp_type, "/")
    dir.create(to_dir, recursive = T, showWarnings = F)
    for (this_dir in c(1, -1)) {
        for (this_bin in c(1, n_extreme_bin)) {
            out <- data[dir == this_dir & bin_disp == this_bin & disp_type == this_disp_type]
            out[, lab := paste0("(yyyymm = ", yyyymm, ", permno = ", permno, ")")]
            out[, cumofi := cumsum(ofi), lab]

            pp <- ggplot(out, aes(x = time_within_quarter, y = cumofi, color = lab)) +
                geom_line(lwd = .75) +
                theme_classic() +
                scale_y_continuous(labels = scales::percent_format(1)) +
                theme(text = element_text(size = 12), plot.title = element_text(hjust = 0.5), legend.position = "none", legend.title = element_blank()) +
                # legend.position = if (this_dir == 1) c(.2, .8) else c(.2, .2),
                geom_hline(yintercept = 0, lty = 2) +
                labs(x = "Time within a quarter", y = "Cumulative OFI")

            ggsave(paste0(to_dir, "dir_", this_dir, "_disp_bin_", this_bin, ".png"), pp, "png", w = 4, h = 3.5)
        }
    }
}
toc()
