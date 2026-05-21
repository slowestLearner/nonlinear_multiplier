# --- Just plot
library(zoo)
library(this.path)
setwd(this.path::this.dir())
source("../../utilities/runmefirst.R")

data <- readRDS("tmp/vols/mkt_and_stock_combined.RDS") %>%
    melt(., id.vars = "yyyymm", variable.name = "type", value.name = "vol") %>%
    mutate(type = as.character(type)) %>%
    mutate(date = as.Date(paste0(as.character(yyyymm), "28"), "%Y%m%d")) %>%
    na.omit() %>%
    setDT()

# beautiful names
data[, type_lab := ifelse(type == "vol_mkt", "Market volatility", ifelse(type == "vol_stock_ew", "Stock volatility (EW)", "Stock volatility (VW)"))]

# 1. Set the multiplier for the secondary axis
mult <- 0.6
data_plot <- copy(data)

# 2. Update transformation logic:
# ONLY scale the EW series. Market and VW remain in original units.
data_plot[, vol_plot := ifelse(type_lab == "Stock volatility (EW)", vol * mult, vol)]

# 3. Generate the Plot
pp <- ggplot(data_plot, aes(x = date, y = vol_plot, color = type_lab)) +
    geom_line(lwd = .75) +
    theme_classic() +
    # geom_hline(yintercept = 0, lty = 2) +
    scale_y_continuous(
        name = "Market Vol and Stock Vol (VW)",
        # The secondary axis transformation (~ . / mult) reverts the
        # scaling only for the right-hand labels.
        sec.axis = sec_axis(~ . / mult, name = "Stock Vol (EW)")
    ) +
    theme(
        text = element_text(size = 10),
        legend.position = c(.75, .14),
        legend.title = element_blank(),
        axis.title.y.right = element_text(vjust = 1.5)
    ) +
    labs(x = element_blank()) +
    coord_cartesian(ylim = c(0, NA))

to_file <- "../../output/figs/referee_responses/r5_minor1/vol_over_time.png"
dir.create(dirname(to_file), showWarnings = FALSE, recursive = TRUE)
ggsave(to_file, pp, "png", w = 4.5, h = 3.5, dpi = 300, units = "in")

# # --- Optional: Regime Plotting for a specific type
# this_type <- "Stock volatility (EW)" # Use the label created above
# out <- data[type_lab == this_type]
# out[, bin := ntile(vol, 2)]

# ggplot(out, aes(x = date, y = vol)) +
#     geom_line(lwd = 2) +
#     geom_point(aes(color = as.factor(bin)), shape = 18, size = 8) +
#     theme_classic() +
#     geom_hline(yintercept = 0, lty = 2) +
#     theme(
#         text = element_text(size = 25),
#         legend.position = "none"
#     ) +
#     labs(x = element_blank(), y = "Volatility") +
#     ggtitle(this_type)
