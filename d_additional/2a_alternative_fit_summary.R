# --- basic facts about alternative construction
library(this.path)
setwd(this.path::this.dir())
source("../utilities/runmefirst.R")

# # get a sense of explanatory power
# data <- readRDS("../../../data/demand_shocks/j_fit/quarterly_flow_residuals.RDS")
# data <- melt(data, id.vars = c("yyyymm", "wficn", "origin"), value.var = "flow", variable.name = "type", value.name = "flow") %>%
#     mutate(type = as.character(type)) %>%
#     setDT()

# # get dispersions
# out <- data[, .(flow_var = var(flow)), .(yyyymm, type, origin)][, .(flow_var = mean(flow_var)), .(type, origin)]
# out[, r2 := 1 - flow_var / out[1, flow_var]]
# out <- out[order(origin)]

# --- get cross-sectional summaries of flows

# flows
data <- readRDS("../../../data/demand_shocks/j_fit/quarterly_flow_residuals.RDS")
data <- melt(data, id.vars = c("yyyymm", "wficn", "origin"), value.var = "flow", variable.name = "type", value.name = "flow") %>%
    mutate(type = as.character(type)) %>%
    rename(id = wficn) %>%
    mutate(var = "flow") %>%
    setDT()

# FIT
tmp <- readRDS("../../../data/demand_shocks/j_fit/quarterly_residuals.RDS")[, .(yyyymm, id = permno, origin, var = "FIT", type = flow_type, flow = fit2shrout_cut01)]
tt <- readRDS("../../../data/demand_shocks/j_fit/quarterly_updated_20250313.RDS")[, .(yyyymm, id = permno, var = "FIT", type = "flow_origin", flow = fit_adj_cut01)]
tt <- merge(tt, unique(tmp[, .(yyyymm, id)]), by = c("yyyymm", "id"))
tt[, origin := "flow"]
tmp <- rbind(tt, tmp)
tmp[, type := as.character(type)]
data <- rbind(data, tmp)
rm(tmp, tt)

# take out zeros by period (because the flow residuals are zero on average). slightly winsorize
data[, flow := flow - mean(flow), .(yyyymm, var, type, origin)]
data[, flow := Winsorize(flow, quantile(flow, probs = c(.005, .995))), .(yyyymm, var, type, origin)]
data <- data %>%
    na.omit() %>%
    setDT()
# data_bk <- copy(data)

# merge with really raw flows
data <- copy(data_bk)
data <- merge(data, data_bk[origin == "flow" & type == "flow_origin", .(yyyymm, id, var, flow_raw = flow)], by = c("yyyymm", "id", "var"), allow.cartesian = T)
data[, bin := ntile(flow_raw, 20), .(yyyymm, var, type, origin)]

out <- data[, .(
    flow_raw = mean(flow_raw),
    flow_mean = mean(flow),
    flow_q25 = quantile(flow, probs = .25),
    flow_q50 = quantile(flow, probs = .50),
    flow_q75 = quantile(flow, probs = .75)
), .(var, type, bin, origin)]

# hmm, there is quite a bit of change here, i think
saveRDS(out, "../tmp/additional/fit_alternative_construction.RDS")

# --- plot to take a look
out_all <- readRDS("../tmp/additional/fit_alternative_construction.RDS")

this_var <- "FIT"
# this_origin <- "flow"
this_origin <- "flow_resid"
out <- out_all[var == this_var & origin == this_origin]

tt <- unique(out[type != "flow_origin", .(type)])[, type_lab := paste0(.I, "_", type)]
out <- merge(out, tt, by = "type")

xx <- range(out[, flow_raw])
ggplot(
    out,
    aes(x = flow_raw, y = flow_mean, color = type_lab)
) +
    geom_line(lwd = 5) +
    geom_point(cex = 5) +
    coord_cartesian(xlim = xx, ylim = xx) +
    geom_abline(slope = 1, intercept = 0, lty = 2) +
    theme_classic() +
    theme(text = element_text(size = 25), legend.position = c(.2, .8)) +
    ggtitle(paste0("var = ", this_var, ", origin = ", this_origin))
