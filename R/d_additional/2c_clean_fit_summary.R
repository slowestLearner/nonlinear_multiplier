# --- Get basic summaries of alternative FIT constructions, for plotting
library(this.path)
setwd(this.path::this.dir())
source("../utilities/runmefirst.R")

# fund flows
data <- readRDS("../tmp/additional/clean_fit/flow_residuals/1_pca_flow_residuals.RDS")
data <- melt(data, id.vars = c("yyyymm", "wficn"), value.var = "flow", variable.name = "type", value.name = "flow") %>%
    mutate(type = as.character(type)) %>%
    rename(id = wficn) %>%
    mutate(var = "flow") %>%
    setDT()

# FIT
tmp <- readRDS("../tmp/additional/clean_fit/flow_residuals/2_fit.RDS")[, .(yyyymm, id = permno, var = "FIT", type = flow_type, flow = fit2shrout_cut01)]
data <- rbind(data, tmp)
rm(tmp)

# center for plotting, slightly winsorize
data[, flow := flow - mean(flow), .(yyyymm, var, type, origin)]
data[, flow := Winsorize(flow, quantile(flow, probs = c(.005, .995))), .(yyyymm, var, type, origin)]
data <- data %>%
    na.omit() %>%
    setDT()

# merge with original raw flows/FIT (before lceaning)
data_bk <- copy(data)
data <- merge(data_bk, data_bk[type == "flow", .(yyyymm, id, var, flow_raw = flow)], by = c("yyyymm", "id", "var"), allow.cartesian = T)
data[, bin := ntile(flow_raw, 20), .(yyyymm, var, type)]

out <- data[, .(
    flow_raw = mean(flow_raw),
    flow_mean = mean(flow),
    flow_q25 = quantile(flow, probs = .25),
    flow_q50 = quantile(flow, probs = .50),
    flow_q75 = quantile(flow, probs = .75)
), .(var, type, bin)]


# save
to_file <- "../tmp/additional/clean_fit/summary.RDS"
dir.create(dirname(to_file), showWarnings = FALSE, recursive = TRUE)
saveRDS(out, to_file)
