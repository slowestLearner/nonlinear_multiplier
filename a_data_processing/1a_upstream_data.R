# === get the unexpected components of OFI in daily data
library(this.path)
setwd(this.path::this.dir())
source("../utilities/runmefirst.R")

# -- russell rdd data
data <- data.table(haven::read_stata("../../../tests/9_russell_rdd/russell_rdd_stock_level_panel.dta"))
to_dir <- "../../../data/demand_shocks/bmi/"
dir.create(to_dir, recursive = T, showWarnings = F)
saveRDS(data, paste0(to_dir, "russell_rdd_stock_level_panel.RDS"))

# -- changes in IO
data <- readRDS("../../../tests/26_bmi_pass_through/tmp/raw_data/io_changes.RDS")[, .(yyyymm, permno, dio)]
data[, dio := Winsorize(dio, quantile(dio, probs = c(0.01, 0.99)))]
saveRDS(data, "../../../data/institutional/s34_io_changes.RDS")

# --- flow residuals
data <- readRDS("../../../tests/27_alternative_fit_construtions/tmp/raw_data/flow_residuals/2_pca_include_raw_origin.RDS")
saveRDS(data, "../../../data/demand_shocks/j_fit/quarterly_flow_residuals.RDS")

# -- FIT residuals (later port over the upstream processing code)
data <- readRDS("../../../tests/27_alternative_fit_construtions/tmp/raw_data/flow_residuals/3_fit_include_raw_origin.RDS")
saveRDS(data, "../../../data/demand_shocks/j_fit/quarterly_residuals.RDS")
