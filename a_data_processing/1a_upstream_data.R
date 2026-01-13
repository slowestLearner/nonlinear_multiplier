# === get the unexpected components of OFI in daily data
library(this.path)
setwd(this.path::this.dir())
source("../utilities/runmefirst.R")

# -- russell rdd data
data <- data.table(haven::read_stata("../../../tests/9_russell_rdd/russell_rdd_stock_level_panel.dta"))
to_dir <- "../../../data/demand_shocks/bmi/"
dir.create(to_dir, recursive = T, showWarnings = F)
saveRDS(data, paste0(to_dir, "russell_rdd_stock_level_panel.RDS"))
