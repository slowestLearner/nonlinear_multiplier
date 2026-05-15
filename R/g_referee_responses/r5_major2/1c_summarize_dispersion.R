# --- Let's check if larger OFI shocks are associated with larger dispersion
library(this.path)
setwd(this.path::this.dir())
source("../../utilities/runmefirst.R")
options(width = 200)

data <- readRDS("../../tmp/raw_data/reg_inputs/reg_table_static.RDS")[type == "OFI"]
data <- data[, .(yyyymm, permno, ofi, bin)]

# merge with dispersion measures
tmp <- readRDS("tmp/raw_files/ofi_dispersion.RDS")[, c("obs", "ofi_sum", "abs_ofi_sum") := NULL]
tmp <- melt(tmp, id.vars = c("yyyymm", "permno"), variable.name = "disp_type", value.name = "disp") %>%
    mutate(disp_type = as.character(disp_type)) %>%
    setDT()
data <- merge(data, tmp, by = c("yyyymm", "permno"), allow.cartesian = TRUE)
rm(tmp)

# look over a few periods of time
tmp <- unique(data[, .(yyyymm)])[, yyyy := floor(yyyymm / 100)]
data <- merge(data, tmp, by = "yyyymm")

out <- data[, .(disp_mean = mean(disp), disp_q25 = quantile(disp, .25), disp_q50 = quantile(disp, .50), disp_q75 = quantile(disp, .75)), .(yyyy, bin, disp_type)]
saveRDS(out, "tmp/raw_files/ofi_dispersion_by_year_and_bin.RDS")
