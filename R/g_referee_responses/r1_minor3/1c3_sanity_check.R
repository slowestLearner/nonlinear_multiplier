# --- Let's see if the PCA residuals look somewhat kosher
library(this.path)
setwd(this.path::this.dir())
source("../../utilities/runmefirst.R")
library(RSpectra) # Faster for extracting top K components
options(width = 200)

# --- quarterly returns
data <- readRDS("../../tmp/raw_data/reg_inputs/reg_table_static.RDS")[type %in% c("FIT", "OFI"), .(yyyymm, permno, ret)] %>% unique()
tmp <- readRDS("tmp/pca_residuals/quarterly_oos.RDS")
data <- merge(data, tmp, by = c("yyyymm", "permno"), all.x = T)
rm(tmp)

round(cor(data[, -c(1:2), with = F]), 3)
round(cov(data[, -c(1:2), with = F]), 4)


# --- monthly returns
data <- readRDS("../../tmp/raw_data/reg_inputs/reg_table_static.RDS")[type == "BMI", .(yyyymm, permno, ret)] %>% unique()
tmp <- readRDS("tmp/pca_residuals/monthly_oos.RDS")
data <- merge(data, tmp, by = c("yyyymm", "permno"), all.x = T)
rm(tmp)

round(cor(data[, -c(1:2), with = F]), 3)
round(cov(data[, -c(1:2), with = F]), 4)
