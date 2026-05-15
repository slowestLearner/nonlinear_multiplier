library(this.path)
setwd(this.path::this.dir())
source("../../utilities/runmefirst.R")
source("../../utilities/regressions.R")
options(width = 200)

data <- readRDS("../../tmp/raw_data/reg_inputs/reg_table_static.RDS")[type == "BMI"]
table(data[, yyyymm])
