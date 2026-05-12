# --- What have we already controlled for?
library(this.path)
setwd(this.path::this.dir())
source("../../utilities/runmefirst.R")

data <- readRDS("../../tmp/raw_data/controls/controls_classification.RDS")

vv <- data[control_type == "return-predictor", var_lab]
paste0(vv, collapse = ", ")

vv <- data[control_type == "liquidity", var_lab]
paste0(vv, collapse = ", ")
