rm(list = ls())
gc()

# load libraries
library(dplyr)
library(tidyr)
library(fixest)
library(ggplot2)
library(parallel)
library(DescTools)
library(tictoc)
library(haven)
library(rstudioapi)
options(scipen = 10)
library(showtext)
library(sandwich)
library(styler)
showtext_auto()
library(data.table)

# num of cores
nc <- detectCores() - 2

# library(arrow)
# library(reshape2)
