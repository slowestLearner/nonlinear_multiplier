# --- compute FIT based on PCA residual fund flows. Need to use holdings data
library(this.path)
setwd(this.path::this.dir())
source("../utilities/runmefirst.R")

# cleaned fund flows
fdata <- readRDS("../tmp/additional/clean_fit/flow_residuals/1_pca_flow_residuals.RDS")
fdata <- melt(fdata, id.vars = c("yyyymm", "wficn"), variable.name = "flow_type", value.name = "flow") %>%
  mutate(flow_type = as.character(flow_type)) %>%
  setDT()

# Produce a few versions with different degrees of winsorization
fdata[, flow_cut01 := Winsorize(flow, val = quantile(flow, probs = c(.005, .995))), .(flow_type)]
fdata[, flow_cut05 := Winsorize(flow, val = quantile(flow, probs = c(.025, .975))), .(flow_type)]
fdata[, flow_cut10 := Winsorize(flow, val = quantile(flow, probs = c(.05, .95))), .(flow_type)]

# get lagged shares outstanding
sdata <- readRDS("../../../../data/stocks/prices/quarterly_return.RDS")[, .(yyyymm, permno, shrout_1)] %>% na.omit()

# helper function to compute FIT for one period
p.get_one <- function(this_file) {
  # this_file <- files[1]

  # read in holdings
  data <- readRDS(paste0(from_dir, this_file))

  # turn holdings to quarterly
  tt <- unique(data[, .(rdate, wficn)])
  tt[, yyyymm := 100 * year(rdate) + month(rdate)]
  tt[, mm := yyyymm - floor(yyyymm / 100) * 100]
  tt[, qq := ceiling(mm / 3) * 3]
  tt[, yyyyqq := yyyymm + qq - mm]
  tt[, c("mm", "qq") := NULL]
  tt <- tt[order(yyyymm)]
  tt <- tt[, .(rdate = last(rdate)), .(wficn, yyyyqq)] %>% rename(yyyymm = yyyyqq)
  data <- merge(data, tt, by = c("rdate", "wficn"))
  rm(tt)

  data <- data[, .(yyyymm, wficn, permno, shares, prc = prc_tfn, shrout = shrout_tfn)] %>%
    na.omit() %>%
    setDT()

  data <- merge(data, data[, .(obs = length(permno), aum = sum(shares * prc)), .(yyyymm, wficn)], by = c("yyyymm", "wficn"))
  data <- data[obs >= 20]
  data[, w := shares * prc / aum]

  # lag
  data[, mm := yyyymm - floor(yyyymm / 100) * 100]
  data <- data[mm %in% c(3, 6, 9, 12)]
  data <- data[, list(yyyymm = ifelse(mm == 12, yyyymm + 100 - 9, yyyymm + 3), wficn, permno, frac_held_1 = shares / shrout, w_1 = w)]
  data <- data[0 == rowSums(is.na(data))]

  # sanity check
  data[frac_held_1 > .1, frac_held_1 := .1]
  data <- merge(data, sdata, by = c("yyyymm", "permno"))

  # get flows
  data <- merge(data, fdata, by = c("yyyymm", "wficn"), allow.cartesian = TRUE)

  # compute FIT and output
  out <- data[, list(
    obs_wficn = length(wficn),
    fit2shrout = sum(flow * frac_held_1),
    fit2shrout_cut01 = sum(flow * frac_held_1),
    fit2shrout_cut05 = sum(flow * frac_held_1),
    fit2shrout_cut10 = sum(flow * frac_held_1),
    frac_held_1 = sum(frac_held_1),
    shrout_1 = last(shrout_1)
  ), list(yyyymm, permno, flow_type)]

  return(out)
}

# loop through holdings by year
from_dir <- "../../../../data/institutional/s12_holdings/"
files <- list.files(from_dir)
tmp <- data.table(ff = files, idx = 1:length(files))
block_size <- nc
tmp[, block_idx := ceiling(idx / block_size)]

# takes 1-2 mins on a 6 core computer
plan(multisession, workers = detectCores() - 2)

data <- data.table()
for (i in unique(tmp[, block_idx])) {
  tic(i / max(tmp[, block_idx]))
  data <- rbind(data, rbindlist(future_lapply(tmp[block_idx == i, ff], p.get_one, future.seed = 123, future.packages = "data.table")))
  gc()
  toc()
}
plan(sequential)
rm(from_dir, i, files, tmp, sdata, fdata)

to_file <- "../tmp/additional/clean_fit/flow_residuals/2_fit.RDS"
dir.create(dirname(to_file), showWarnings = FALSE, recursive = TRUE)
saveRDS(data, to_file)

# # --- SANITY CHECK: not identical to earlier, due to the various differences

# data <- readRDS("../tmp/additional/clean_fit/flow_residuals/2_fit.RDS")
# data_old <- readRDS("../../../../tests/27_alternative_fit_construtions/tmp/raw_data/flow_residuals/3_fit_include_raw_origin.RDS")
# data_old <- data_old[origin == 'flow'][, origin := NULL]
# data <- merge(data, data_old, by = c("yyyymm", "permno", "flow_type"), all = T)
# data[, cor(fit2shrout.x, fit2shrout.y, use = "complete.obs")]
