# --- Modified from 2a in main code
library(this.path)
setwd(this.path::this.dir())
source("../../utilities/runmefirst.R")
source("../../utilities/regressions.R")
options(width = 200)
nc <- detectCores() - 2

n_lags <- 4

# load regression data
tic("preparing data")
# data
data <- readRDS("../../tmp/raw_data/reg_inputs/reg_table_static.RDS")[type %in% c("FIT", "OFI")]

# get lagged demand, turn into unif[-0.5, 0.5]
tmp <- data[, .(yyyymm, permno, type, absofi = abs(ofi))]
tmp[, absofi := frank(absofi, ties.method = "dense") / .N - .5, .(yyyymm, type)] # nnormalize

# some lags
tmp[, idx := frank(yyyymm, ties.method = "dense")]
for (i in 1:n_lags) {
  tmp <- merge(tmp, tmp[, .(idx = idx + i, type, permno, xx = absofi)], by = c("idx", "type", "permno"), all.x = T)
  setnames(tmp, "xx", paste0("absofi_", i))
}
tmp[, c("idx", "absofi") := NULL]
tmp <- tmp %>% na.omit()

# put together
data[, c("bin", "ofi_absofi", "ofi_bin1", "ofi_bin2", "ofi_bin3") := NULL]
data <- merge(data, tmp, by = c("yyyymm", "permno", "type"))
rm(tmp)

# create interaced variables
data[, ofi_absofi_1 := ofi * absofi_1]
data[, ofi_absofi_2 := ofi * absofi_2]
data[, ofi_absofi_3 := ofi * absofi_3]
data[, ofi_absofi_4 := ofi * absofi_4]
data[, c("absofi_1", "absofi_2", "absofi_3", "absofi_4") := NULL]
data_all <- copy(data)
rm(data)
gc()

# get control variable names
cdata <- readRDS("../../tmp/raw_data/controls/controls_classification.RDS")
controls_char <- cdata[control_type == "return-predictor", var]
controls_liq <- cdata[control_type == "liquidity", var]
controls_list <- c(controls_char, controls_liq)
rm(cdata)

toc()


# ###########################################################
# Fama-MacBeth
# ###########################################################


# worker function to estimate regression with one type of demand variable. reg_spec = "nonlinear" or "stdev"
# data <- data_list[[1]]
p.process_one_type <- function(data) {
  # parse
  this_type <- data[1, type]

  # control variables for different regression specifications
  controls <- paste0(c(controls_char, controls_liq), collapse = "+")

  specs <- list(
    "ret ~ ofi + ofi_absofi_1 + ofi_absofi_2 + ofi_absofi_3 + ofi_absofi_4",
    "ret ~ ofi + ofi_absofi_1",
    "ret ~ ofi + ofi_absofi_2",
    "ret ~ ofi + ofi_absofi_3",
    "ret ~ ofi + ofi_absofi_4"
  )

  out_all <- data.table() # save results here

  # introduce direct controls
  for (spec_idx in 1:length(specs)) {
    # print(spec_idx)
    ff <- paste0(specs[spec_idx], " + ", controls)
    out <- p.fama_macbeth(data, ff, compare_coefs = ifelse(spec_idx == 1, TRUE, FALSE))
    out[, spec_idx := spec_idx]
    out_all <- rbind(out_all, out)
  }

  # return
  out_all[, type := this_type]
  return(out_all)
}

tic("regression")
data_list <- split(data_all, by = "type")
out <- rbindlist(mclapply(data_list, p.process_one_type, mc.cores = nc))
toc()

to_file <- "tmp/reg_results/fm.RDS"
dir.create(dirname(to_file), recursive = T, showWarnings = F)
saveRDS(out, to_file)

# # --- take a look at results

# out_all <- readRDS("tmp/reg_results/fm.RDS")[var %in% c("ofi", paste0("ofi_absofi_", 1:4))]
# out_all[, coef_round := round(coef, 2)]
# out_all[, tstat_round := round(coef / se, 2)]

# dcast(out_all, var ~ type + spec_idx, value.var = "coef_round")
# dcast(out_all, var ~ type + spec_idx, value.var = "tstat_round")
