# --- moves data from local files to this folder
library(this.path)
setwd(this.path::this.dir())
source("../utilities/runmefirst.R")

# -- russell rdd data (from Anna)
data <- data.table(haven::read_stata("../../../../tests/9_russell_rdd/russell_rdd_stock_level_panel.dta"))
to_dir <- "../../../../data/demand_shocks/bmi/"
dir.create(to_dir, recursive = T, showWarnings = F)
saveRDS(data, paste0(to_dir, "russell_rdd_stock_level_panel.RDS"))

# -- changes in S34IO. No need to provide detailed processing code
data <- readRDS("../../../../tests/26_bmi_pass_through/tmp/raw_data/io_changes.RDS")[, .(yyyymm, permno, dio)]
data[, dio := Winsorize(dio, quantile(dio, probs = c(0.01, 0.99)))]
to_file <- "../../../../data/institutional/s34_io_changes.RDS"
dir.create(dirname(to_file), recursive = T, showWarnings = F)
saveRDS(data, to_file)

# -- raw domestic equity fund flows by wficn
data <- readRDS("~/Desktop/J-Leaves/data/institutions/crspMutual/processed/fund_flows/all_fund_flows/quarterly_by_wficn/202405.RDS")[obj2 == "ED", .(yyyymm, wficn, tna_1, flow)]
to_file <- "../../../../data/institutional/fund_flows/quarterly_flow_by_wficn.RDS"
dir.create(dirname(to_file), recursive = T, showWarnings = F)
saveRDS(data, to_file)

# -- s12 mutual fund raw holdings
from_dir <- "~/Desktop/J-Leaves/data/institutions/13f/processed/s12/holdings_first_fdate_for_each_rdate/holdings/"
files <- list.files(from_dir)
files <- files[files >= "1980.RDS"]
to_dir <- "../../../../data/institutional/s12_holdings/"
dir.create(to_dir, recursive = T, showWarnings = F)
for (file in files) {
    tic(file)
    data <- readRDS(paste0(from_dir, file))
    saveRDS(data, paste0(to_dir, file))
    toc()
}

# -- quarterly FIT (used to be called quarterly_updated_20250313.RDS)
tmp <- readRDS("~/Desktop/J-Leaves/data/institutions/crspMutual/processed/FIT/redone_2025/earliestFdate/quarterly_holdingsFilledForward.RDS")
tmp <- tmp[frac_held_1 < 1] # sanity check
tmp <- tmp[, .(yyyymm, permno,
    fit = fit2shrout,
    fit_adj = fit2shrout_adj,
    fit_adj_cut01 = fit2shrout_adj_cut01,
    fit_adj_cut05 = fit2shrout_adj_cut05,
    fit_adj_cut10 = fit2shrout_adj_cut10
)]
tmp <- tmp[0 == rowSums(is.na(tmp))]

# tiny winsorization, almost none
vv <- names(tmp)
vv <- vv[grepl("fit", vv)]
for (this_v in vv) {
    setnames(tmp, this_v, "xx")
    tmp[, xx := Winsorize(xx, val = quantile(xx, probs = c(.0001, .9999)))]
    setnames(tmp, "xx", this_v)
}
rm(this_v, vv)

to_file <- "../../../../data/demand_shocks/j_fit/quarterly.RDS"
dir.create(dirname(to_file), recursive = T, showWarnings = F)
saveRDS(tmp, to_file)

# --- liquidity-related variables, downloaded from TAQ and processed
data <- readRDS("~/Desktop/J-Leaves/data/stockprices/processed/taq/liquidity/monthly_file.RDS")


tmp <- readRDS("~/Desktop/J-Leaves/data/stockprices/processed/taq/liquidity/daily_file.RDS")
