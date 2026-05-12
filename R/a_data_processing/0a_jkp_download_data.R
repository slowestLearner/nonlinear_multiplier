# --- download needed data from WRDS
library(this.path)
setwd(this.path::this.dir())
source("../utilities/runmefirst.R")
library(getPass)
library(RPostgres)

# --- JPK characteristics

# ----- LOGIN TO WRDS
user <- getPass("wrds username")
pass <- getPass("wrds password")

tic()
wrds <- dbConnect(Postgres(),
    host = "wrds-pgdata.wharton.upenn.edu",
    port = 9737,
    dbname = "wrds",
    user = user,
    password = pass,
    sslmode = "require"
)
rm(pass, user)
toc()

# === Download all characteristics from JKP

# gosh this can take 50 mins to download
tic()
data <- dbSendQuery(
    conn = wrds, statement =
        "SELECT *
                      FROM contrib.global_factor WHERE primary_sec = 1 AND common = 1 AND exch_main = 1 AND obs_main = 1 AND excntry = 'USA' AND permno IS NOT NULL AND date >= '1991-01-01'"
) %>%
    dbFetch(n = -1) %>%
    setDT()
toc()
gc()

# these are not characteristics
vars_todel <- c("id", "permco", "gvkey", "iid", "excntry", "exch_main", "common", "primary_sec", "bidask", "crsp_shrcd", "crsp_exchcd", "comp_tpci", "comp_exchg", "curcd", "fx", "date", "adjfct", "shares", "me", "me_company", "prc", "prc_local", "prc_high", "prc_low", "dolvol", "tvol", "ret", "ret_local", "ret_exc", "ret_lag_dif", "source_crsp", "ret_exc_lead1m", "obs_main", "gics", "sic", "naics", "ff49", "size_grp", "market_equity", "div1m_me")
data[, (vars_todel) := NULL]

# change month column, get quarterly
tmp <- unique(data[, .(eom)])[order(eom)]
tmp[, yyyymm := 100 * year(eom) + month(eom)]
data <- merge(data, tmp, by = "eom")
data[, eom := yyyymm][, yyyymm := NULL]
setnames(data, "eom", "yyyymm")

# just quarterly
tmp <- unique(data[, .(yyyymm)])
tmp[, mm := yyyymm - 100 * floor(yyyymm / 100)]
tmp <- tmp[mm %in% c(3, 6, 9, 12)][, mm := NULL]
data <- merge(data, tmp, by = "yyyymm")
rm(tmp)
gc()

# save into 10 pieces
vv <- names(data)
vv <- setdiff(vv, c("yyyymm", "permno"))

tmp <- data.table(var = vv)
tmp[, idx := .I]
tmp[, bin := ntile(idx, 10)]

to_dir <- "../tmp/raw_data/jkp_chars_not_lagged/0_raw/"
dir.create(to_dir, recursive = T, showWarnings = F)

for (i in 1:10) {
    tic(i)
    to_file <- paste0(to_dir, "part_", i, ".RDS")
    saveRDS(data[, c("yyyymm", "permno", tmp[bin == i, var]), with = F], to_file)
    toc()
}
