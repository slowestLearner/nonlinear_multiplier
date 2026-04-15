# --- compute unexpected components of OFI in daily data (pre-whitenening)
# TODO: we don't want this in the final version. Another thing is - this pre-whitenening step is sort of gangster. It runs a separate regression for each day...
library(this.path)
setwd(this.path::this.dir())
source("../utilities/runmefirst.R")

# daily OFI data
data <- readRDS("../../../../data/demand_shocks/ofi/daily.RDS")

# get 5 daily lags
data[, idx := frank(date, ties.method = "dense")] # create an index for trading day
for (i in 1:5) {
    print(i)
    data <- merge(data, data[, .(idx = idx + i, permno, xx = ofi)], by = c("idx", "permno"), all.x = T)
    setnames(data, "xx", paste0("ofi_", i))
}

# get weekly lags (each week is 5 trading days)
data[, ofi_1w := ofi_1 + ofi_2 + ofi_3 + ofi_4 + ofi_5]
for (i in 1:3) {
    print(i)
    data <- merge(data, data[, .(idx = idx + 5 * i, permno, xx = ofi_1w)], by = c("idx", "permno"), all.x = T)
    setnames(data, "xx", paste0("ofi_", (i + 1), "w"))
}

# get monthly lags (each month is 21 trading days)
data[, ofi_1m := ofi + ofi_1w + ofi_2w + ofi_3w + ofi_4w]
data[, ofi_1w := NULL]
for (i in 1:11) {
    print(i)
    data <- merge(data, data[, .(idx = idx + 21 * i, permno, xx = ofi_1m)], by = c("idx", "permno"), all.x = T)
    setnames(data, "xx", paste0("ofi_", (i + 1), "m"))
}
rm(i)
data[, idx := NULL]
data <- data %>% na.omit()
gc()

# --- estimate OFI residuals by day


# worker function to estimate by day. return both residuals and coefficients
p.get_one <- function(subdata) {
    # regression formula
    ff <- paste0("ofi ~ ", paste0(paste0("ofi_", 1:5), collapse = " + "))
    ff <- paste0(ff, " + ", paste0(paste0("ofi_", 2:4, "w"), collapse = " + "))
    ff <- as.formula(paste0(ff, " + ", paste0(paste0("ofi_", 2:12, "m"), collapse = " + ")))

    mm <- lm(ff, subdata)
    subdata[, ofi_resid := mm$residuals]

    return(list(
        data.table(date = subdata[1, date], var = names(coef(mm)), coef = coef(mm)),
        subdata[, .(date, permno, ofi, ofi_resid = mm$residuals)]
    ))
}

# process by date
data_list <- split(data, by = "date")
rm(data)
gc()

nc <- parallel::detectCores() - 2
block_size <- 1000
tmp <- data.table(idx = 1:length(data_list))
tmp[, block_idx := ceiling(idx / block_size)]

# why is this initially so slow?

# on a computer with 6 cores, this takes around a min
out <- list() # To store all results
for (i in unique(tmp[, block_idx])) {
    print(i / max(tmp[, block_idx]))
    tic()
    # start_idx <- (i - 1) * block_size + 1
    # end_idx <- min(i * block_size, length(data_list))
    # out <- c(out, parallel::mclapply(data_list[start_idx:end_idx], p.get_one, mc.cores = nc))
    out <- c(out, parallel::mclapply(data_list[tmp[block_idx == i, idx]], p.get_one, mc.cores = nc))
    gc()
    toc()
}
rm(i, block_size, p.get_one, tmp)
stopifnot(length(out) == length(data_list))

# combine coef estimates and residuals
coefdata <- rbindlist(lapply(out, `[[`, 1))
ofidata <- rbindlist(lapply(out, `[[`, 2))
ofidata[, ofi := NULL]

# save coef estimates
to_file <- "../tmp/raw_data/cleaning/ofi/fm_regression_coefs.RDS"
dir.create(dirname(to_file), showWarnings = FALSE, recursive = TRUE)
saveRDS(coefdata, to_file)

# save ofi residuals
to_file <- "../tmp/raw_data/cleaning/ofi/residual_daily_ofi.RDS"
dir.create(dirname(to_file), showWarnings = FALSE, recursive = TRUE)
saveRDS(ofidata, to_file)

# # # === SANITY: sufficiently similar

# # coefs: sufficient similar

# # residuals - end up with slightly more data now, after fixing a bug
# new <- readRDS('../tmp/raw_data/cleaning/ofi/residual_daily_ofi.RDS')
# old <- readRDS('../tmp/raw_data/cleaning/ofi_todel/residual_daily_ofi.RDS')
# new <- new[date %in% unique(old$date)]
# tt <- merge(new, old, by = c('date','permno'))
# tt[, cor(ofi_resid.x, ofi_resid.y)]
