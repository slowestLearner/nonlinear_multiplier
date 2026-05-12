# --- for each quarter, use previous 1y data to estimate
library(this.path)
setwd(this.path::this.dir())
source("../../utilities/runmefirst.R")
library(RSpectra) # Faster for extracting top K components
options(width = 200)

# these are the stocks we use.
stock_data <- readRDS("../../../../../data/stocks/prices/quarterly_return.RDS")[, .(yyyymm, permno, ret)]
stock_data <- stock_data[yyyymm %in% 199303:202212]

# get daily returns
tic()
from_dir <- "../../../../../../../../Desktop/J-Leaves/data/stockprices/raw/daily/by_decade/"
files <- list.files(from_dir)
files <- files[files >= "1990s.RDS"]
p.get_one <- function(this_file) {
    data <- readRDS(paste0(from_dir, this_file))
    data <- data[, .(date, permno, ret)]
    return(data)
}

data <- rbindlist(mclapply(files, p.get_one, mc.cores = detectCores() - 2))
rm(from_dir, files, p.get_one)
toc()

# turn into yyyymm and idx
tmp <- unique(data[, .(date)])
tmp[, yyyymm := 100 * year(date) + 3 * quarter(date)]
tmp[, idx := frank(yyyymm, ties.method = "dense")]
stock_data <- merge(stock_data, unique(tmp[, .(yyyymm, idx)]), by = "yyyymm")
data <- merge(data, tmp, by = "date")
rm(tmp)

# ---- stopped here

# 1. Prepare a clean daily dataset for the worker
# We still need the grid/imputation to handle the PCA requirements
# (Note: Doing this globally or per-block is a memory vs. speed trade-off)
message("Building grid and imputing...")
grid <- stock_data[, .(idx, permno)]
full_grid <- copy(grid)
for (i in 1:4) {
    grid[, idx := idx - 1]
    full_grid <- rbind(full_grid, grid)
}
grid <- unique(full_grid)
rm(full_grid, i)

ym_to_date <- unique(data[, .(yyyymm, date, idx)])
grid <- merge(grid, ym_to_date, by = c("idx"), all.x = T, allow.cartesian = T)
data <- merge(grid, data, by = c("yyyymm", "date", "idx", "permno"), all.x = T)
gc()

data[, mean_ret := mean(ret, na.rm = T), .(date)] # Daily cross-sectional mean
data[is.na(ret), ret := mean_ret]
data[, mean_ret := NULL]
rm(grid)



# grid <- stock_data[, .(yyyymm, idx, permno)]
# ym_to_date <- unique(data[, .(yyyymm, date, idx)])
# grid <- merge(grid, ym_to_date, by = c("yyyymm", "idx"), all = T, allow.cartesian = T)
# data <- merge(grid, data, by = c("yyyymm", "date", "idx", "permno"), all.x = T)
# data[, mean_ret := mean(ret, na.rm = T), .(date)] # Daily cross-sectional mean
# data[is.na(ret), ret := mean_ret]
# data[, mean_ret := NULL]
# rm(grid)

# 2. The Rolling Worker Function
p.get_rolling_residuals <- function(curr_idx, k_list = c(1, 3, 5, 10, 15, 20)) {
    # A. Get Training Data (Previous 4 quarters)
    # We use these to estimate the PCA LOADINGS
    train_data <- data[idx %in% (curr_idx - 4):(curr_idx - 1)]

    wide_train <- dcast(train_data, date ~ permno, value.var = "ret")
    wide_train[is.na(wide_train)] <- 0 # FIX
    mat_train <- as.matrix(wide_train[, -1, with = FALSE])
    permnos_train <- as.integer(colnames(mat_train))

    # B. Get Current Data (The target quarter)
    # We use these to find the current FACTOR REALIZATIONS
    test_data <- data[idx == curr_idx]
    wide_test <- dcast(test_data, date ~ permno, value.var = "ret")
    mat_test <- as.matrix(wide_test[, -1, with = FALSE])
    permnos_test <- as.integer(colnames(mat_test))

    # C. Alignment: Only keep stocks present in both periods
    common_permnos <- intersect(permnos_train, permnos_test)
    mat_train <- mat_train[, which(permnos_train %in% common_permnos)]
    mat_test <- mat_test[, which(permnos_test %in% common_permnos)]


    # D. Estimate PCA on Training Data
    max_k <- max(k_list)
    # V are the loadings (N_stocks x K)
    svd_train <- RSpectra::svds(mat_train, k = max_k)
    loadings_H <- svd_train$v

    # E. Project Current Returns onto Historical Loadings
    # Factors_Daily = Returns_Current %*% Loadings_Historical
    # This tells us how the "Historical Factor Portfolios" performed TODAY
    f_daily_curr <- mat_test %*% loadings_H
    f_quarterly_curr <- colSums(f_daily_curr)

    # F. Compute Residuals for the Current Quarter
    q_rets_curr <- stock_data[idx == curr_idx]
    idx_match <- match(common_permnos, q_rets_curr$permno)
    actual_q_ret <- q_rets_curr$ret[idx_match]

    this_ym <- unique(test_data$yyyymm)
    results <- data.table(permno = common_permnos, yyyymm = this_ym, idx = curr_idx)

    for (k in k_list) {
        # Systematic = Loadings_H * Factors_Q_Realized
        predicted_ret <- as.vector(loadings_H[, 1:k, drop = FALSE] %*% f_quarterly_curr[1:k])
        results[, (paste0("res_pc", k)) := actual_q_ret - predicted_ret]
    }

    results[, idx := NULL]
    return(results[!is.na(actual_q_ret)])
}

# 3. Execution
tic()
target_indices <- sort(unique(stock_data$idx))
out <- rbindlist(mclapply(target_indices, p.get_rolling_residuals, mc.cores = detectCores() - 2))
toc()

stopifnot(nrow(out) == nrow(stock_data))

# cor(out[, -c(1:2), with = F])
# cov(out[, -c(1:2), with = F])

# save
to_file <- "tmp/pca_residuals/quarterly_oos.RDS"
dir.create(dirname(to_file), showWarnings = FALSE, recursive = TRUE)
saveRDS(out, to_file)
