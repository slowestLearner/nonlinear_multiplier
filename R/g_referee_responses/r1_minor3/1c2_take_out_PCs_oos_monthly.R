# --- BMI Monthly Rolling Inference
library(this.path)
setwd(this.path::this.dir())
source("../../utilities/runmefirst.R")
library(RSpectra)
options(width = 200)

# 1. Prepare Stock Data (Monthly)
stock_data <- readRDS("../../../../../data/stocks/prices/monthly_return.RDS")[, .(yyyymm, permno, ret)]

# Filter for BMI months and merge BMI-specific returns
bmi_ref <- readRDS("../../tmp/raw_data/reg_inputs/reg_table_static.RDS")[type == "BMI"]
tmp_months <- unique(bmi_ref[, .(yyyymm)])
stock_data <- merge(stock_data, tmp_months, by = "yyyymm")

# Update with BMI-specific returns where applicable
tmp_bmi <- unique(bmi_ref[, .(yyyymm, permno, ret_bmi = ret)])
stock_data <- merge(stock_data, tmp_bmi, by = c("yyyymm", "permno"), all = T)
stock_data[, ret := ifelse(!is.na(ret_bmi), ret_bmi, ret)]
stock_data[, ret_bmi := NULL]
rm(bmi_ref, tmp_months, tmp_bmi)

# 2. Get Daily Returns & Standardize Index
tic("Loading Daily Data")
from_dir <- "../../../../../../../../Desktop/J-Leaves/data/stockprices/raw/daily/by_decade/"
files <- list.files(from_dir, pattern = ".RDS")
files <- files[files >= "1990s.RDS"]

p.get_one <- function(this_file) {
    data <- readRDS(paste0(from_dir, this_file))
    return(data[, .(date, permno, ret)])
}
data <- rbindlist(mclapply(files, p.get_one, mc.cores = detectCores() - 2))
toc()

# Create Monthly Index (idx)
tmp_dates <- unique(data[, .(date)])
tmp_dates[, yyyymm := 100 * year(date) + month(date)] # True monthly buckets
tmp_dates[, idx := frank(yyyymm, ties.method = "dense")]

# Apply index to both datasets
data <- merge(data, tmp_dates, by = "date")
stock_data <- merge(stock_data, unique(tmp_dates[, .(yyyymm, idx)]), by = "yyyymm")
rm(tmp_dates)

# 3. Clean Imputation (Daily Cross-Sectional Mean)
message("Imputing daily data...")
data[, mean_ret := mean(ret, na.rm = T), .(date)]
data[is.na(ret), ret := mean_ret]
data[, mean_ret := NULL]
setkey(data, idx, permno)

# 4. The Monthly Rolling Worker
p.get_monthly_rolling_residuals <- function(curr_idx, k_list = c(1, 3, 5, 10, 15, 20)) {
    # A. Training Data: Last 12 Months
    train_data <- data[idx %in% (curr_idx - 12):(curr_idx - 1)]
    if (nrow(train_data) == 0) {
        return(NULL)
    }

    wide_train <- dcast(train_data, date ~ permno, value.var = "ret")
    mat_train <- as.matrix(wide_train[, -1, with = FALSE])
    mat_train[is.na(mat_train)] <- 0
    permnos_train <- as.integer(colnames(mat_train))

    # B. Test Data: Current Month
    test_data <- data[idx == curr_idx]
    wide_test <- dcast(test_data, date ~ permno, value.var = "ret")
    mat_test <- as.matrix(wide_test[, -1, with = FALSE])
    mat_test[is.na(mat_test)] <- 0
    permnos_test <- as.integer(colnames(mat_test))

    # C. Alignment & SVD
    common_permnos <- intersect(permnos_train, permnos_test)
    max_k <- max(k_list)
    if (length(common_permnos) <= max_k) {
        return(NULL)
    }

    mat_train <- mat_train[, match(common_permnos, permnos_train)]
    mat_test <- mat_test[, match(common_permnos, permnos_test)]

    # D. PCA on Trailing Window
    svd_train <- RSpectra::svds(mat_train, k = max_k)
    loadings_H <- svd_train$v

    # E. Projection to find Monthly Factor Realization
    f_daily_curr <- mat_test %*% loadings_H
    f_monthly_curr <- colSums(f_daily_curr) # The "Systematic" return for this month

    # F. Compute Residuals
    q_rets_curr <- stock_data[idx == curr_idx]
    idx_match <- match(common_permnos, q_rets_curr$permno)
    actual_m_ret <- q_rets_curr$ret[idx_match]

    this_ym <- unique(test_data$yyyymm)
    results <- data.table(permno = common_permnos, yyyymm = this_ym)

    for (k in k_list) {
        predicted_ret <- as.vector(loadings_H[, 1:k, drop = FALSE] %*% f_monthly_curr[1:k])
        results[, (paste0("res_pc", k)) := actual_m_ret - predicted_ret]
    }

    return(results[!is.na(actual_m_ret)])
}

# 5. Execution
tic("Rolling PCA Execution")
target_indices <- sort(unique(stock_data$idx))
out <- rbindlist(mclapply(target_indices, p.get_monthly_rolling_residuals, mc.cores = detectCores() - 2))
toc()

# 6. Save
to_file <- "tmp/pca_residuals/monthly_oos.RDS"
dir.create(dirname(to_file), showWarnings = FALSE, recursive = TRUE)
saveRDS(out, to_file)
