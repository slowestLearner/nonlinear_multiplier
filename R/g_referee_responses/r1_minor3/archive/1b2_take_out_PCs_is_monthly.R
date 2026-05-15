# --- BMI uses monthly data. Trick: I am going to use the whole quarter to do PCA, but then just subtract in monthly data
library(this.path)
setwd(this.path::this.dir())
source("../../utilities/runmefirst.R")
library(RSpectra) # Faster for extracting top K components
options(width = 200)

# these are the stocks we use. append with bmi returns if needed
stock_data <- readRDS("../../../../../data/stocks/prices/monthly_return.RDS")[, .(yyyymm, permno, ret)]

# restrict to bmi months
tmp <- readRDS("../../tmp/raw_data/reg_inputs/reg_table_static.RDS")[type == "BMI", .(yyyymm)] %>% unique()
stock_data <- merge(stock_data, tmp, by = "yyyymm")
rm(tmp)

# append bmi returns if applicable
tmp <- readRDS("../../tmp/raw_data/reg_inputs/reg_table_static.RDS")[type == "BMI", .(yyyymm, permno, ret_bmi = ret)] %>% unique()
stock_data <- merge(stock_data, tmp, by = c("yyyymm", "permno"), all = T)
stock_data[, ret := ifelse(!is.na(ret_bmi), ret_bmi, ret)]
stock_data[, ret_bmi := NULL]
rm(tmp)

# get daily returns
from_dir <- "../../../../../../../../Desktop/J-Leaves/data/stockprices/raw/daily/by_decade/"
files <- list.files(from_dir)
files <- files[files >= "1990s.RDS"]
p.get_one <- function(this_file) {
    # this_file <- files[1]
    data <- readRDS(paste0(from_dir, this_file))
    data <- data[, .(date, permno, ret)]
    return(data)
}

tic()
data <- rbindlist(mclapply(files, p.get_one, mc.cores = detectCores() - 2))
rm(from_dir, files, p.get_one)
toc()

tmp <- unique(data[, .(date)])
tmp[, yyyymm := 100 * year(date) + 3 * quarter(date)]
tmp <- tmp[yyyymm %in% stock_data[, unique(yyyymm)]]
data <- merge(data, tmp, by = "date")
rm(tmp)

# create a grid
grid <- stock_data[, .(yyyymm, permno)]
ym_to_date <- unique(data[, .(yyyymm, date)])
grid <- merge(grid, ym_to_date, by = "yyyymm", all = T, allow.cartesian = T)
data <- merge(grid, data, by = c("yyyymm", "date", "permno"), all.x = T)
rm(grid, ym_to_date)

# fill in some missing returns
data[, mean_ret := mean(ret, na.rm = T), .(yyyymm)]
data[is.na(ret), ret := mean_ret]
data[, mean_ret := NULL]

# for each year-quarter, use daily data to estimate a PCA, and then use quarter data in stock_data to get the residual returns
p.get_one_period <- function(this_data, k_list = c(1, 3, 5, 10, 15, 20)) {
    # 1. Identity the quarter
    this_ym <- this_data[1, yyyymm]

    # 2. Pivot to Wide Matrix (T days x N stocks)
    # dcast is efficient here; it will create a 'date' column and one column per permno
    wide_data <- dcast(this_data, date ~ permno, value.var = "ret")

    # Convert to matrix, dropping the 'date' column (column 1)
    # We use raw returns (Option 1) to allow PCs to capture the quarterly drift
    mat <- as.matrix(wide_data[, -1, with = FALSE])

    # Store permnos as integers to ensure matching works with stock_data
    permnos_in_pca <- as.integer(colnames(mat))

    # 3. PCA via SVD (RSpectra handles the N > T case efficiently)
    # We extract the max k needed
    max_k <- max(k_list)
    svd_decomp <- RSpectra::svds(mat, k = max_k)

    # f_daily (T x K): The realized daily factor returns (PC scores)
    # These are derived from U %*% Sigma
    f_daily <- svd_decomp$u %*% diag(svd_decomp$d)

    # f_quarterly (1 x K): The total realized factor return for the quarter
    # Since we didn't center, these will now be non-zero values
    f_quarterly <- colSums(f_daily)

    # loadings (N x K): The beta/exposure of each stock to the factors
    # This is the V matrix from the SVD
    loadings <- svd_decomp$v

    # 4. Align with actual quarterly returns from stock_data
    # Pull the target quarterly returns for this specific yyyymm
    q_rets <- stock_data[yyyymm == this_ym]

    # Use match to ensure the vector of quarterly returns aligns 1:1 with the PCA matrix columns
    idx <- match(permnos_in_pca, q_rets$permno)
    actual_ret_vec <- q_rets$ret[idx]

    # 5. Calculate Residuals
    # Initialize the results table
    results <- data.table(permno = permnos_in_pca, yyyymm = this_ym)

    for (k in k_list) {
        # Predicted Quarterly Return = Loadings * Quarterly Factor Realizations
        # We subset the loadings and factor returns to the first k components
        predicted_ret <- as.vector(loadings[, 1:k, drop = FALSE] %*% f_quarterly[1:k])

        # Residual = Actual Quarterly - Predicted Quarterly
        col_name <- paste0("res_pc", k)
        results[, (col_name) := actual_ret_vec - predicted_ret]
    }

    # Return results, filtering out stocks that didn't have a quarterly return record
    return(results[!is.na(actual_ret_vec)])
}

tic()
out <- rbindlist(mclapply(split(data, by = "yyyymm"), p.get_one_period, mc.cores = detectCores() - 2))
toc()

stopifnot(nrow(out) == nrow(stock_data))

# save
to_file <- "tmp/pca_residuals/monthly.RDS"
dir.create(dirname(to_file), showWarnings = FALSE, recursive = TRUE)
saveRDS(out, to_file)
