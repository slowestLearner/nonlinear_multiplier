# --- Use alternative least squares to do PCA on fund flows
# TODEL: different from earlier where i required the characteristics-based residuals to exist. We no longer do that, so the resulting panel of flows became larger
library(softImpute)
library(this.path)
setwd(this.path::this.dir())
source("../utilities/runmefirst.R")

# load raw fund flows
data <- readRDS("../../../../data/institutional/fund_flows/quarterly_flow_by_wficn.RDS")

# worker function. Do PCA with "this_rank"
p.get_one <- function(this_rank) {
  fit <- softImpute(flow_matrix, rank.max = this_rank, lambda = 0, type = "als", maxit = 500)
  flows_imputed <- fit$u %*% (fit$d * t(fit$v))
  flows_imputed <- data.table(flows_imputed)
  names(flows_imputed) <- colnames(flow_matrix)
  flows_imputed[, yyyymm := yms]
  flows_imputed <- melt(flows_imputed, id.vars = "yyyymm", variable.name = "wficn", value.name = "flow_resid_imputed")
  flows_imputed[, wficn := as.integer(as.character(wficn))]
  flows_imputed <- merge(flows_imputed, data[, .(yyyymm, wficn, flow)], by = c("yyyymm", "wficn"))
  flows_imputed <- flows_imputed[, .(yyyymm, wficn, resid = flow - flow_resid_imputed)]
  setnames(flows_imputed, "resid", paste0("flow_took_out_", this_rank, "_pcs"))
  return(flows_imputed)
}

# do PCA with different number of PCs
out_all <- data.table()
flow_matrix <- dcast(data, yyyymm ~ wficn, value.var = "flow")
yms <- flow_matrix[, yyyymm]
flow_matrix[, yyyymm := NULL]
flow_matrix <- as.matrix(flow_matrix)

# define PCA ranks
ranks <- c(1, 2, 3, 5, 7, 10)

# run PCA in parallel. On a computer w/ 6 cores, takes around 20 secs total
tic()
plan(multisession, workers = detectCores() - 2)
results_list <- future_lapply(ranks, function(this_rank) {
  out <- p.get_one(this_rank)
  return(out)
}, future.seed = 123, future.packages = "data.table")
plan(sequential)

# Merge the results back
for (i in seq_along(ranks)) {
  data <- merge(data, results_list[[i]], by = c("yyyymm", "wficn"), all.x = TRUE)
}
toc()

# save
to_file <- "../tmp/additional/clean_fit/flow_residuals/1_pca_flow_residuals.RDS"
dir.create(dirname(to_file), showWarnings = FALSE, recursive = TRUE)
saveRDS(data, to_file)
