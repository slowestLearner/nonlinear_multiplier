# --- Compute SD, HHI, and cumulative past demand for dynamic demand over the previous H quarters
library(this.path)
setwd(this.path::this.dir())
source("../../utilities/runmefirst.R")
source("../../utilities/regressions.R")

data <- readRDS("../../tmp/raw_data/reg_inputs/reg_table_static.RDS")[type %in% c("OFI", "FIT"), .(yyyymm, permno, type, ofi)]

# get lags
data[, idx := frank(yyyymm, ties.method = "dense")]
for (i in 1:12) {
  data <- merge(data, data[, .(idx = idx + i, permno, type, xx = ofi)], by = c("idx", "permno", "type"), all.x = T)
  setnames(data, "xx", paste0("ofi_", i))
}
data[, c("idx", "ofi") := NULL]

# make long table
data <- melt(data, id.vars = c("yyyymm", "permno", "type"), variable.name = "lag", value.name = "ofi") %>%
  mutate(lag = as.integer(gsub("ofi_", "", lag))) %>%
  setDT()

# function to compute modified hhi
p.hhi <- function(x) {
  x_sum <- sum(x)
  dir <- sign(x_sum)
  x_tilde <- x * (sign(x) == dir)
  w_tilde <- x_tilde / sum(x_tilde)
  hhi <- sum(w_tilde^2)
  return(hhi)
}

# compute concentration by lookback horizon
out <- data.table()
for (this_lag in c(4, 8, 12)) {
  tic(this_lag)
  tmp <- data[lag %in% 1:this_lag, .(yyyymm, permno, type, lag, ofi)] %>% na.omit()
  tmp <- tmp[, .(lag = this_lag, cumofi_lag = sum(ofi), obs = .N, ofi_sd = sd(ofi), hhi_pow2 = p.hhi(ofi)), .(yyyymm, type, permno)] %>%
    filter(obs == this_lag) %>%
    select(-obs) %>%
    na.omit()
  out <- rbind(out, tmp)
  toc()
}
# out_bk <- copy(out)

out <- melt(out, id.vars = c("yyyymm", "permno", "type", "lag", "cumofi_lag"), variable.name = "disp_type", value.name = "disp") %>%
  mutate(disp_type = as.character(disp_type)) %>%
  setDT()

to_file <- "tmp/raw_files/concentration_of_lagged_demand.RDS"
dir.create(dirname(to_file), showWarnings = FALSE, recursive = TRUE)
saveRDS(out, to_file)

# # --- sanity: computation is correct

# new <- readRDS("tmp/raw_files/concentration_of_lagged_demand.RDS")
# new <- new[, .(new = last(cumofi_lag)), .(yyyymm, permno, type, lag)]
# new <- new[lag == 4][, lag := NULL]

# old <- readRDS("../../tmp/raw_data/reg_inputs/reg_table_dynamic.RDS")
# old <- old[, .(yyyymm, type, permno, old = cumofi_4)] %>% na.omit()
# old <- old[type %in% c("OFI", "FIT")]

# out <- merge(new, old, by = c("yyyymm", "permno", "type"))
# feols(new ~ old, out)
