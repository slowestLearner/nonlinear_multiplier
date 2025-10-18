# === get the unexpected components of OFI in daily data

# load libraries
source('utilities/runmefirst.R')

# daily OFI data
data = readRDS('../../data/demand_shocks/ofi/daily.RDS')

# get daily lags
data[, idx := frank(date, ties.method = 'dense')] # create an index for trading day
for (i in 1:5){
  print(i)
  data = merge(data, data[, .(idx = idx+i, permno, xx = ofi)], by = c('idx','permno'))
  setnames(data, 'xx', paste0('ofi_', i))
}

# get weekly lags (assumed 5 days)
data[, ofi_1w := ofi_1 + ofi_2 + ofi_3 + ofi_4 + ofi_5]
for (i in 1:3){
  print(i)
  data = merge(data, data[, .(idx = idx+5*i, permno, xx = ofi_1w)], by = c('idx','permno'), all.x = T)
  setnames(data, 'xx', paste0('ofi_', (i+1), 'w'))
}

# get monthly lags (assumed 21 days)
data[, ofi_1m := ofi + ofi_1w + ofi_2w + ofi_3w + ofi_4w]
data[, ofi_1w := NULL]
for (i in 1:11){
  print(i)
  data = merge(data, data[, .(idx = idx+21*i, permno, xx = ofi_1m)], by = c('idx','permno'), all.x = T)
  setnames(data, 'xx', paste0('ofi_', (i+1), 'm'))
}
rm(i); gc()
data = data[0 == rowSums(is.na(data))]

# # == take a look at full sample
# ff = paste0('ofi ~ ', paste0(paste0('ofi_', 1:5), collapse = ' + '))
# ff = paste0(ff, ' + ', paste0(paste0('ofi_', 2:4, 'w'), collapse = ' + '))
# ff = paste0(ff, ' + ', paste0(paste0('ofi_', 2:12, 'm'), collapse = ' + '), ' | date')

# === estimate OFI residuals in a rolling window

# rergession formula
ff = paste0('ofi ~ ', paste0(paste0('ofi_', 1:5), collapse = ' + '))
ff = paste0(ff, ' + ', paste0(paste0('ofi_', 2:4, 'w'), collapse = ' + '))
ff = as.formula(paste0(ff, ' + ', paste0(paste0('ofi_', 2:12, 'm'), collapse = ' + ')))

# function to estimate residual for a subset of data
p.get_one = function(subdata){
  mm = lm(ff, subdata)
  subdata[, ofi_resid := mm$residuals]
  
  return(list(data.table(date = subdata[1, date], var = names(coef(mm)), coef = coef(mm)), 
              subdata[, .(date, permno, ofi, ofi_resid = mm$residuals)]))
}

data_list = split(data, by = 'date')

nc = parallel::detectCores() - 2
block_size = 500 # Adjust this number based on your memory/performance needs
num_blocks = ceiling(length(data_list) / block_size)

# on a computer with 6 cores, this takes around a min
out = list() # To store all results
for (i in 2:num_blocks) {
  print(i/num_blocks)
  tic()
  start_idx = (i - 1) * block_size + 1
  end_idx = min(i * block_size, length(data_list))
  out = c(out, parallel::mclapply(data_list[start_idx:end_idx], p.get_one, mc.cores = nc)); gc()
  toc()
}
rm(i, num_blocks, block_size, start_idx, end_idx, ff)
stopifnot(length(out) == length(data_list))

# coef estimates
coefdata = rbindlist(lapply(out, `[[`, 1))
ofidata = rbindlist(lapply(out, `[[`, 2))
ofidata[, ofi := NULL]

dir.create('tmp/raw_data/cleaning/ofi/', showWarnings = F, recursive = T)
saveRDS(coefdata, 'tmp/raw_data/cleaning/ofi/fm_regression_coefs.RDS')
saveRDS(ofidata, 'tmp/raw_data/cleaning/ofi/residual_daily_ofi.RDS')

# # === SANITY: check with earlier data

# # coefs
# old = readRDS('../20250117_quarterly/tmp/raw_data/cleaning/ofi/fm_regression_coefs.RDS')
# new = readRDS('tmp/raw_data/cleaning/ofi/fm_regression_coefs.RDS')
# mean(old == new)

# # residuals - no problem here
# old = readRDS('../20250117_quarterly/tmp/raw_data/cleaning/ofi/residual_daily_ofi.RDS')
# new = readRDS('tmp/raw_data/cleaning/ofi/residual_daily_ofi.RDS')
# dim(old) == dim(new)
# mean(old == new)
# compare = merge(old, new, by = c('date','permno')); rm(old,new)
# compare[(date == '1994-03-25') & (permno == 10119)]
