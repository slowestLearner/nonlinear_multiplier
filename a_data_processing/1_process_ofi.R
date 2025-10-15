# === get the unexpected components of OFI

# load libraries
library(data.table)

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



# PAUSE

data = data[0 == rowSums(is.na(data))]

# == take a look at full sample
ff = paste0('ofi ~ ', paste0(paste0('ofi_', 1:5), collapse = ' + '))
ff = paste0(ff, ' + ', paste0(paste0('ofi_', 2:4, 'w'), collapse = ' + '))
ff = paste0(ff, ' + ', paste0(paste0('ofi_', 2:12, 'm'), collapse = ' + '), ' | date')

# === estimate in a rolling window
ff = paste0('ofi ~ ', paste0(paste0('ofi_', 1:5), collapse = ' + '))
ff = paste0(ff, ' + ', paste0(paste0('ofi_', 2:4, 'w'), collapse = ' + '))
ff = as.formula(paste0(ff, ' + ', paste0(paste0('ofi_', 2:12, 'm'), collapse = ' + ')))

p.get_one = function(subdata){
  mm = lm(ff, subdata)
  subdata[, ofi_resid := mm$residuals]
  
  return(list(data.table(date = subdata[1, date], var = names(coef(mm)), coef = coef(mm)), 
              subdata[, .(date, permno, ofi, ofi_resid = mm$residuals)]))
}

data_list = split(data, by = 'date')

nc = detectCores() - 2

block_size = 500 # Adjust this number based on your memory/performance needs
num_blocks = ceiling(length(data_list) / block_size)

# the whole thing takes around a min
out = list() # To store all results
for (i in 1:num_blocks) {
  print(i/num_blocks)
  tic()
  start_idx = (i - 1) * block_size + 1
  end_idx = min(i * block_size, length(data_list))
  out = c(out, mclapply(data_list[start_idx:end_idx], p.get_one, mc.cores = nc)); gc()
  toc()
}
stopifnot(length(out) == length(data_list))

# coef estimates
coefdata = rbindlist(lapply(out, `[[`, 1))
ofidata = rbindlist(lapply(out, `[[`, 2))
ofidata[, ofi := NULL]

dir.create('tmp/raw_data/cleaning/ofi/', showWarnings = F, recursive = T)
saveRDS(coefdata, 'tmp/raw_data/cleaning/ofi/fm_regression_coefs.RDS')
saveRDS(ofidata, 'tmp/raw_data/cleaning/ofi/residual_daily_ofi.RDS')

# === let's take a look at the coefs
coefdata = readRDS('tmp/raw_data/ofi_processed/fm_regression_coefs.RDS')
coefdata = coefdata[var != '(Intercept)']
out = coefdata[, .(coef = mean(coef), se = sd(coef)/sqrt(length(.I))), var]
out[, hor := c(1:5, 5*2:4 - 3, 21*2:12 - 11)]

ggplot(out, aes(x = hor, y = coef)) + geom_point() + geom_line() + 
  geom_ribbon(aes(ymin = coef-2*se, ymax = coef+2*se), alpha = .2) + 
  scale_x_continuous(trans = 'log10')

