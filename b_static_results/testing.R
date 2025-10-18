# formulas for controls
control_formulas <- c(
  "1",
  paste0(c(controls_char), collapse = "+"),
  paste0(c(controls_char, controls_liq), collapse = "+")
)

# process one type of data
p.process_one_type <- function(data) {
  print(this_type)

  data <- copy(data_all[type == this_type])
  this_type <- data[1, type]

  out_all <- data.table() # save all results here

  # introduce direct controls
  for (spec_idx in 1:length(control_formulas)) {
    ff <- paste0("ret ~ ", control_formulas[spec_idx], " + ofi_bin1 + ofi_bin2 + ofi_bin3")
    if (this_type == "BMI") {
      ff <- paste0(ff, "+", paste0(controls_bmi, collapse = "+"))
    }
    out <- p.fama_macbeth(data, ff, compare_coefs = TRUE)
    out[, spec_idx := spec_idx]
    out_all <- rbind(out_all, out)
  }

  # then add interactions with ofi
  for (this_v in c(controls_char, controls_liq)) {
    # print(this_v)

    setnames(data, this_v, "xx")
    data[, yy := xx * ofi]
    setnames(data, c("xx", "yy"), c(this_v, paste0("ofi_", this_v)))
    ff <- paste0(ff, " + ofi_", this_v)
    spec_idx <- spec_idx + 1

    out <- p.fama_macbeth(data, ff, compare_coefs = TRUE)
    out[, spec_idx := spec_idx]
    out_all <- rbind(out_all, out)
  }

  # name the control variables being added
  tmp <- data.table(
    idx = 1:(length(controls_list) + 4),
    var_added = c("none_init", "controls_char", "controls_char+controls_liq", "none", controls_list),
    var_type = c(
      rep("", 4), rep("return-predicting chars", length(controls_char)),
      rep("liquidity", length(controls_liq))
    )
  )
  out_all <- merge(out_all, tmp, by = "idx", all.x = T)
  rm(tmp)
  out_all[, type := this_type]

  return(out_all)
}


# process all types in parallel
types <- unique(data_all[, type])
out_all <- rbindlist(mclapply(types, p.process_one_type, mc.cores = nc))

dir.create("tmp/price_impact/regression_contemp/", recursive = T, showWarnings = F)
saveRDS(out_all, "tmp/price_impact/regression_contemp/full_sample.RDS")

# check coefs
range(out_all[(var == "ofi_absofi") & (type != "OFI"), coef])
range(out_all[(var == "ofi_absofi") & (type != "OFI"), coef / se])

# == FM with stdev-based bin interactions. Also keep track of the cutoffs

out_all <- data.table()
for (this_type in unique(data_all[, type])) {
  print(this_type)

  tic()
  data <- copy(data_all[type == this_type])

  # bins
  data <- merge(data, data[, .(ss = sd(ofi)), yyyymm], by = "yyyymm")
  data[, bin := 1]
  data[abs(ofi) > ss, bin := 2]
  data[abs(ofi) > 2 * ss, bin := 3]
  data[, ofi_bin1 := ifelse(bin == 1, ofi, 0)]
  data[, ofi_bin2 := ifelse(bin == 2, ofi, 0)]
  data[, ofi_bin3 := ifelse(bin == 3, ofi, 0)]

  # function to do Fama-MacBeth, and also to record coef differences
  p.get_one_fm <- function(ff) {
    p.get_one_period <- function(this_ym) {
      ols <- lm(paste0("ret ~ ofi_bin1 + ofi_bin2 + ofi_bin3 + ", ff), data[yyyymm == this_ym])
      ols_no_ofi <- lm(paste0("ret ~ ", ff), data[yyyymm == this_ym])
      out <- data.table(
        yyyymm = this_ym, var = names(coef(ols)), coef = ols$coef,
        r2 = var(ols$fitted.values) / var(ols$model$ret),
        r2_no_ofi = var(ols_no_ofi$fitted.values) / var(ols_no_ofi$model$ret),
        obs = nrow(data[yyyymm == this_ym])
      )
      return(out)
    }

    # first regress all together
    nc <- detectCores() - 2
    tmp <- rbindlist(mclapply(unique(data[, yyyymm]), p.get_one_period, mc.cores = nc))[var %in% paste0("ofi_bin", 1:3)]
    gc()

    cc <- as.matrix(tmp[, mean(coef), var][, V1])
    C <- cov(cast(tmp[, list(yyyymm, var, coef)], yyyymm ~ var, value = "coef")[, 2:4]) / nrow(tmp)

    b_12 <- matrix(c(-1, 1, 0))
    b_23 <- matrix(c(0, -1, 1))
    b_13 <- matrix(c(-1, 0, 1))

    out <- data.table(
      type = this_type,
      r2 = mean(tmp[, r2]),
      r2_no_ofi = mean(tmp[, r2_no_ofi]),
      obs = unique(tmp[, .(yyyymm, obs)])[, sum(obs)],
      var = c("M1", "M2", "M3", "M2-M1", "M3-M2", "M3-M1"),
      coef = c(
        cc,
        (t(b_12) %*% cc)[1],
        (t(b_23) %*% cc)[1],
        (t(b_13) %*% cc)[1]
      ),
      se = c(
        sqrt(diag(C)),
        sqrt((t(b_12) %*% C %*% b_12)[1]),
        sqrt((t(b_23) %*% C %*% b_23)[1]),
        sqrt((t(b_13) %*% C %*% b_13)[1])
      )
    )
    out[, tstat := coef / se]
    out[, idx := idx]
    return(out)
  }

  # ()
  idx <- 1
  ff <- copy(ff0_1)
  if (this_type == "BMI") {
    ff <- paste0(ff, "+", paste0(controls_bmi, collapse = "+"))
  }
  out <- p.get_one_fm(ff)

  # (char)
  idx <- 2
  ff <- copy(ff0_2)
  if (this_type == "BMI") {
    ff <- paste0(ff, "+", paste0(controls_bmi, collapse = "+"))
  }
  out <- rbind(out, p.get_one_fm(ff))

  # (char, liq)
  idx <- 3
  ff <- copy(ff0_3)
  if (this_type == "BMI") {
    ff <- paste0(ff, "+", paste0(controls_bmi, collapse = "+"))
  }
  out <- rbind(out, p.get_one_fm(ff))

  # then add interactions
  for (this_v in c(controls_char, controls_liq)) {
    # print(this_v)

    setnames(data, this_v, "xx")
    data[, yy := xx * ofi]
    setnames(data, c("xx", "yy"), c(this_v, paste0("ofi_", this_v)))
    ff <- paste0(ff, " + ofi_", this_v)
    idx <- idx + 1

    out <- rbind(out, p.get_one_fm(ff))
  }

  # name the variables being added
  tmp <- data.table(
    idx = 1:(length(controls_list) + 4),
    var_added = c("none_init", "controls_char", "controls_char+controls_liq", "none", controls_list),
    var_type = c(
      rep("", 4), rep("return-predicting chars", length(controls_char)),
      rep("liquidity", length(controls_liq))
    )
  )
  out <- merge(out, tmp, by = "idx", all.x = T)
  rm(tmp)
  out[, type := this_type]

  # also keep track of sd(OFI)
  out[, stdev_ofi := data[, .(ss = last(ss)), yyyymm][, mean(ss)]]
  out_all <- rbind(out_all, out)
  toc()
}

dir.create("tmp/price_impact/multiplier_by_shock_size_quarterly/", showWarnings = F, recursive = T)
saveRDS(out_all, "tmp/price_impact/multiplier_by_shock_size_quarterly/fm_by_stdev_based_bins.RDS")


# # === FM with pos/neg (not clear)
#
# # takes around X mins
# out_all = data.table()
# for (this_type in unique(data_all[, type])){
#   print(this_type)
#
#   tic()
#   data = copy(data_all[type == this_type])
#
#   # function to do Fama-MacBeth
#   p.get_one_fm = function(ff){
#
#     p.get_one_period = function(this_ym){
#       # ols = lm(paste0('ret ~ ofi + I(ofi_absofi * (ofi >= 0)) + I(ofi_absofi * (ofi < 0)) + ', ff), data[yyyymm == this_ym])
#       ols = lm(paste0('ret ~ ofi + ofi_absofi + I(ofi_absofi * (ofi < 0)) + ', ff), data[yyyymm == this_ym])
#       ols_no_ofi = lm(paste0('ret ~ ', ff), data[yyyymm == this_ym])
#       out = data.table(yyyymm = this_ym, var = names(coef(ols)), coef = ols$coef,
#                        r2 = var(ols$fitted.values)/var(ols$model$ret),
#                        r2_no_ofi = var(ols_no_ofi$fitted.values)/var(ols_no_ofi$model$ret))
#       return(out)
#     }
#
#     # first regress all together
#     nc = detectCores() - 2
#     out = rbindlist(mclapply(unique(data[, yyyymm]), p.get_one_period, mc.cores = nc)); gc()
#
#     tt = unique(out[, .(yyyymm, r2, r2_no_ofi)])
#     tt = tt[, .(r2 = mean(r2), r2_no_ofi = mean(r2_no_ofi))]
#
#     out = out[, list(coef = mean(coef, na.rm = T),
#                      se = sd(coef, na.rm = T)/sqrt(sum(!is.na(coef)))), var]
#     out[, r2 := tt[, r2]]
#     out[, r2_no_ofi := tt[, r2_no_ofi]]; rm(tt)
#     out[, obs := nrow(data)]
#     out[, idx := idx]
#     return(out)
#   }
#
#   # ()
#   idx = 1
#   ff = copy(ff0_1)
#   if (this_type == 'BMI'){
#     ff = paste0(ff, '+', paste0(controls_bmi, collapse = '+'))}
#   out = p.get_one_fm(ff)
#
#   # (char)
#   idx = 2
#   ff = copy(ff0_2)
#   if (this_type == 'BMI'){
#     ff = paste0(ff, '+', paste0(controls_bmi, collapse = '+'))}
#   out = rbind(out, p.get_one_fm(ff))
#
#   # (char, liq)
#   idx = 3
#   ff = copy(ff0_3)
#   if (this_type == 'BMI'){
#     ff = paste0(ff, '+', paste0(controls_bmi, collapse = '+'))}
#   out = rbind(out, p.get_one_fm(ff))
#
#   # then add interactions
#   for (this_v in c(controls_char, controls_liq)){
#     # print(this_v)
#
#     setnames(data, this_v, 'xx')
#     data[, yy := xx * ofi]
#     setnames(data, c('xx','yy'), c(this_v, paste0('ofi_', this_v)))
#     ff = paste0(ff, ' + ofi_', this_v)
#     idx = idx + 1
#
#     out = rbind(out, p.get_one_fm(ff))
#   }
#
#   # name the variables being added
#   tmp = data.table(idx = 1:(length(controls_list)+4),
#                    var_added = c('none_init','controls_char','controls_char+controls_liq','none',controls_list),
#                    var_type = c(rep('',4), rep('return-predicting chars', length(controls_char)),
#                                 rep('liquidity', length(controls_liq))))
#   out = merge(out, tmp, by = 'idx', all.x = T); rm(tmp)
#   out[, type := this_type]
#
#   out_all = rbind(out_all, out)
#   toc()
# }
# out_all[, tstat := coef/se]
#
# out_all = out_all[type != 'OFI']
# out_all[type == 'OFI_resid', type := 'OFI']
#
# # vv = c('ofi','I(ofi_absofi * (ofi >= 0))','I(ofi_absofi * (ofi < 0))')
# vv = c('ofi','ofi_absofi','I(ofi_absofi * (ofi < 0))')
# dcast(out_all[(var %in% vv) & (idx %in% 1:3)], var ~ paste0(type, '_', idx), value.var = 'coef')
# dcast(out_all[(var %in% vv) & (idx %in% 1:3)], var ~ paste0(type, '_', idx), value.var = 'tstat')

# === panel regression as robustness

out_all <- data.table()
for (this_type in unique(data_all[, type])) {
  print(this_type)

  tic()
  data <- copy(data_all[type == this_type])
  data[, ret := ret / 100]
  data[, ofi := ofi / 100]
  data <- data[0 == rowSums(is.na(data))]

  out <- data.table()

  # ()
  idx <- 1
  ff <- copy(ff0_1)
  if (this_type == "BMI") {
    ff <- paste0(ff, "+", paste0(controls_bmi, collapse = "+"))
  }
  ff <- paste0(ff, " | yyyymm")

  ols <- feols(as.formula(paste0("ret ~ ofi + ofi_absofi + ", ff)), data, cluster = c("yyyymm", "permno"))
  ols_no_ofi <- feols(as.formula(paste0("ret ~ ", ff)), data, cluster = "yyyymm")
  out <- rbind(out, data.table(
    spec_idx = idx, type = this_type, var = names(coef(ols)), coef = coef(ols),
    se = sqrt(diag(vcov(ols))), obs = nrow(data),
    r2 = r2(ols)["ar2"], r2_no_ofi = r2(ols_no_ofi)["ar2"]
  ))


  # (char)
  idx <- 2
  ff <- copy(ff0_2)
  if (this_type == "BMI") {
    ff <- paste0(ff, "+", paste0(controls_bmi, collapse = "+"))
  }
  ff <- paste0(ff, " | yyyymm")

  ols <- feols(as.formula(paste0("ret ~ ofi + ofi_absofi + ", ff)), data, cluster = c("yyyymm", "permno"))
  ols_no_ofi <- feols(as.formula(paste0("ret ~ ", ff)), data, cluster = "yyyymm")
  out <- rbind(out, data.table(
    spec_idx = idx, type = this_type, var = names(coef(ols)), coef = coef(ols),
    se = sqrt(diag(vcov(ols))), obs = nrow(data),
    r2 = r2(ols)["ar2"], r2_no_ofi = r2(ols_no_ofi)["ar2"]
  ))


  # (char, liq)
  idx <- 3
  ff <- copy(ff0_3)
  if (this_type == "BMI") {
    ff <- paste0(ff, "+", paste0(controls_bmi, collapse = "+"))
  }
  ff <- paste0(ff, " | yyyymm")

  ols <- feols(as.formula(paste0("ret ~ ofi + ofi_absofi + ", ff)), data, cluster = c("yyyymm", "permno"))
  ols_no_ofi <- feols(as.formula(paste0("ret ~ ", ff)), data, cluster = "yyyymm")
  out <- rbind(out, data.table(
    spec_idx = idx, type = this_type, var = names(coef(ols)), coef = coef(ols),
    se = sqrt(diag(vcov(ols))), obs = nrow(data),
    r2 = r2(ols)["ar2"], r2_no_ofi = r2(ols_no_ofi)["ar2"]
  ))


  # # then add interactions
  # for (this_v in c(controls_char, controls_liq)){
  #   # print(this_v)
  #
  #   setnames(data, this_v, 'xx')
  #   data[, yy := xx * ofi]
  #   setnames(data, c('xx','yy'), c(this_v, paste0('ofi_', this_v)))
  #   ff = paste0(ff, ' + ofi_', this_v)
  #   idx = idx + 1
  #
  #   out = rbind(out, p.get_one_fm(ff))
  # }
  #
  # # name the variables being added
  # tmp = data.table(idx = 1:(length(controls_list)+4),
  #                  var_added = c('none_init','controls_char','controls_char+controls_liq','none',controls_list),
  #                  var_type = c(rep('',4), rep('return-predicting chars', length(controls_char)),
  #                               rep('liquidity', length(controls_liq))))
  # out = merge(out, tmp, by = 'idx', all.x = T); rm(tmp)
  # out[, type := this_type]
  #
  out_all <- rbind(out_all, out)
  toc()
}

saveRDS(out_all, "tmp/price_impact/regression_contemp/full_sample_panel.RDS")

# === panel with std based-bins

out_all <- data.table()
for (this_type in unique(data_all[, type])) {
  print(this_type)

  tic()
  data <- copy(data_all[type == this_type])
  data[, ret := ret / 100]
  data[, ofi := ofi / 100]
  data <- data[0 == rowSums(is.na(data))]

  # bins
  data <- merge(data, data[, .(ss = sd(ofi)), yyyymm], by = "yyyymm")
  data[, bin := 1]
  data[abs(ofi) > ss, bin := 2]
  data[abs(ofi) > 2 * ss, bin := 3]
  data[, ofi_bin1 := ifelse(bin == 1, ofi, 0)]
  data[, ofi_bin2 := ifelse(bin == 2, ofi, 0)]
  data[, ofi_bin3 := ifelse(bin == 3, ofi, 0)]

  out <- data.table()

  # ()
  ff <- copy(ff0_1)
  if (this_type == "BMI") {
    ff <- paste0(ff, "+", paste0(controls_bmi, collapse = "+"))
  }
  ff <- paste0(ff, " | yyyymm")

  ols <- feols(as.formula(paste0("ret ~ ofi_bin1 + ofi_bin2 + ofi_bin3 + ", ff)), data, cluster = c("yyyymm", "permno"))
  ols_no_ofi <- feols(as.formula(paste0("ret ~ ", ff)), data, cluster = "yyyymm")

  # get coef differences
  cc <- as.matrix(coef(ols)[1:3])
  C <- vcov(ols)[1:3, 1:3]

  b_12 <- matrix(c(-1, 1, 0))
  b_23 <- matrix(c(0, -1, 1))
  b_13 <- matrix(c(-1, 0, 1))

  out <- rbind(out, data.table(
    spec_idx = 1, type = this_type,
    r2 = r2(ols)["ar2"],
    r2_no_ofi = r2(ols_no_ofi)["ar2"],
    obs = nrow(data),
    var = c("M1", "M2", "M3", "M2-M1", "M3-M2", "M3-M1"),
    coef = c(
      cc,
      (t(b_12) %*% cc)[1],
      (t(b_23) %*% cc)[1],
      (t(b_13) %*% cc)[1]
    ),
    se = c(
      sqrt(diag(C)),
      sqrt((t(b_12) %*% C %*% b_12)[1]),
      sqrt((t(b_23) %*% C %*% b_23)[1]),
      sqrt((t(b_13) %*% C %*% b_13)[1])
    )
  )[, tstat := coef / se])


  # (char)
  ff <- copy(ff0_2)
  if (this_type == "BMI") {
    ff <- paste0(ff, "+", paste0(controls_bmi, collapse = "+"))
  }
  ff <- paste0(ff, " | yyyymm")

  ols <- feols(as.formula(paste0("ret ~ ofi_bin1 + ofi_bin2 + ofi_bin3 + ", ff)), data, cluster = c("yyyymm", "permno"))
  ols_no_ofi <- feols(as.formula(paste0("ret ~ ", ff)), data, cluster = "yyyymm")

  # get coef differences
  cc <- as.matrix(coef(ols)[1:3])
  C <- vcov(ols)[1:3, 1:3]

  b_12 <- matrix(c(-1, 1, 0))
  b_23 <- matrix(c(0, -1, 1))
  b_13 <- matrix(c(-1, 0, 1))

  out <- rbind(out, data.table(
    spec_idx = 2, type = this_type,
    r2 = r2(ols)["ar2"],
    r2_no_ofi = r2(ols_no_ofi)["ar2"],
    obs = nrow(data),
    var = c("M1", "M2", "M3", "M2-M1", "M3-M2", "M3-M1"),
    coef = c(
      cc,
      (t(b_12) %*% cc)[1],
      (t(b_23) %*% cc)[1],
      (t(b_13) %*% cc)[1]
    ),
    se = c(
      sqrt(diag(C)),
      sqrt((t(b_12) %*% C %*% b_12)[1]),
      sqrt((t(b_23) %*% C %*% b_23)[1]),
      sqrt((t(b_13) %*% C %*% b_13)[1])
    )
  )[, tstat := coef / se])


  # (char, liq)
  ff <- copy(ff0_3)
  if (this_type == "BMI") {
    ff <- paste0(ff, "+", paste0(controls_bmi, collapse = "+"))
  }
  ff <- paste0(ff, " | yyyymm")

  ols <- feols(as.formula(paste0("ret ~ ofi_bin1 + ofi_bin2 + ofi_bin3 + ", ff)), data, cluster = c("yyyymm", "permno"))
  ols_no_ofi <- feols(as.formula(paste0("ret ~ ", ff)), data, cluster = "yyyymm")

  # get coef differences
  cc <- as.matrix(coef(ols)[1:3])
  C <- vcov(ols)[1:3, 1:3]

  b_12 <- matrix(c(-1, 1, 0))
  b_23 <- matrix(c(0, -1, 1))
  b_13 <- matrix(c(-1, 0, 1))

  out <- rbind(out, data.table(
    spec_idx = 3, type = this_type,
    r2 = r2(ols)["ar2"],
    r2_no_ofi = r2(ols_no_ofi)["ar2"],
    obs = nrow(data),
    var = c("M1", "M2", "M3", "M2-M1", "M3-M2", "M3-M1"),
    coef = c(
      cc,
      (t(b_12) %*% cc)[1],
      (t(b_23) %*% cc)[1],
      (t(b_13) %*% cc)[1]
    ),
    se = c(
      sqrt(diag(C)),
      sqrt((t(b_12) %*% C %*% b_12)[1]),
      sqrt((t(b_23) %*% C %*% b_23)[1]),
      sqrt((t(b_13) %*% C %*% b_13)[1])
    )
  )[, tstat := coef / se])


  # # then add interactions
  # for (this_v in c(controls_char, controls_liq)){
  #   # print(this_v)
  #
  #   setnames(data, this_v, 'xx')
  #   data[, yy := xx * ofi]
  #   setnames(data, c('xx','yy'), c(this_v, paste0('ofi_', this_v)))
  #   ff = paste0(ff, ' + ofi_', this_v)
  #   idx = idx + 1
  #
  #   out = rbind(out, p.get_one_fm(ff))
  # }
  #
  # # name the variables being added
  # tmp = data.table(idx = 1:(length(controls_list)+4),
  #                  var_added = c('none_init','controls_char','controls_char+controls_liq','none',controls_list),
  #                  var_type = c(rep('',4), rep('return-predicting chars', length(controls_char)),
  #                               rep('liquidity', length(controls_liq))))
  # out = merge(out, tmp, by = 'idx', all.x = T); rm(tmp)
  # out[, type := this_type]
  #
  out_all <- rbind(out_all, out)
  toc()
}

saveRDS(out_all, "tmp/price_impact/regression_contemp/full_sample_panel_with_stdev_bins.RDS")
