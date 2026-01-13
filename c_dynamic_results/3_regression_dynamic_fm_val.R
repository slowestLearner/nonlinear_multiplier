# Contemp + dynamic regression, inspired by Val's question
# code modified from code/c_dynamic_results/2_regression_dynamic_fm.R

rm(list = ls())
source("../utilities/runmefirst.R")
source("../utilities/regressions.R")

# load regression data
tic("Loading regression data")
data <- readRDS("../tmp/raw_data/reg_inputs/reg_table_dynamic.RDS")

# separate controls and main variables
vars_id <- c("yyyymm", "permno", "type")
vars_reg <- c("ret", "ofi", "cumofi_1", "cumofi_2", "cumofi_3", "cumofi_4")
vars_controls <- setdiff(names(data), c(vars_id, vars_reg))

# put regression-related data into long format
data_controls <- data[, c(vars_id, vars_controls), with = F]
setnames(data_controls, "type", "demand_type")

data_reg <- data[, c(vars_id, vars_reg), with = F]
data_reg <- melt(data_reg, id.vars = c("yyyymm", "permno", "type", "ret", "ofi"), variable.name = "hor", value.name = "cumofi_lag")
data_reg <- data_reg[0 == rowSums(is.na(data_reg))]
data_reg[, hor := as.integer(gsub("cumofi_", "", hor))]
rm(data, vars_id, vars_reg, vars_controls)

# NEW: also have contemp bins
data_reg[, demand_type := type]
data_reg[, type := paste0(demand_type, "_", hor, "lag")]
data_reg[, sd_ofi := sd(cumofi_lag), .(yyyymm, type)]
data_reg[, bin := 1]
data_reg[abs(cumofi_lag) > sd_ofi, bin := 2]
data_reg[abs(cumofi_lag) > 2 * sd_ofi, bin := 3]
data_reg[, ofi_bin1 := ifelse(bin == 1, ofi, 0)]
data_reg[, ofi_bin2 := ifelse(bin == 2, ofi, 0)]
data_reg[, ofi_bin3 := ifelse(bin == 3, ofi, 0)]
data_reg[, sd_ofi := NULL]

data_reg[, sd_ofi := sd(ofi), .(yyyymm, type)]
data_reg[, cbin := 1]
data_reg[abs(ofi) > sd_ofi, cbin := 2]
data_reg[abs(ofi) > 2 * sd_ofi, cbin := 3]
data_reg[, ofi_cbin1 := ifelse(cbin == 1, ofi, 0)]
data_reg[, ofi_cbin2 := ifelse(cbin == 2, ofi, 0)]
data_reg[, ofi_cbin3 := ifelse(cbin == 3, ofi, 0)]
data_reg[, sd_ofi := NULL]
# data_reg[, ofi_absofi := ofi * abs(cumofi_lag)]

# get control variable names
cdata <- readRDS("../tmp/raw_data/controls/controls_classification.RDS")
controls_char <- cdata[control_type == "return-predictor", var]
controls_liq <- cdata[control_type == "liquidity", var]
controls_list <- c(controls_char, controls_liq)
rm(cdata)
toc()

# ###########################################################
# Fama-MacBeth
# ###########################################################

# control variables for different specifications
control_formulas <- c(
  "1",
  paste0(c(controls_char), collapse = "+"),
  paste0(c(controls_char, controls_liq), collapse = "+")
)

# -- process one type of data. reg_spec = "nonlinear" or "stdev"
# This version omits bin3 (lagged bin 3)
p.process_one_type_no_bin3 <- function(data, reg_spec = "nonlinear") {
  data <- merge(data, data_controls, by = c("yyyymm", "permno", "demand_type"), all.x = T)
  data <- data[0 == rowSums(is.na(data))]
  this_type <- data[1, type]
  
  out_all <- data.table() # save all results here
  
  # introduce direct controls
  for (spec_idx in 1:length(control_formulas)) {
    
    # put together formula
    if (reg_spec == "nonlinear") {
      ff <- paste0("ret ~ ", control_formulas[spec_idx], " + ofi + ofi_absofi")
    } else if (reg_spec == "stdev") {
      ff <- paste0("ret ~ ", control_formulas[spec_idx], " + ofi_cbin1 + ofi_cbin2 + ofi_cbin3 + ofi_bin1 + ofi_bin2")
    }
    if (this_type == "BMI") {
      ff <- paste0(ff, "+", paste0(controls_bmi, collapse = "+"))
    }
    
    # helper function for FM
    p.fama_macbeth <- function(data, ff, compare_coefs = FALSE) {
      # regression for one period
      p.get_one_period <- function(this_ym) {
        ols <- lm(ff, data[yyyymm == this_ym])
        # Get the dependent variable name dynamically
        dep_var_name <- all.vars(as.formula(ff))[1]
        return(data.table(
          yyyymm = this_ym, var = names(coef(ols)), coef = ols$coef,
          r2 = var(ols$fitted.values) / var(ols$model[[dep_var_name]])
        ))
      }
      
      # regression by period
      out <- rbindlist(lapply(unique(data[, yyyymm]), p.get_one_period))
      
      # extract R2
      r2_data <- unique(out[, .(yyyymm, r2)])[, .(r2 = mean(r2))]
      out[, r2 := NULL]
      
      # let's do Newey-West
      yms <- sort(unique(out[, yyyymm]))
      
      # compare coefficient differences
      if (compare_coefs) {
        out <- rbind(out, data.table(
          yyyymm = yms, var = "ofi_cbin2 - ofi_cbin1",
          coef = out[var == "ofi_cbin2", coef] - out[var == "ofi_cbin1", coef]
        ))
        out <- rbind(out, data.table(
          yyyymm = yms, var = "ofi_cbin3 - ofi_cbin1",
          coef = out[var == "ofi_cbin3", coef] - out[var == "ofi_cbin1", coef]
        ))
        out <- rbind(out, data.table(
          yyyymm = yms, var = "ofi_cbin3 - ofi_cbin2",
          coef = out[var == "ofi_cbin3", coef] - out[var == "ofi_cbin2", coef]
        ))
      }
      
      # newey-west lag
      if ("hor" %in% names(data)) {
        this_hor <- data[1, hor]
      } else {
        this_hor <- 1
      }
      
      coef_data <- data.table()
      for (this_var in unique(out[, var])) {
        mm <- lm(coef ~ 1, out[var == this_var])
        coef_data <- rbind(coef_data, data.table(
          var = this_var, coef = mm$coef[1],
          se = sqrt(vcov(mm)[1, 1]),
          se_nw = sqrt(NeweyWest(mm, this_hor)[1, 1])
        ))
      }
      
      coef_data[, r2 := r2_data[, r2]]
      coef_data[, obs := nrow(data)]
      coef_data[, type := data[1, type]]
      coef_data[, nw_lag := this_hor]
      return(coef_data)
    }
    
    out <- p.fama_macbeth(data, ff, compare_coefs = ifelse(reg_spec == "nonlinear", FALSE, TRUE))
    # out <- p.fama_macbeth(data, ff)
    out[, spec_idx := spec_idx]
    out_all <- rbind(out_all, out)
    
    
    # then add interactions with ofi
    # for (this_v in c(controls_char, controls_liq)) {
    #   setnames(data, this_v, "xx")
    #   data[, yy := xx * ofi]
    #   setnames(data, c("xx", "yy"), c(this_v, paste0("ofi_", this_v)))
    #   ff <- paste0(ff, " + ofi_", this_v)
    #   spec_idx <- spec_idx + 1
    #
    #   out <- p.fama_macbeth(data, ff, compare_coefs = ifelse(reg_spec == "nonlinear", FALSE, TRUE))
    #   out[, spec_idx := spec_idx]
    #   out_all <- rbind(out_all, out)
    # }
    #
    # # name the control variables being added
    # tmp <- data.table(
    #   spec_idx = 1:(length(controls_list) + 3),
    #   var_added = c("none_init", "controls_char", "controls_char+controls_liq", controls_list),
    #   var_type = c(
    #     rep("", 3), rep("return-predicting chars", length(controls_char)),
    #     rep("liquidity", length(controls_liq))
    #   )
    # )
    # out_all <- merge(out_all, tmp, by = "spec_idx", all.x = T)
    out_all[, type := this_type]
    
    
  }
  return(out_all)
}


# stdev-based specification---
tic("dynamic fm: stdev - omit bin3")
out <- rbindlist(mclapply(split(data_reg, by = "type"), function(x) {
  p.process_one_type_no_bin3(x, reg_spec = "stdev")
}, mc.cores = nc))
to_dir = '../tmp/price_impact/regression_contemp_and_dynamic/'
dir.create(to_dir, recursive = T)
saveRDS(out, paste0(to_dir, 'fm_stdev_missing_bin3.RDS'))
toc()


p.process_one_type_no_cbin3 <- function(data, reg_spec = "nonlinear") {
  data <- merge(data, data_controls, by = c("yyyymm", "permno", "demand_type"), all.x = T)
  data <- data[0 == rowSums(is.na(data))]
  this_type <- data[1, type]
  
  out_all <- data.table() # save all results here
  
  # introduce direct controls
  for (spec_idx in 1:length(control_formulas)) {
    
    # put together formula
    if (reg_spec == "nonlinear") {
      ff <- paste0("ret ~ ", control_formulas[spec_idx], " + ofi + ofi_absofi")
    } else if (reg_spec == "stdev") {
      ff <- paste0("ret ~ ", control_formulas[spec_idx], " + ofi_cbin1 + ofi_cbin2 + ofi_bin1 + ofi_bin2 + ofi_bin3")
    }
    if (this_type == "BMI") {
      ff <- paste0(ff, "+", paste0(controls_bmi, collapse = "+"))
    }
    
    # helper function for FM
    p.fama_macbeth <- function(data, ff, compare_coefs = FALSE) {
      # regression for one period
      p.get_one_period <- function(this_ym) {
        ols <- lm(ff, data[yyyymm == this_ym])
        # Get the dependent variable name dynamically
        dep_var_name <- all.vars(as.formula(ff))[1]
        return(data.table(
          yyyymm = this_ym, var = names(coef(ols)), coef = ols$coef,
          r2 = var(ols$fitted.values) / var(ols$model[[dep_var_name]])
        ))
      }
      
      # regression by period
      out <- rbindlist(lapply(unique(data[, yyyymm]), p.get_one_period))
      
      # extract R2
      r2_data <- unique(out[, .(yyyymm, r2)])[, .(r2 = mean(r2))]
      out[, r2 := NULL]
      
      # let's do Newey-West
      yms <- sort(unique(out[, yyyymm]))
      
      # compare coefficient differences
      if (compare_coefs) {
        out <- rbind(out, data.table(
          yyyymm = yms, var = "ofi_bin2 - ofi_bin1",
          coef = out[var == "ofi_bin2", coef] - out[var == "ofi_bin1", coef]
        ))
        out <- rbind(out, data.table(
          yyyymm = yms, var = "ofi_bin3 - ofi_bin1",
          coef = out[var == "ofi_bin3", coef] - out[var == "ofi_bin1", coef]
        ))
        out <- rbind(out, data.table(
          yyyymm = yms, var = "ofi_bin3 - ofi_bin2",
          coef = out[var == "ofi_bin3", coef] - out[var == "ofi_bin2", coef]
        ))
      }
      
      # newey-west lag
      if ("hor" %in% names(data)) {
        this_hor <- data[1, hor]
      } else {
        this_hor <- 1
      }
      
      coef_data <- data.table()
      for (this_var in unique(out[, var])) {
        mm <- lm(coef ~ 1, out[var == this_var])
        coef_data <- rbind(coef_data, data.table(
          var = this_var, coef = mm$coef[1],
          se = sqrt(vcov(mm)[1, 1]),
          se_nw = sqrt(NeweyWest(mm, this_hor)[1, 1])
        ))
      }
      
      coef_data[, r2 := r2_data[, r2]]
      coef_data[, obs := nrow(data)]
      coef_data[, type := data[1, type]]
      coef_data[, nw_lag := this_hor]
      return(coef_data)
    }
    
    out <- p.fama_macbeth(data, ff, compare_coefs = ifelse(reg_spec == "nonlinear", FALSE, TRUE))
    out[, spec_idx := spec_idx]
    out_all <- rbind(out_all, out)
    
    
    out_all[, type := this_type]
    
    
  }
  return(out_all)
}


# stdev-based specification---
tic("dynamic fm: stdev - omit cbin3")
out <- rbindlist(mclapply(split(data_reg, by = "type"), function(x) {
  p.process_one_type_no_cbin3(x, reg_spec = "stdev")
}, mc.cores = nc))
to_dir = '../tmp/price_impact/regression_contemp_and_dynamic/'
dir.create(to_dir, recursive = T)
saveRDS(out, paste0(to_dir, 'fm_stdev_missing_cbin3.RDS'))
toc()

# --- taking a look

# out <- readRDS('../tmp/price_impact/regression_contemp_and_dynamic/fm_stdev_missing_bin3.RDS')
out <- readRDS('../tmp/price_impact/regression_contemp_and_dynamic/fm_stdev_missing_cbin3.RDS')

# take a look
tmp <- copy(out[spec_idx == 3])
tmp <- tmp[grepl("bin", var)]
tmp[, tstat := round(coef / se, 2)]

# lagged
dcast(tmp[var %in% paste0('ofi_bin', 1:3)], var ~ type, value.var = "coef")
dcast(tmp[var %in% paste0('ofi_bin', 1:3)], var ~ type, value.var = "tstat")

# current
dcast(tmp[var %in% paste0('ofi_cbin', 1:3)], var ~ type, value.var = "coef")
dcast(tmp[var %in% paste0('ofi_cbin', 1:3)], var ~ type, value.var = "tstat")

dcast(tmp[grepl(' - ', var, fixed = T)], var ~ type, value.var = "coef")
dcast(tmp[grepl(' - ', var, fixed = T)], var ~ type, value.var = "tstat")

