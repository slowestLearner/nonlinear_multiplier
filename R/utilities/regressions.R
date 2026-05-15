# --- function to estimate Fama-MacBeth regression
#
# the code loops through various specifications of the regression formula
# and reports the coefficients and standard errors. Reports both conventiional
# and NW standard errors. For the latter, if data[1, hor] is present, it will use that
# as the lag length; otherwise assumed to be 1.
#
# - inputs
# data: the data frame containing the data
# ff: the regression formula
# compare_coefs: whether to compare the coefficients of the different specifications
#
# - outputs
# returns: a data frame containing the coefficients and standard errors
#
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
  
  if (compare_coefs) {
    # 1. Identify all variables containing 'ofi_'
    all_vars <- unique(out$var)
    target_vars <- sort(grep("ofi_", all_vars, value = TRUE))
    
    # 2. Generate all permutations (v_j, v_i)
    # expand.grid creates every possible combination of the two vectors
    pairs_grid <- expand.grid(var_j = target_vars, var_i = target_vars, stringsAsFactors = FALSE)
    
    # 3. Filter out self-comparisons (where var_j == var_i)
    pairs_grid <- pairs_grid[pairs_grid$var_j != pairs_grid$var_i, ]
    
    # 4. Use lapply to iterate through the rows of the grid
    diff_list <- lapply(1:nrow(pairs_grid), function(idx) {
      var_j <- pairs_grid$var_j[idx]
      var_i <- pairs_grid$var_i[idx]
      
      # Extract coefficients for each variable
      dt_j <- out[var == var_j, .(yyyymm, coef_j = coef)]
      dt_i <- out[var == var_i, .(yyyymm, coef_i = coef)]
      
      # Merge by yyyymm to ensure temporal alignment
      res <- merge(dt_j, dt_i, by = "yyyymm")
      
      # Calculate the directional difference: v_j - v_i
      res[, coef := coef_j - coef_i]
      res[, var := paste0(var_j, " - ", var_i)]
      
      return(res[, .(yyyymm, var, coef)])
    })
    
    # Optional: combine into a single data.table
    out_diffs <- rbindlist(diff_list)
    out <- rbind(out, out_diffs)
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

# utility function: panel regression
p.panel_regression <- function(data, ff, compare_coefs = FALSE) {
    ols <- feols(as.formula(paste0(ff, " | yyyymm")), data, cluster = c("yyyymm", "permno"))
    out <- data.table(
        var = names(coef(ols)),
        coef = coef(ols),
        se = sqrt(diag(vcov(ols))),
        obs = ols$nobs, r2 = r2(ols)["ar2"]
    )

    # also report coef differences
    if (compare_coefs == TRUE) {
        var_indices <- names(coef(ols)) %in% paste0("ofi_bin", 1:3)

        cc <- matrix(coef(ols)[var_indices])
        C <- vcov(ols)[var_indices, var_indices]

        b_12 <- matrix(c(-1, 1, 0))
        b_23 <- matrix(c(0, -1, 1))
        b_13 <- matrix(c(-1, 0, 1))

        out <- rbind(out, data.table(
            var = c("ofi_bin2 - ofi_bin1", "ofi_bin3 - ofi_bin2", "ofi_bin3 - ofi_bin1"),
            coef = c(
                (t(b_12) %*% cc)[1],
                (t(b_23) %*% cc)[1],
                (t(b_13) %*% cc)[1]
            ),
            se = c(
                sqrt((t(b_12) %*% C %*% b_12)[1]),
                sqrt((t(b_23) %*% C %*% b_23)[1]),
                sqrt((t(b_13) %*% C %*% b_13)[1])
            ),
            obs = ols$nobs, r2 = r2(ols)["ar2"]
        ))
    }

    return(out)
}
