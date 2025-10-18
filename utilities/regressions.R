# function to estimate Fama-MacBeth regression
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
    yms <- out[var == first(var), yyyymm]

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
            se = sqrt(vcov(mm))[1, 1],
            se_nw = sqrt(NeweyWest(mm, this_hor)[1, 1])
        ))
    }
    coef_data[, r2 := r2_data[, r2]]
    coef_data[, obs := nrow(data)]
    coef_data[, type := data[1, type]]
    coef_data[, nw_lag := this_hor]
    return(coef_data)
}
