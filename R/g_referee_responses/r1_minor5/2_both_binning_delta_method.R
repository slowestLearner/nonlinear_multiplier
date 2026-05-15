# --- Binning by cross-section vs full-sample
library(this.path)
setwd(this.path::this.dir())
source("../../../../formal_tests/code/R/utilities/runmefirst.R")
library(numDeriv)
options(width = 150)

# load regression data
tic("preparing data")
data_all <- readRDS("../../../../formal_tests/code/R/tmp/raw_data/reg_inputs/reg_table_static.RDS")[type != "OFI_pre_whitened"]

# get control variable names
cdata <- readRDS("../../../../formal_tests/code/R/tmp/raw_data/controls/controls_classification.RDS")
controls_char <- cdata[control_type == "return-predictor", var]
controls_liq <- cdata[control_type == "liquidity", var]
controls_list <- c(controls_char, controls_liq)
rm(cdata)

# get names for bmi controls
tmp <- readRDS("../../../../formal_tests/code/R/tmp/raw_data/controls/controls_for_BMI.RDS")
controls_bmi <- setdiff(names(tmp), c("yyyymm", "permno"))
rm(tmp)
toc()

# add full-sample bins
data_all[, sd_ofi_full := sd(ofi), by = .(type)]
data_all[, bin_full := ifelse(abs(ofi) > 2 * sd_ofi_full, 3, ifelse(abs(ofi) > sd_ofi_full, 2, 1))]
data_all[, ofi_bin1_full := ifelse(bin_full == 1, ofi, 0)]
data_all[, ofi_bin2_full := ifelse(bin_full == 2, ofi, 0)]
data_all[, ofi_bin3_full := ifelse(bin_full == 3, ofi, 0)]

# ###########################################################
# Fama-MacBeth
# ###########################################################

# NEW fama-macbeth function that returns the VCOV of the AVERAGED coefficients
p.fama_macbeth <- function(data, ff, compare_coefs = FALSE, return_list = FALSE) {
    library(sandwich)
    library(data.table)

    # 1. Regression for one period
    p.get_one_period <- function(this_ym) {
        # Subset data for the period
        sub_data <- data[yyyymm == this_ym]
        if (nrow(sub_data) == 0) {
            return(NULL)
        }

        ols <- lm(ff, data = sub_data)
        dep_var_name <- all.vars(as.formula(ff))[1]

        return(data.table(
            yyyymm = this_ym,
            var = names(coef(ols)),
            coef = as.numeric(coef(ols)),
            r2 = var(ols$fitted.values) / var(sub_data[[dep_var_name]])
        ))
    }

    # 2. Run regressions by period
    out <- rbindlist(lapply(unique(data[, yyyymm]), p.get_one_period))

    # Extract Mean R2
    r2_mean <- mean(unique(out[, .(yyyymm, r2)])$r2, na.rm = TRUE)
    out[, r2 := NULL]

    # 3. Add coefficient differences if requested
    if (compare_coefs) {
        yms <- sort(unique(out[, yyyymm]))
        diff_vars <- list(
            c("ofi_bin2 - ofi_bin1", "ofi_bin2", "ofi_bin1"),
            c("ofi_bin3 - ofi_bin1", "ofi_bin3", "ofi_bin1"),
            c("ofi_bin3 - ofi_bin2", "ofi_bin3", "ofi_bin2")
        )

        for (dv in diff_vars) {
            dt_diff <- out[var %in% dv[2:3]]
            dt_diff <- dcast(dt_diff, yyyymm ~ var, value.var = "coef")
            dt_diff[, coef := get(dv[2]) - get(dv[3])]
            dt_diff[, var := dv[1]]
            out <- rbind(out, dt_diff[, .(yyyymm, var, coef)])
        }
    }

    # 4. Reshape to Wide format for VCV calculation (T x K)
    wide_out <- dcast(out, yyyymm ~ var, value.var = "coef")
    var_names <- setdiff(names(wide_out), "yyyymm")

    # Lag for Newey-West
    this_hor <- if ("hor" %in% names(data)) data[1, hor] else 1

    # 5. Multivariate intercept-only regression to get Full VCV
    # We use as.matrix to run an 'mlm' (multivariate linear model)
    fit_all <- lm(as.matrix(wide_out[, ..var_names]) ~ 1)

    # Full Newey-West VCV matrix
    # Note: vcovHAC/NeweyWest on mlm returns a KxK matrix for the intercepts
    full_vcv <- NeweyWest(fit_all, lag = this_hor, prewhite = FALSE)

    # The matrix comes out with names like "Intercept:var_name". Let's clean them.
    colnames(full_vcv) <- gsub(":\\(Intercept\\)|\\(Intercept\\):", "", colnames(full_vcv))
    rownames(full_vcv) <- colnames(full_vcv)

    # 6. Build the summary table
    coef_data <- data.table(
        var = var_names,
        coef = colMeans(wide_out[, ..var_names], na.rm = TRUE),
        se_nw = sqrt(diag(full_vcv))
    )

    coef_data[, `:=`(
        r2 = r2_mean,
        obs = nrow(data),
        type = data[1, type],
        nw_lag = this_hor
    )]

    if (return_list) {
        return(list(table = coef_data, vcv = full_vcv))
    } else {
        return(coef_data)
    }
}

# control variables for different regression specifications
controls <- paste0(c(controls_char, controls_liq), collapse = "+")

# worker function to estimate regression with one type of demand variable. reg_spec = "nonlinear" or "stdev"
# data <- data_list[[1]]
p.process_one_type <- function(data) {
    this_type <- data[1, type]

    # 1. Estimate Regression
    # Note: Ensure p.fama_macbeth returns the VCOV of the AVERAGED coefficients
    # and the mean coefficients.
    ff <- paste0("ret ~ ", controls, " + ofi_bin1 + ofi_bin2 + ofi_bin3 + ofi_bin2_full + ofi_bin3_full")
    if (this_type == "BMI") {
        ff <- paste0(ff, "+", paste0(controls_bmi, collapse = "+"))
    }

    # Assuming p.fama_macbeth returns a list with $coef and $vcov
    fm_results <- p.fama_macbeth(data, ff, return_list = TRUE)
    fm_results$table

    # Isolate the coefficients for the 5 ofi bins
    # Indices might change depending on how your 'controls' are ordered
    target_vars <- c("ofi_bin1", "ofi_bin2", "ofi_bin3", "ofi_bin2_full", "ofi_bin3_full")
    theta_hat <- fm_results$table[var %in% target_vars, coef]
    sigma_hat <- fm_results$vcv[target_vars, target_vars]

    # 2. Define the Share Function for the Delta Method
    # This function must only take the coefficient vector as the first argument
    calc_share_ts <- function(theta, df) {
        # Extract params
        b1 <- theta[1]
        b2 <- theta[2]
        b3 <- theta[3]
        g2 <- theta[4]
        g3 <- theta[5]

        # Construct components based on fixed indicators in the data
        # M_bin uses b1, b2, b3 based on the 'bin' column
        m_bin <- ifelse(df$bin == 1, b1, ifelse(df$bin == 2, b2, b3))

        # M_bin_full uses g1=0, g2, g3 based on the 'bin_full' column
        m_bin_full <- ifelse(df$bin_full == 1, 0, ifelse(df$bin_full == 2, g2, g3))

        m_imputed <- m_bin + m_bin_full

        # Calculate Share = Cov(M_bin, M_imputed) / Var(M_imputed)
        share <- cov(m_bin, m_imputed) / var(m_imputed)
        return(share)
    }

    # 3. Calculate Point Estimate
    share_point <- calc_share_ts(theta_hat, data)

    # 4. Calculate Gradient using numDeriv
    # We treat 'data' as fixed/constant for the differentiation
    grad_s <- grad(func = calc_share_ts, x = theta_hat, df = data)

    # 5. Delta Method SE: sqrt( grad' * Sigma * grad )
    se_share <- sqrt(as.numeric(t(grad_s) %*% sigma_hat %*% grad_s))

    # 1. Re-extract the point estimates for the final table construction
    # Using named indexing to be 100% safe against R's alphabetical sorting
    # b1 <- theta_hat[1]
    # b2 <- theta_hat[2]
    # b3 <- theta_hat[3]
    # g2 <- theta_hat[4]
    # g3 <- theta_hat[5]

    b1 <- fm_results$table[var == "ofi_bin1", coef]
    b2 <- fm_results$table[var == "ofi_bin2", coef]
    b3 <- fm_results$table[var == "ofi_bin3", coef]
    g2 <- fm_results$table[var == "ofi_bin2_full", coef]
    g3 <- fm_results$table[var == "ofi_bin3_full", coef]

    # 2. Re-create the component series on the actual data
    m_ts <- ifelse(data$bin == 1, b1, ifelse(data$bin == 2, b2, b3))
    m_xs <- ifelse(data$bin_full == 1, 0, ifelse(data$bin_full == 2, g2, g3))
    m_total <- m_ts + m_xs

    # 3. Calculate Raw Covariances and Variance
    cov_ts <- cov(m_ts, m_total)
    cov_xs <- cov(m_xs, m_total)
    var_m <- var(m_total)

    M_decomp <- data.table(
        type           = this_type,
        # Raw Covariances
        cov_M_ts       = cov_ts,
        cov_M_xs       = cov_xs,
        var_M_total    = var_m,
        # Shares
        share_ts       = share_point, # From your calc_share_ts call
        share_xs       = 1 - share_point, # Mathematically guaranteed
        # Inference (SE is identical for both shares)
        share_se       = se_share,
        tstat_ts       = share_point / se_share,
        tstat_xs       = (1 - share_point) / se_share
    )

    # M_decomp <- data.table(
    #     type = this_type,
    #     share_ts = share_point,
    #     share_ts_se = se_share,
    #     share_ts_tstat = share_point / se_share,
    #     cov_M_cx = cov(
    #         ifelse(data$bin == 1, theta_hat[1], ifelse(data$bin == 2, theta_hat[2], theta_hat[3])),
    #         (ifelse(data$bin == 1, theta_hat[1], ifelse(data$bin == 2, theta_hat[2], theta_hat[3])) +
    #             ifelse(data$bin_full == 1, 0, ifelse(data$bin_full == 2, theta_hat[4], theta_hat[5])))
    #     )
    # )

    return(list(multipliers = fm_results$table[grepl("ofi_", var)], M_decomp = M_decomp))
}

tic("regressions")
data_list <- split(data_all, by = "type")
results_list <- mclapply(split(data_all, by = "type"), p.process_one_type, mc.cores = detectCores() - 2)
multipliers <- rbindlist(lapply(results_list, `[[`, "multipliers"))
M_decomp <- rbindlist(lapply(results_list, `[[`, "M_decomp"))
toc()

to_dir <- "tmp/binned_spec/"
dir.create(to_dir, recursive = T, showWarnings = F)
saveRDS(multipliers, paste0(to_dir, "multipliers.RDS"))
saveRDS(M_decomp, paste0(to_dir, "M_decomp.RDS"))

multipliers
M_decomp
