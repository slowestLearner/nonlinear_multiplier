# --- Do binning vs magnitude-based controls
library(this.path)
setwd(this.path::this.dir())
source("../../utilities/runmefirst.R")
library(numDeriv)
options(width = 120)

# load data
tic("preparing data")
data_all <- readRDS("../../tmp/raw_data/reg_inputs/reg_table_static.RDS")[type != "OFI_pre_whitened"]

# get control variable names
cdata <- readRDS("../../tmp/raw_data/controls/controls_classification.RDS")
controls_char <- cdata[control_type == "return-predictor", var]
controls_liq <- cdata[control_type == "liquidity", var]
controls_list <- c(controls_char, controls_liq)
rm(cdata)

# get names for bmi controls
tmp <- readRDS("../../tmp/raw_data/controls/controls_for_BMI.RDS")
controls_bmi <- setdiff(names(tmp), c("yyyymm", "permno"))
rm(tmp)
toc()

# do some size interactions
data_all[, abs_ofi := abs(ofi)]
data_all[, abs_ofi_x_ofi := abs_ofi * ofi]

# ###########################################################
# Fama-MacBeth
# ###########################################################

# new fama-macbeth function that returns the VCOV of the AVERAGED coefficients
p.fama_macbeth <- function(data, ff, return_list = FALSE) {
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

    # 4. Reshape to Wide format for VCV calculation (T x K)
    wide_out <- dcast(out, yyyymm ~ var, value.var = "coef") %>% setDT()
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
        se = sqrt(diag(vcov(fit_all))),
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

# worker function for each type
# data <- data_list[[1]]
p.process_one_type <- function(data) {
    this_type <- data[1, type]

    # control variables for different regression specifications
    specs <- list(
        "1",
        paste0(controls_char, collapse = "+"),
        paste0(c(controls_char, controls_liq), collapse = "+")
    )
    if (this_type == "BMI") {
        specs <- paste0(specs, "+", paste0(controls_bmi, collapse = "+"))
    }

    # We append each spec to lists defined outside the loop
    decomp_list <- list()
    mult_list <- list()

    for (spec_idx in 1:length(specs)) {
        # spec_idx <- 1

        # 1. Regression Specification
        ff <- paste0("ret ~ ofi_bin1 + ofi_bin2 + ofi_bin3 + abs_ofi_x_ofi + ", specs[[spec_idx]])

        # 2. Estimate via Fama-MacBeth
        # Ensure return_list = TRUE to get the VCV matrix
        fm_results <- p.fama_macbeth(data, ff, return_list = TRUE)

        # 3. Extract parameters using NAMED indexing
        # This bypasses any alphabetical sorting issues from dcast
        target_vars <- fm_results$table$var[grepl("ofi", fm_results$table$var)]
        theta_table <- fm_results$table[var %in% target_vars]

        # Create a named vector: c("ofi" = 0.1, "ofi_rank_inter" = 0.05, ...)
        theta_hat <- setNames(theta_table$coef, theta_table$var)
        sigma_hat <- fm_results$vcv[target_vars, target_vars]

        # 4. Define the Share Function for the Delta Method
        # We calculate the share of variance attributable to the Bin components vs Polynomial components
        # theta_vec <- theta_hat
        # df <- data
        calc_share_bin <- function(theta_vec, df) {
            # Safely extract coefficients (they might not exist in lower-degree specs, default to 0)
            b_bin1 <- if ("ofi_bin1" %in% names(theta_vec)) theta_vec["ofi_bin1"] else 0
            b_bin2 <- if ("ofi_bin2" %in% names(theta_vec)) theta_vec["ofi_bin2"] else 0
            b_bin3 <- if ("ofi_bin3" %in% names(theta_vec)) theta_vec["ofi_bin3"] else 0

            b_pow <- if ("abs_ofi_x_ofi" %in% names(theta_vec)) theta_vec["abs_ofi_x_ofi"] else 0

            # Calculate the imputed multiplier components. TODO: FIXED!
            m_bin <- b_bin1 * df$ofi_bin1 + b_bin2 * df$ofi_bin2 + b_bin3 * df$ofi_bin3
            m_pow <- b_pow * df$abs_ofi
            m_total <- m_bin + m_pow


            # Share = Cov(Bin Component, Total) / Var(Total)
            share <- cov(m_bin, m_total) / var(m_total)
            return(share)
        }

        # 5. Point Estimate & Delta Method Standard Errors
        share_point <- calc_share_bin(theta_hat, data)

        # Gradient of the share with respect to the estimated coefficients
        grad_s <- grad(func = calc_share_bin, x = theta_hat, df = data)

        # SE = sqrt( g' * Sigma * g )
        se_share <- sqrt(as.numeric(t(grad_s) %*% sigma_hat %*% grad_s))

        # 6. Construct Decomposition Table
        # Extract variables statically to calculate raw covariances
        b_bin1 <- if ("ofi_bin1" %in% names(theta_hat)) theta_hat["ofi_bin1"] else 0
        b_bin2 <- if ("ofi_bin2" %in% names(theta_hat)) theta_hat["ofi_bin2"] else 0
        b_bin3 <- if ("ofi_bin3" %in% names(theta_hat)) theta_hat["ofi_bin3"] else 0
        b_pow <- if ("abs_ofi_x_ofi" %in% names(theta_hat)) theta_hat["abs_ofi_x_ofi"] else 0

        m_bin_vec <- b_bin1 * data$ofi_bin1 + b_bin2 * data$ofi_bin2 + b_bin3 * data$ofi_bin3
        m_pow_vec <- b_pow * data$abs_ofi
        m_total_vec <- m_bin_vec + m_pow_vec


        decomp_list[[spec_idx]] <- data.table(
            type       = this_type,
            spec_idx   = spec_idx,
            # Raw variance stats
            cov_bin    = cov(m_bin_vec, m_total_vec),
            cov_pow    = cov(m_pow_vec, m_total_vec),
            var_total  = var(m_total_vec),
            # Shares and Inference
            share_bin  = share_point,
            share_pow  = 1 - share_point,
            share_se   = se_share,
            tstat_bin  = share_point / se_share,
            tstat_pow  = (1 - share_point) / se_share
        )

        # Append multipliers and tag with spec_idx so they don't blend together
        tmp_mult <- fm_results$table[var %in% target_vars]
        tmp_mult[, spec_idx := spec_idx]
        mult_list[[spec_idx]] <- tmp_mult
    }

    # Return both the detailed multiplier table and the summary decomposition
    return(list(
        multipliers = rbindlist(mult_list),
        M_decomp    = rbindlist(decomp_list)
    ))
}

tic("regressions")
data_list <- split(data_all, by = "type")
results_list <- mclapply(split(data_all, by = "type"), p.process_one_type, mc.cores = detectCores() - 2)
multipliers <- rbindlist(lapply(results_list, `[[`, "multipliers"), fill = T)
M_decomp <- rbindlist(lapply(results_list, `[[`, "M_decomp"), fill = T)
toc()

to_dir <- "tmp/continuous_spec/"
dir.create(to_dir, recursive = T, showWarnings = F)
saveRDS(multipliers, paste0(to_dir, "multipliers.RDS"))
saveRDS(M_decomp, paste0(to_dir, "M_decomp.RDS"))
