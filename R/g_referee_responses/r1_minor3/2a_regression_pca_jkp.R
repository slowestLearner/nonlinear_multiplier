# --- Use PCA residuals to do regressions. Also control for JKP characteristics
# this is a time-consuming script to run
# NEW: made sure that when num of controls exceed the num of data points, we do not report results or obs
library(this.path)
setwd(this.path::this.dir())
source("../../utilities/runmefirst.R")
source("../../utilities/regressions.R")
options(width = 200)
library(patchwork)

# regression data. Can also use BMI
tic("preparing data")
data_all <- readRDS("../../tmp/raw_data/reg_inputs/reg_table_static.RDS")[type %in% c("FIT", "OFI", "BMI")]

# get alternative, PCA-residualized returns
data_all[, freq := ifelse(type == "BMI", "monthly", "quarterly")]
tmp <- readRDS("tmp/pca_residuals/quarterly_oos.RDS")[, freq := "quarterly"]
tmp <- rbind(tmp, readRDS("tmp/pca_residuals/monthly_oos.RDS")[, freq := "monthly"])
tmp[, ret := NULL]
data_all <- merge(data_all, tmp, by = c("yyyymm", "permno", "freq"), all.x = T)
setnames(data_all, "ret", "res_pc0")
rm(tmp)

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

# merge with jkp chars
from_dir <- "../../tmp/raw_data/jkp_chars_not_lagged/1_unif/"
files <- list.files(from_dir, pattern = "*.RDS")
jkp_vars <- data.table()
for (file in files) {
  # file <- files[1]
  idx <- as.integer(gsub("part_|.RDS", "", file))
  tic(file)
  tmp <- readRDS(paste0(from_dir, file))
  jkp_vars <- rbind(jkp_vars, unique(tmp[, .(var)])[, idx := idx])
  tmp <- dcast(tmp, yyyymm + permno ~ var, value.var = "char")
  if (file == files[1]) {
    jkp_data <- copy(tmp)
  } else {
    jkp_data <- merge(jkp_data, tmp, by = c("yyyymm", "permno"), all = T)
  }
  toc()
}
rm(from_dir, files, idx)

# need to lag these variables
tmp <- unique(jkp_data[, .(yyyymm)])
tmp[, mm := yyyymm %% 100]
tmp[, yyyymm_next := ifelse(mm == 12, yyyymm + 100 - 9, yyyymm + 3)][, mm := NULL]
jkp_data <- merge(jkp_data, tmp, by = "yyyymm")
jkp_data[, yyyymm := yyyymm_next][, yyyymm_next := NULL]
rm(tmp)

# merge together
data_all <- merge(data_all, jkp_data, by = c("yyyymm", "permno")) %>% na.omit()
rm(jkp_data)
gc()

# ###########################################################
# Fama-MacBeth
# ###########################################################

# NEW worker function that checks if the number of controls exceeds the number of data points, and if so, does not report results or obs
p.fama_macbeth <- function(data, ff, compare_coefs = FALSE) {
  # 1. Determine number of regressors (K) including intercept
  n_regressors <- length(attr(terms(as.formula(ff)), "term.labels")) + 1

  # regression for one period
  p.get_one_period <- function(this_ym) {
    sub_data <- data[yyyymm == this_ym]
    n_obs_this_month <- nrow(sub_data)

    # REQUIREMENT: Skip if observations (N) <= regressors (K)
    if (n_obs_this_month <= n_regressors) {
      return(NULL)
    }

    ols <- lm(ff, data = sub_data)
    dep_var_name <- all.vars(as.formula(ff))[1]

    return(data.table(
      yyyymm = this_ym,
      var = names(coef(ols)),
      coef = as.numeric(coef(ols)),
      r2 = var(ols$fitted.values) / var(sub_data[[dep_var_name]]),
      n_obs = n_obs_this_month # Store this to count later
    ))
  }

  # 2. Run regressions by period (rbindlist automatically handles NULLs by skipping them)
  out <- rbindlist(lapply(unique(data[, yyyymm]), p.get_one_period))

  if (nrow(out) == 0) {
    return(NULL)
  }

  # 3. Calculate Total Valid Observations
  # We look at unique (yyyymm, n_obs) to avoid double-counting vars within a month
  total_valid_obs <- sum(unique(out[, .(yyyymm, n_obs)])$n_obs)

  # Extract average R2 across valid periods
  r2_data <- unique(out[, .(yyyymm, r2)])[, .(r2 = mean(r2, na.rm = TRUE))]

  # Cleanup intermediate columns
  out[, `:=`(r2 = NULL, n_obs = NULL)]

  # --- Comparison Logic (Calculated based on valid out table) ---
  if (compare_coefs) {
    all_vars <- unique(out$var)
    target_vars <- sort(grep("ofi_", all_vars, value = TRUE))

    pairs_grid <- expand.grid(var_j = target_vars, var_i = target_vars, stringsAsFactors = FALSE)
    pairs_grid <- pairs_grid[pairs_grid$var_j != pairs_grid$var_i, ]

    diff_list <- lapply(1:nrow(pairs_grid), function(idx) {
      var_j <- pairs_grid$var_j[idx]
      var_i <- pairs_grid$var_i[idx]

      dt_j <- out[var == var_j, .(yyyymm, coef_j = coef)]
      dt_i <- out[var == var_i, .(yyyymm, coef_i = coef)]
      res <- merge(dt_j, dt_i, by = "yyyymm")

      if (nrow(res) > 0) {
        res[, coef := coef_j - coef_i]
        res[, var := paste0(var_j, " - ", var_i)]
        return(res[, .(yyyymm, var, coef)])
      } else {
        return(NULL)
      }
    })

    out <- rbind(out, rbindlist(diff_list))
  }

  # Newey-West lag
  this_hor <- if ("hor" %in% names(data)) data[1, hor] else 1

  # Summarize Results
  coef_data <- data.table()
  for (this_var in unique(out[, var])) {
    sub_out <- out[var == this_var & !is.na(coef)]

    if (nrow(sub_out) > 0) {
      mm <- lm(coef ~ 1, data = sub_out)
      se_val <- sqrt(vcov(mm)[1, 1])

      se_nw_val <- tryCatch(
        {
          sqrt(NeweyWest(mm, lag = this_hor, prewhite = FALSE)[1, 1])
        },
        error = function(e) {
          return(NA_real_)
        }
      )

      res_row <- data.table(
        var = this_var,
        coef = as.numeric(mm$coef[1]),
        se = se_val,
        se_nw = se_nw_val
      )
    } else {
      res_row <- data.table(var = this_var, coef = NA_real_, se = NA_real_, se_nw = NA_real_)
    }
    coef_data <- rbind(coef_data, res_row)
  }

  # 4. Final Meta-data
  coef_data[, `:=`(
    r2 = r2_data$r2,
    obs = total_valid_obs, # Now reflects ONLY used cross-sections
    type = data[1, type],
    nw_lag = this_hor
  )]

  return(coef_data)
}

# worker function to estimate regression with one type of demand variable. reg_spec = "nonlinear" or "stdev"
# data <- data_list[[3]]
p.process_one_type <- function(data) {
  this_type <- data[1, type] # parse

  # baseline controls
  controls <- paste0(c(controls_char, controls_liq), collapse = "+")

  out_all <- data.table() # save results here

  # change LHS variable
  vv <- names(data)[grepl("res_pc", names(data))]
  for (this_v in vv) {
    # this_v <- vv[1]

    ff <- paste0(this_v, " ~ ", controls, " + ofi_bin1 + ofi_bin2 + ofi_bin3")
    if (this_type == "BMI") {
      ff <- paste0(ff, "+", paste0(controls_bmi, collapse = "+"))
    }

    out <- p.fama_macbeth(data, ff, compare_coefs = T)[grepl("ofi_bin", var)]
    out[, y_var := this_v]
    out[, jkp_vars := 0]

    out_all <- rbind(out_all, out)
  }

  # add jkp chars progressively
  for (i in 1:10) {
    # i <- 1
    print(i)
    tic(i)
    this_jkp_vars <- jkp_vars[idx == i, var]
    ff <- paste0(ff, " + ", paste0(this_jkp_vars, collapse = " + "))

    out <- p.fama_macbeth(data, ff, compare_coefs = T)[grepl("ofi_bin", var)]
    out[, y_var := this_v]
    out[, jkp_vars := i]

    out_all <- rbind(out_all, out)
    gc()
    toc()
  }
  out_all[, type := this_type]

  return(out_all)
}

data_list <- split(data_all, by = "type")
rm(data_all)
out <- rbindlist(mclapply(data_list, p.process_one_type, mc.cores = 3))

to_file <- "tmp/reg_results/reg_results_fm_pca_oos_jkp_new.RDS"
dir.create(dirname(to_file), recursive = T, showWarnings = F)
saveRDS(out, to_file)

# # --- SANITY: check with earlier reg results
# new <- readRDS("tmp/reg_results/reg_results_fm_pca_oos_jkp_new.RDS")
# old <- readRDS("tmp/reg_results/reg_results_fm_pca_oos_jkp.RDS")

# data <- merge(new, old, by = c("var", "type", "y_var", "jkp_vars"))
# data[type == "BMI", cor(coef.x, coef.y, use = "complete.obs"), .(jkp_vars)]
# data[type == "BMI", .(obs.x = last(obs.x), obs.y = last(obs.y)), jkp_vars]

# old <- readRDS("../../tmp/price_impact/regression_contemp/fm_stdev.RDS")[spec_idx == 3 & grepl("ofi_bin", var) & type %in% c("BMI", "FIT", "OFI")]
# out <- merge(new, old, by = c("var", "type"))
# feols(coef.x ~ coef.y, data = out)
# feols(se.x ~ se.y, data = out)


# # === plot results

# out_all <- readRDS("tmp/reg_results/reg_results_fm_pca_oos_jkp.RDS")
# out_all[, num_pcs := as.integer(gsub("res_pc", "", y_var))]
# out_all[, var_type := ifelse(var %in% paste0("ofi_bin", 1:3), "coef", "diff")]

# # rank specifications
# tmp <- out_all[, .(y_var, jkp_vars)] %>% unique()
# tmp[, spec_idx := .I]
# tmp[, spec_lab := paste0(y_var, ifelse(jkp_vars > 0, paste0(" + jkp_", jkp_vars), ""))]
# out_all <- merge(out_all, tmp, by = c("y_var", "jkp_vars"))
# rm(tmp)

# this_type <- "FIT"
# this_type <- "OFI"
# this_type <- "BMI"
# spec_cut_bmi <- 20

# out <- copy(out_all[var_type == "coef" & type == this_type]) %>% setorder(var, spec_idx)
# if (this_type == "BMI") {
#   out <- out[spec_idx <= spec_cut_bmi]
# }
# p1 <- ggplot(
#   out,
#   aes(x = spec_idx, y = coef, color = var)
# ) +
#   geom_line(lwd = 2) +
#   geom_point(cex = 5) +
#   geom_ribbon(aes(ymin = coef - se, ymax = coef + se, fill = var), alpha = 0.2, color = NA) +
#   theme_classic() +
#   scale_x_continuous(breaks = out[, spec_idx], labels = out[, spec_lab]) +
#   labs(x = element_blank(), y = "Coefficient") +
#   geom_hline(yintercept = 0, lty = 3) +
#   theme(text = element_text(size = 25), legend.position = "bottom", legend.title = element_blank(), axis.text.x = element_text(angle = 45, vjust = .6)) +
#   ggtitle(paste0(this_type, " coefs"))

# # Plot 2: Differences
# out <- copy(out_all[var_type == "diff" & type == this_type]) %>% setorder(var, spec_idx)
# if (this_type == "BMI") {
#   out <- out[spec_idx <= spec_cut_bmi]
# }
# p2 <- ggplot(
#   out,
#   aes(x = spec_idx, y = coef, color = var)
# ) +
#   geom_line(lwd = 2) +
#   geom_point(cex = 5) +
#   geom_ribbon(aes(ymin = coef - se, ymax = coef + se, fill = var), alpha = 0.2, color = NA) +
#   theme_classic() +
#   scale_x_continuous(breaks = out[, spec_idx], labels = out[, spec_lab]) +
#   labs(x = element_blank(), y = "Coefficient") +
#   geom_hline(yintercept = 0, lty = 3) +
#   theme(text = element_text(size = 25), legend.position = "bottom", legend.title = element_blank(), axis.text.x = element_text(angle = 45, vjust = .6)) +
#   ggtitle(paste0(this_type, " diffs"))

# # plot
# p1 / p2
