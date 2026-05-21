# --- Estimate rolling multipliers over 3y windows
library(this.path)
setwd(this.path::this.dir())
source("../../utilities/runmefirst.R")
source("../../utilities/regressions.R")
options(width = 120)
nc <- detectCores() - 2

# load regression data
tic("preparing data")
data_all <- readRDS("../../tmp/raw_data/reg_inputs/reg_table_static.RDS")[type %in% c("BMI", "OFI", "FIT")]

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

# window length
data_all[, idx := frank(yyyymm, ties.method = "dense")]
window_length <- 12

# ###########################################################
# Fama-MacBeth
# ###########################################################

# worker function to estimate regression with one type of demand variable. reg_spec = "nonlinear" or "stdev"
# data <- copy(this_data)
p.process_one_type <- function(data) {
  this_type <- data[1, type] # parse
  last_yyyymm <- as.integer(data[.N, yyyymm])

  # control variables for different regression specifications
  control_formulas <- c(
    "1",
    paste0(c(controls_char), collapse = "+"),
    paste0(c(controls_char, controls_liq), collapse = "+")
  )

  out_all <- data.table() # save results here

  # loop over specifications
  for (spec_idx in 1:length(control_formulas)) {
    # spec_idx <- 1
    ff_no_ofi <- paste0("ret ~ ", control_formulas[spec_idx])
    ff <- paste0("ret ~ ", control_formulas[spec_idx], " + ofi_bin1 + ofi_bin2 + ofi_bin3")
    if (this_type == "BMI") {
      ff <- paste0(ff, "+", paste0(controls_bmi, collapse = "+"))
      ff_no_ofi <- paste0(ff_no_ofi, "+", paste0(controls_bmi, collapse = "+"))
    }

    out <- p.fama_macbeth(data, ff, compare_coefs = TRUE)
    out_no_ofi <- p.fama_macbeth(data, ff_no_ofi, compare_coefs = FALSE)

    out[, r2_no_ofi := out_no_ofi[1, r2]]
    out[, spec_idx := spec_idx]

    # keep track
    out_all <- rbind(out_all, out)
  }

  # return
  out_all[, type := this_type]
  out_all[, yyyymm := last_yyyymm]
  return(out_all)
}

# --- estimate
out_all <- data.table()
types <- unique(data_all[, type])
for (this_type in types) {
  # this_type <- types[1]
  data <- data_all[type == this_type]
  tmp <- unique(data[, .(idx)])
  first_idx <- tmp[1, idx]
  tmp <- tmp[idx >= (first_idx + window_length - 1)]

  # estimate regressios for ech window
  p.get_one_wrapper <- function(this_idx) {
    # this_idx <- tmp[1, idx]
    this_data <- data[idx %in% c(this_idx - 0:(window_length - 1))]
    out <- p.process_one_type(this_data)
    return(out)
  }

  tic()
  out <- rbindlist(mclapply(tmp[, idx], p.get_one_wrapper, mc.cores = nc))[grepl("ofi", var)][, date := as.Date(paste0(yyyymm, "01"), "%Y%m%d")]
  out_all <- rbind(out_all, out)
  toc()
}

# save
to_file <- "tmp/reg_results/rolling_multipliers.RDS"
dir.create(dirname(to_file), recursive = T, showWarnings = F)
saveRDS(out_all, to_file)


# # plot multiplier
# ggplot(out[var %in% paste0("ofi_bin", 1:3) & spec_idx == 3], aes(x = date, y = coef, fill = var)) +
#   geom_line(lwd = 2, aes(color = var)) +
#   geom_ribbon(aes(ymin = coef - 1.96 * se, ymax = coef + 1.96 * se), alpha = 0.2) +
#   theme_classic() +
#   theme(text = element_text(size = 30), legend.position = "bottom", legend.title = element_blank()) +
#   labs(x = element_blank(), y = element_blank()) +
#   ggtitle(paste0(this_type, " multipliers, ", window_length, "Q window")) +
#   geom_hline(yintercept = 0, lty = 2)

# # plot differences
# ggplot(out[!(var %in% paste0("ofi_bin", 1:3)) & spec_idx == 3], aes(x = date, y = coef, fill = var)) +
#   geom_line(lwd = 2, aes(color = var)) +
#   geom_ribbon(aes(ymin = coef - 1.96 * se, ymax = coef + 1.96 * se), alpha = 0.2) +
#   theme_classic() +
#   theme(text = element_text(size = 30), legend.position = "bottom", legend.title = element_blank()) +
#   labs(x = element_blank(), y = element_blank()) +
#   ggtitle(paste0(this_type, " multipliers differences, ", window_length, "Q window")) +
#   geom_hline(yintercept = 0, lty = 2)
