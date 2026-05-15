# --- Use PCA residuals to do regressions. Also control for JKP characteristics
# this is a time-consuming script to run
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


# worker function to estimate regression with one type of demand variable. reg_spec = "nonlinear" or "stdev"
# data <- data_all[type == this_type]
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

# later can try to see if I can speed up the code
# types <- unique(data_all[, type])
# out_all <- data.table()

# for (this_type in types) {
#   tic(this_type)
#   # this_type <- types[1]
#   out <- p.process_one_type(data_all[type == this_type])
#   out_all <- rbind(out_all, out)
#   gc()
#   toc()
# }

to_file <- "tmp/reg_results/reg_results_fm_pca_oos_jkp.RDS"
dir.create(dirname(to_file), recursive = T, showWarnings = F)
saveRDS(out, to_file)

# # --- SANITY: check with earlier reg results
# new <- readRDS("tmp/reg_results/reg_results_fm_pca_oos_jkp.RDS")[y_var == "res_pc0" & jkp_vars == 0]
# old <- readRDS("../../tmp/price_impact/regression_contemp/fm_stdev.RDS")[spec_idx == 3 & grepl("ofi_bin", var) & type %in% c("BMI", "FIT", "OFI")]
# out <- merge(new, old, by = c("var", "type"))
# feols(coef.x ~ coef.y, data = out)
# feols(se.x ~ se.y, data = out)


# # === plot results

out_all <- readRDS("tmp/reg_results/reg_results_fm_pca_oos_jkp.RDS")
out_all[, num_pcs := as.integer(gsub("res_pc", "", y_var))]
out_all[, var_type := ifelse(var %in% paste0("ofi_bin", 1:3), "coef", "diff")]

# rank specifications
tmp <- out_all[, .(y_var, jkp_vars)] %>% unique()
tmp[, spec_idx := .I]
tmp[, spec_lab := paste0(y_var, ifelse(jkp_vars > 0, paste0(" + jkp_", jkp_vars), ""))]
out_all <- merge(out_all, tmp, by = c("y_var", "jkp_vars"))
rm(tmp)

this_type <- "FIT"
this_type <- "OFI"
this_type <- "BMI"
spec_cut_bmi <- 20

out <- copy(out_all[var_type == "coef" & type == this_type]) %>% setorder(var, spec_idx)
if (this_type == "BMI") {
  out <- out[spec_idx <= spec_cut_bmi]
}
p1 <- ggplot(
  out,
  aes(x = spec_idx, y = coef, color = var)
) +
  geom_line(lwd = 2) +
  geom_point(cex = 5) +
  geom_ribbon(aes(ymin = coef - se, ymax = coef + se, fill = var), alpha = 0.2, color = NA) +
  theme_classic() +
  scale_x_continuous(breaks = out[, spec_idx], labels = out[, spec_lab]) +
  labs(x = element_blank(), y = "Coefficient") +
  geom_hline(yintercept = 0, lty = 3) +
  theme(text = element_text(size = 25), legend.position = "bottom", legend.title = element_blank(), axis.text.x = element_text(angle = 45, vjust = .6)) +
  ggtitle(paste0(this_type, " coefs"))

# Plot 2: Differences
out <- copy(out_all[var_type == "diff" & type == this_type]) %>% setorder(var, spec_idx)
if (this_type == "BMI") {
  out <- out[spec_idx <= spec_cut_bmi]
}
p2 <- ggplot(
  out,
  aes(x = spec_idx, y = coef, color = var)
) +
  geom_line(lwd = 2) +
  geom_point(cex = 5) +
  geom_ribbon(aes(ymin = coef - se, ymax = coef + se, fill = var), alpha = 0.2, color = NA) +
  theme_classic() +
  scale_x_continuous(breaks = out[, spec_idx], labels = out[, spec_lab]) +
  labs(x = element_blank(), y = "Coefficient") +
  geom_hline(yintercept = 0, lty = 3) +
  theme(text = element_text(size = 25), legend.position = "bottom", legend.title = element_blank(), axis.text.x = element_text(angle = 45, vjust = .6)) +
  ggtitle(paste0(this_type, " diffs"))

# plot
p1 / p2
