# --- Estimate static Fama-MacBeth with an increasing number of controls
library(this.path)
setwd(this.path::this.dir())
source("../../utilities/runmefirst.R")
source("../../utilities/regressions.R")
options(width = 120)

# regression data. Do not use BMI, not enough data
tic("preparing data")
data_all <- readRDS("../../tmp/raw_data/reg_inputs/reg_table_static.RDS")[type %in% c("FIT", "OFI")]

# get control variable names
cdata <- readRDS("../../tmp/raw_data/controls/controls_classification.RDS")
controls_char <- cdata[control_type == "return-predictor", var]
controls_liq <- cdata[control_type == "liquidity", var]
controls_list <- c(controls_char, controls_liq)
rm(cdata)

# # get names for bmi controls
# tmp <- readRDS("../tmp/raw_data/controls/controls_for_BMI.RDS")
# controls_bmi <- setdiff(names(tmp), c("yyyymm", "permno"))
# rm(tmp)
toc()

# load industries
tmp <- readRDS("../../../../../data/stocks/controls/ff_refined_industry_dummies.RDS")
data_all <- merge(data_all, tmp, by = c("yyyymm", "permno"))

# jkp chars
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


# create control formulas
ff_list <- c(
  "1",
  paste0(c(controls_char), collapse = "+"),
  paste0(c(controls_char, controls_liq), collapse = "+")
)

# add industries
ff_ind_nums <- c(5, 12, 17, 30, 49)
for (ind_num in ff_ind_nums) {
  ff_list <- c(ff_list, paste(c(controls_char, controls_liq, paste0("as.factor(ind", ind_num, ")")), collapse = " + "))
}

# add jkp chars
ff <- ff_list[length(ff_list)]
for (jkp_idx in sort(unique(jkp_vars[, idx]))) {
  ff <- paste0(ff, " + ", paste0(jkp_vars[idx == jkp_idx, var], collapse = " + "))
  ff_list <- c(ff_list, ff)
}

specs <- data.table(
  ff = ff_list, spec_idx = 1:length(ff_list),
  var_added = c("none", "chars", "liq", paste0("ff_ind", ff_ind_nums), paste0("jkp_", 1:10))
)


# worker function to estimate regression with one type of demand variable. reg_spec = "nonlinear" or "stdev"
# data <- data_all[type == this_type]
# reg_spec <- 'stdev'
p.process_one_type <- function(data) {
  this_type <- data[1, type] # parse

  out_all <- data.table() # save results here


  # introduce direct controls
  for (this_spec_idx in specs[, spec_idx]) {
    # this_spec_idx <- 1
    tic(paste0("type = ", this_type, " spec_idx = ", this_spec_idx))


    ff <- paste0("ret ~ ", specs[spec_idx == this_spec_idx, ff], " + ofi_bin1 + ofi_bin2 + ofi_bin3")


    out <- p.fama_macbeth(data, ff, compare_coefs = T)
    out <- out[grepl("ofi_bin", var)]
    out[, spec_idx := this_spec_idx]

    # keep track
    out_all <- rbind(out_all, out)
    gc()
    toc()
  }

  out_all <- merge(out_all, specs[, .(spec_idx, var_added)], by = "spec_idx")
  out_all[, type := this_type]

  # return
  return(out_all)
}

# stdev-based specification---
types <- unique(data_all[, type])
out_all <- data.table()
for (this_type in types) {
  # this_type <- types[1]
  out <- p.process_one_type(data_all[type == this_type])
  out_all <- rbind(out_all, out)
  gc()
}

to_file <- "tmp/reg_results/reg_results_fm.RDS"
dir.create(dirname(to_file), recursive = T, showWarnings = F)
saveRDS(out_all, to_file)



# # === SANITY check. plot

out_all <- readRDS("tmp/reg_results/reg_results_fm.RDS")
out_all[, var_type := ifelse(var %in% paste0("ofi_bin", 1:3), "coef", "diff")]

this_type <- "FIT"
this_type <- "OFI"

out <- out_all[var_type == "coef" & type == this_type]
ggplot(out, aes(x = spec_idx, y = coef, fill = var)) +
  geom_line(aes(color = var), lwd = 2) +
  geom_point(aes(color = var), cex = 5) +
  geom_ribbon(aes(ymin = coef - se, ymax = coef + se), alpha = 0.2) +
  theme_classic() +
  labs(x = "Specification", y = "Coefficient") +
  geom_hline(yintercept = 0, lty = 3) +
  theme(text = element_text(size = 35), legend.position = c(.8, .8), legend.title = element_blank()) +
  ggtitle(this_type) +
  geom_vline(xintercept = c(3.5, 8.5), lty = 3, lwd = 2)

out <- out_all[var_type == "diff" & type == this_type]
ggplot(out, aes(x = spec_idx, y = coef, fill = var)) +
  geom_line(aes(color = var)) +
  geom_point(aes(color = var)) +
  geom_ribbon(aes(ymin = coef - se, ymax = coef + se), alpha = 0.2) +
  theme_classic() +
  labs(x = "Specification", y = "Coefficient") +
  geom_hline(yintercept = 0, lty = 3) +
  theme(text = element_text(size = 35), legend.position = c(.8, .2), legend.title = element_blank()) +
  ggtitle(this_type) +
  geom_vline(xintercept = c(3.5, 8.5), lty = 3, lwd = 2)


# # --- SANITY: check with earlier reg results

# out_all <- readRDS("tmp/reg_results/reg_results_fm.RDS")
# out_all[, var_type := ifelse(var %in% paste0("ofi_bin", 1:3), "coef", "diff")]
# out_all <- out_all[spec_idx %in% 1:3]

# tmp <- readRDS('../../tmp/price_impact/regression_contemp/fm_stdev.RDS')
# tmp <- tmp[spec_idx %in% 1:3, .(spec_idx, var, coef, se, type)]
# out <- merge(out_all, tmp, by = c("spec_idx", "var", "type"))

# out[, cor(coef.x, coef.y, use = "complete.obs"), type]
# out[, cor(se.x, se.y, use = "complete.obs"), type]
