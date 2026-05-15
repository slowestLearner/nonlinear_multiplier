# --- modified from 2a_regression_fm.R
# estimate dispersion by num of bins
library(this.path)
setwd(this.path::this.dir())
source("../../utilities/runmefirst.R")
source("../../utilities/regressions.R")
options(width = 200)

# load regression data
tic("preparing data")
data_all <- readRDS("../../tmp/raw_data/reg_inputs/reg_table_static.RDS")[type == "OFI"]

# get control variable names
cdata <- readRDS("../../tmp/raw_data/controls/controls_classification.RDS")
controls_char <- cdata[control_type == "return-predictor", var]
controls_liq <- cdata[control_type == "liquidity", var]
controls_list <- c(controls_char, controls_liq)
rm(cdata)

# merge with some dispersion measurse
tmp <- readRDS("tmp/raw_files/ofi_dispersion.RDS")[, .(yyyymm, permno, ofi_sd, hhi_pow2)]
tmp <- merge(tmp, unique(data_all[, .(yyyymm, permno)]), by = c("yyyymm", "permno"))
tmp <- melt(tmp, id.vars = c("yyyymm", "permno"), variable.name = "disp_type", value.name = "disp") %>%
  mutate(disp_type = as.character(disp_type)) %>%
  setDT()

data_all <- merge(data_all, tmp, by = c("yyyymm", "permno"), allow.cartesian = TRUE)
rm(tmp)
toc()


# ###########################################################
# Fama-MacBeth
# ###########################################################


# worker function
# data <- data_list[[5]]
p.process_one_type <- function(data) {
  this_type <- data[1, type] # parse
  this_disp_type <- data[1, disp_type]
  this_bin_disp <- data[1, bin_disp]

  # control variables for different regression specifications
  control_formulas <- c(
    "1",
    paste0(c(controls_char), collapse = "+"),
    paste0(c(controls_char, controls_liq), collapse = "+")
  )

  out_all <- data.table() # save results here

  max_spec_idx <- length(control_formulas)
  # cov_all <- vector("list", max_spec_idx)
  # names(cov_all) <- paste0("spec_", 1:max_spec_idx)

  # introduce direct controls
  for (spec_idx in 1:length(control_formulas)) {
    # spec_idx <- 1
    ff_no_ofi <- paste0("ret ~ ", control_formulas[spec_idx])
    ff <- paste0("ret ~ ", control_formulas[spec_idx], " + ofi_bin1 + ofi_bin2 + ofi_bin3")

    out <- p.fama_macbeth(data, ff, compare_coefs = TRUE)[grepl("ofi_", var)]
    out_no_ofi <- p.fama_macbeth(data, ff_no_ofi, compare_coefs = FALSE)

    out[, r2_no_ofi := out_no_ofi[1, r2]]
    out[, spec_idx := spec_idx]

    # keep track
    out_all <- rbind(out_all, out)
  }


  # name the control variables being added
  tmp <- data.table(
    spec_idx = 1:3,
    var_added = c("none_init", "controls_char", "controls_char+controls_liq")
  )
  out_all <- merge(out_all, tmp, by = "spec_idx", all.x = T)

  # mark and return
  out_all[, type := this_type]
  out_all[, disp_type := this_disp_type]
  out_all[, bin_disp := this_bin_disp]
  return(out_all)
}

# nbins_list <- c(3, 5, 10)
nbins_list <- c(3)

to_dir <- "tmp/reg_by_dispersion/"
dir.create(to_dir, recursive = T, showWarnings = F)

for (nbins in nbins_list) {
  # nbins <- nbins_list[1]
  tic(paste0("regressions, nbins = ", nbins))


  # sort by dispersion
  data <- copy(data_all)
  data[, bin_disp := ntile(disp, nbins), .(yyyymm, type, disp_type)]

  # # optional: resort bins
  # data[, sd_ofi := sd(ofi), .(yyyymm, type, disp_type, bin_disp)]
  # data[, bin := ifelse(abs(ofi) > 2 * sd_ofi, 3, ifelse(abs(ofi) > sd_ofi, 2, 1)), .(yyyymm, type, disp_type, bin_disp)]
  # data[, ofi_bin1 := ifelse(bin == 1, ofi, 0)]
  # data[, ofi_bin2 := ifelse(bin == 2, ofi, 0)]
  # data[, ofi_bin3 := ifelse(bin == 3, ofi, 0)]
  # data[, sd_ofi := NULL]

  data_list <- split(data, by = c("type", "disp_type", "bin_disp"))
  out <- rbindlist(mclapply(data_list, p.process_one_type, mc.cores = nc))

  # # -- take a look
  # out <- copy(out_all)[spec_idx == 3]
  # out[, coef_round := round(coef, 2)]
  # out[, tstat_round := round(coef / se, 2)]
  # tmp <- unique(out[, .(var)])[c(1:3, 4, 7, 5)][, var_lab := paste0(.I, "_", var)]
  # out <- merge(out, tmp, by = "var")
  # rm(tmp)

  # dcast(out, var_lab ~ disp_type + bin_disp, value.var = "coef_round")
  # dcast(out, var_lab ~ disp_type + bin_disp, value.var = "tstat_round")

  # # get dispersion numbers
  # data <- copy(data_all)
  # data[, bin_disp := ntile(disp, nbins), .(yyyymm, type, disp_type)][, dummy := ""]
  # dcast(data[, mean(disp), .(dummy, disp_type, bin_disp)], dummy ~ disp_type + bin_disp, value.var = "V1")

  saveRDS(out, paste0(to_dir, "fm_coefs_nbins_", nbins, ".RDS"))
  toc()
}


# # --- SANITY: take a look at plotting results

# data <- readRDS("tmp/reg_by_dispersion/fm_coefs_nbins_5.RDS")[var %in% paste0('ofi_bin', 1:3) & spec_idx == 3]

# this_disp_type <- 'ofi_sd'
# this_disp_type <- 'hhi_pow2'
# ggplot(data[disp_type == this_disp_type], aes(x = bin_disp, y = coef, color = var)) +
#   geom_line() +
#   geom_point() +
#   theme_classic() +
#   labs(x = "Dispersion bin", y = "Coefficient") +
#   geom_hline(yintercept = 0, lty = 3) +
#   theme(text = element_text(size = 35), legend.position = c(.8, .8), legend.title = element_blank()) +
#   ggtitle(this_disp_type)
