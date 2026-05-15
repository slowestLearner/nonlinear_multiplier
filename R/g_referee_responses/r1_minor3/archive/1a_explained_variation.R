# --- So it seems we still have to output the rethat. And later assess...
library(this.path)
setwd(this.path::this.dir())
source("../../utilities/runmefirst.R")

# Quarterly stock returns
data <- readRDS("../../tmp/raw_data/reg_inputs/reg_table_static.RDS")[, .(yyyymm, permno, ret)] %>% unique()

# get sizes
tmp <- readRDS("../../../../../data/stocks/prices/quarterly_return.RDS")[, .(yyyymm, permno, me_1)]
data <- merge(data, tmp, by = c("yyyymm", "permno"))
rm(tmp)

# make next period
data[, mm := yyyymm %% 100]
data <- data[, .(yyyymm = ifelse(mm == 3, yyyymm - 100 + 9, yyyymm - 3), permno, me_1, ret)]

# get industries
tmp <- readRDS("../../../../../data/stocks/controls/ff_refined_industry_dummies.RDS")
data <- merge(data, tmp, by = c("yyyymm", "permno"))

# existing chars
char_data <- readRDS("../../../../../data/stocks/controls/monthly_characteristics_not_lagged.RDS")
char_data[is.na(char_data)] <- 0
data <- merge(data, char_data, by = c("yyyymm", "permno"))
vv_chars <- setdiff(names(char_data), c("yyyymm", "permno"))

liq_data <- readRDS("../../../../../data/stocks/controls/monthly_liquidity_measures_not_lagged.RDS")[, c("realized_vol", "me") := NULL]
liq_data[is.na(liq_data)] <- 0
data <- merge(data, liq_data, by = c("yyyymm", "permno"))
vv_liq <- setdiff(names(liq_data), c("yyyymm", "permno"))
rm(char_data, liq_data)

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

# merge together
data <- merge(data, jkp_data, by = c("yyyymm", "permno")) %>% na.omit()
rm(jkp_data)
gc()

# create 3 types of weights
data[, w_ew := 1 / .N, yyyymm]
data[, w_sqrt_me := sqrt(me_1) / sum(sqrt(me_1)), yyyymm]
data[, w_me := me_1 / sum(me_1), yyyymm]

# process one month
p.get_one <- function(this_data) {
    # this_data <- data_list[[1]]

    tmp <- copy(this_data[, .(yyyymm, permno, ret, me_1)])

    num_vars <- data.table()

    weighting <- c("ew", "sqrt_me", "me")
    ff_ind_nums <- c(5, 12, 17, 30, 49)

    # create list of formulas
    ff_list <- c()
    ff_list <- c(ff_list, paste0("ret ~ ", paste(vv_chars, collapse = " + ")))
    ff_list <- c(ff_list, paste0("ret ~ ", paste(c(vv_chars, vv_liq), collapse = " + ")))
    for (ind_num in ff_ind_nums) {
        ff_list <- c(ff_list, paste0("ret ~ ", paste(c(vv_chars, vv_liq, paste0("as.factor(ind", ind_num, ")")), collapse = " + ")))
    }

    ff <- ff_list[length(ff_list)]
    for (jkp_idx in sort(unique(jkp_vars[, idx]))) {
        ff <- paste0(ff, " + ", paste0(jkp_vars[idx == jkp_idx, var], collapse = " + "))
        ff_list <- c(ff_list, ff)
    }

    specs <- data.table(
        ff = ff_list, spec_idx = 1:length(ff_list),
        var_added = c("chars", "liq", paste0("ff_ind", ff_ind_nums), paste0("jkp_", 1:10))
    )

    # keep track of predictions
    out <- data.table()
    for (this_spec_idx in specs[, spec_idx]) {
        # this_spec_idx <- 1
        # tic(this_spec_idx)
        formula <- specs[spec_idx == this_spec_idx, ff] %>% as.formula()
        mm_ew <- lm(formula, data = this_data)
        mm_sqrt_me <- lm(formula, data = this_data, weights = this_data[, sqrt(me_1)])
        mm_me <- lm(formula, data = this_data, weights = this_data[, me_1])

        tt <- copy(tmp)
        tt[, rethat_ew := mm_ew$fitted.values]
        tt[, rethat_sqrt_me := mm_sqrt_me$fitted.values]
        tt[, rethat_me := mm_me$fitted.values]
        tt[, spec_idx := this_spec_idx]
        tt[, num_vars := summary(mm_ew)$df[1]]
        out <- rbind(out, tt)
        # toc()
    }

    return(out)
}

data_list <- split(data, by = "yyyymm")
rm(data)
gc()

nc <- detectCores() - 2

tmp <- data.table(idx = 1:length(data_list))
block_size <- 2 * nc
tmp[, block_idx := ceiling(idx / block_size)]

# takes around 12 mins on a 8 core machine
out <- data.table()
for (this_block_idx in sort(unique(tmp[, block_idx]))) {
    # this_block_idx <- 1
    tic(this_block_idx)
    out <- rbind(out, rbindlist(mclapply(data_list[tmp[block_idx == this_block_idx, idx]], p.get_one, mc.cores = nc)))
    gc()
    toc()
}

# save locally for now
to_file <- "tmp/variance_explained/rethat.RDS"
dir.create(dirname(to_file), recursive = T, showWarnings = F)
saveRDS(out, to_file)

# --- take a look at variance explained, total

# this is all in-sample
out <- readRDS("tmp/variance_explained/rethat.RDS")
out[, w_ew := 1 / .N, yyyymm]
out[, w_sqrt_me := sqrt(me_1) / sum(sqrt(me_1)), yyyymm]
out[, w_me := me_1 / sum(me_1), yyyymm]

tt <- out[, .(type = "1_ew", r2 = 1 - sum((ret - rethat_ew)^2 * w_ew) / sum(ret^2 * w_ew)), spec_idx]
tt <- rbind(tt, out[, .(type = "2_sqrt_me", r2 = 1 - sum((ret - rethat_sqrt_me)^2 * w_sqrt_me) / sum(ret^2 * w_sqrt_me)), spec_idx])
tt <- rbind(tt, out[, .(type = "3_me", r2 = 1 - sum((ret - rethat_me)^2 * w_me) / sum(ret^2 * w_me)), spec_idx])

# get num of obs, etc
tt2 <- out[, .(obs_permno = length(unique(permno)), num_vars = last(num_vars)), .(yyyymm, spec_idx)]
tt2 <- tt2[, .(obs_permno = mean(obs_permno), num_vars = last(num_vars)), spec_idx]
tt <- merge(tt, tt2, by = "spec_idx")
rm(tt2)
tt[, adj_r2 := 1 - (1 - r2) * (obs_permno - 1) / (obs_permno - num_vars - 1)] # This is right?

ggplot(tt, aes(x = num_vars, y = adj_r2, color = type)) +
    geom_line(lwd = 2) +
    geom_hline(yintercept = 0, lty = 3) +
    theme(text = element_text(size = 35))


# --- can also look at average R2 by period

out <- readRDS("tmp/variance_explained/rethat.RDS")
out[, w_ew := 1 / .N, yyyymm]
out[, w_sqrt_me := sqrt(me_1) / sum(sqrt(me_1)), yyyymm]
out[, w_me := me_1 / sum(me_1), yyyymm]

tt <- out[, .(type = "1_ew", r2 = 1 - sum((ret - rethat_ew)^2 * w_ew) / sum(ret^2 * w_ew)), spec_idx]
tt <- rbind(tt, out[, .(type = "2_sqrt_me", r2 = 1 - sum((ret - rethat_sqrt_me)^2 * w_sqrt_me) / sum(ret^2 * w_sqrt_me)), spec_idx])
tt <- rbind(tt, out[, .(type = "3_me", r2 = 1 - sum((ret - rethat_me)^2 * w_me) / sum(ret^2 * w_me)), spec_idx])

# get num of obs, etc
tt2 <- out[, .(obs_permno = length(unique(permno)), num_vars = last(num_vars)), .(yyyymm, spec_idx)]
tt2 <- tt2[, .(obs_permno = mean(obs_permno), num_vars = last(num_vars)), spec_idx]
tt <- merge(tt, tt2, by = "spec_idx")
rm(tt2)
tt[, adj_r2 := 1 - (1 - r2) * (obs_permno - 1) / (obs_permno - num_vars - 1)] # This is right?

ggplot(tt, aes(x = num_vars, y = adj_r2, color = type)) +
    geom_line(lwd = 2) +
    geom_hline(yintercept = 0, lty = 3) +
    theme_classic() +
    theme(text = element_text(size = 35), legend.position = c(.8, .2), legend.title = element_blank()) +
    geom_vline(xintercept = c(22, 80), lty = 3, lwd = 2)
