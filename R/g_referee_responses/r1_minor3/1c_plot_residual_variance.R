# --- Plot residual variance
library(this.path)
setwd(this.path::this.dir())
source("../../utilities/runmefirst.R")
options(width = 200)

# --- quarterly ones
data <- readRDS("tmp/pca_residuals/quarterly_oos.RDS")
setnames(data, "ret", "res_pc0")
data <- melt(data, id.vars = c("yyyymm", "permno"), variable.name = "num_pc", value.name = "res") %>%
    mutate(num_pc = as.integer(sub("res_pc", "", num_pc))) %>%
    setDT()
data <- data[, .(vv = var(res)), .(num_pc)][order(num_pc)]

pp <- ggplot(data, aes(x = num_pc, y = vv)) +
    geom_point() +
    geom_line() +
    labs(x = "Number of PCs", y = "Residual return variance") +
    theme_classic() +
    theme(text = element_text(size = 12)) +
    geom_hline(yintercept = 0, lty = 2)

to_file <- "../../output/figs/referee_responses/r1_minor3/pca/residual_variance_quarterly.png"
dir.create(dirname(to_file), showWarnings = FALSE, recursive = TRUE)
ggsave(to_file, pp, "png", w = 4, h = 3.5)
