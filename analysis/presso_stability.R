# MR-PRESSO stability: rerun with different seeds / numbers of simulations and record the outlier set and corrected estimate.
suppressMessages(library(MRPRESSO))
a <- commandArgs(TRUE); code <- a[1]; nsim <- as.integer(a[2]); seed <- as.integer(a[3])
d <- read.csv(sprintf("mr_v4_inputs/forward_selected_r2_0.001_%s.csv", code))
set.seed(seed); t0 <- Sys.time()
pr <- mr_presso(BetaOutcome = "by", BetaExposure = "bx", SdOutcome = "sy", SdExposure = "sx", OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
                data = data.frame(bx = d$bx, sx = d$sx, by = d$by, sy = d$sy), NbDistribution = nsim, SignifThreshold = 0.05)
mr <- pr$`Main MR results`; out <- pr$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
row <- if (!is.na(mr[2, "Causal Estimate"])) 2 else 1
res <- data.frame(code = code, nsim = nsim, seed = seed, n_outliers = length(out), outliers = paste(d$rsid[out], collapse = ";"),
                  b = mr[row, "Causal Estimate"], se = mr[row, "Sd"], p = mr[row, "P-value"], global_p = as.character(pr$`MR-PRESSO results`$`Global Test`$Pvalue),
                  minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
write.csv(res, sprintf("presso_stab/%s_%d_%d.csv", code, nsim, seed), row.names = FALSE); print(res[, c(1:4, 6, 10)])
