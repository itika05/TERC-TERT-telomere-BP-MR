# v7: MRPRESSO 1.0 run with its own code on the primary forward inputs, capturing the simulated data it generates so that the
# Python re-implementation can be tested on exactly the same random arrays (step9_presso_equivalence.py).
# The package function is copied unchanged except for two assign() calls that save randomData and RSSexp; results are those of the package.
# Usage: Rscript presso_equivalence.R <code> <NbDistribution> <seed>
suppressMessages(library(MRPRESSO))
a <- commandArgs(TRUE); code <- a[1]; B <- as.integer(a[2]); seed <- as.integer(a[3])
f <- MRPRESSO::mr_presso; bd <- as.list(body(f))
i_rand <- which(sapply(bd, function(e) any(grepl("^randomData <-", deparse(e)))))
i_rss <- which(sapply(bd, function(e) any(grepl("^RSSexp <-", deparse(e)))))
stopifnot(length(i_rand) == 1, length(i_rss) == 1)
bd <- append(bd, list(quote(assign("PRESSO_RSSexp", RSSexp, envir = .GlobalEnv))), after = i_rss)
bd <- append(bd, list(quote(assign("PRESSO_random", randomData, envir = .GlobalEnv))), after = i_rand)
body(f) <- as.call(bd); environment(f) <- asNamespace("MRPRESSO")
d <- read.csv(sprintf("mr_v4_inputs/forward_selected_r2_0.001_%s.csv", code))
dat <- data.frame(by = d$by, bx = d$bx, sy = d$sy, sx = d$sx)
t0 <- Sys.time()
r <- suppressWarnings(f(BetaOutcome = "by", BetaExposure = "bx", SdOutcome = "sy", SdExposure = "sx", data = dat, OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
                        SignifThreshold = 0.05, NbDistribution = B, seed = seed))
el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
stopifnot(length(PRESSO_random) == B)
arr <- do.call(rbind, lapply(PRESSO_random, function(m) c(m[, "bx"], m[, "by"])))   # B x (2k): exposure then outcome
dir.create("presso_equiv", showWarnings = FALSE)
con <- file(sprintf("presso_equiv/%s_B%d_s%d_random.bin", code, B, seed), "wb"); writeBin(as.vector(t(arr)), con, size = 8); close(con)
rss <- PRESSO_RSSexp; if (is.matrix(rss)) rss <- unlist(rss[1, ])
writeBin(as.numeric(rss), sprintf("presso_equiv/%s_B%d_s%d_RSSexp.bin", code, B, seed), size = 8)
ot <- r$`MR-PRESSO results`$`Outlier Test`; gt <- r$`MR-PRESSO results`$`Global Test`; mr <- r$`Main MR results`
dt <- r$`MR-PRESSO results`$`Distortion Test`
write.csv(data.frame(rsid = d$rsid, sign_bx = sign(d$bx), Dif2 = ot$RSSobs, Pvalue = as.character(ot$Pvalue)),
          sprintf("presso_equiv/%s_B%d_s%d_outlier_test.csv", code, B, seed), row.names = FALSE)
write.csv(data.frame(code = code, B = B, seed = seed, seconds = el, RSSobs = gt$RSSobs, global_P = as.character(gt$Pvalue),
                     n_outliers = if (is.character(dt$`Outliers Indices`)) 0 else length(dt$`Outliers Indices`),
                     outliers = if (is.character(dt$`Outliers Indices`)) "" else paste(d$rsid[dt$`Outliers Indices`], collapse = ";"),
                     b_corrected = mr[2, "Causal Estimate"], se_corrected = mr[2, "Sd"], p_corrected = mr[2, "P-value"],
                     distortion_coef = if (is.null(dt$`Distortion Coefficient`)) NA else dt$`Distortion Coefficient`, distortion_P = as.character(dt$Pvalue)),
          sprintf("presso_equiv/%s_B%d_s%d_package.csv", code, B, seed), row.names = FALSE)
cat(code, B, seed, "done in", round(el), "s\n")
