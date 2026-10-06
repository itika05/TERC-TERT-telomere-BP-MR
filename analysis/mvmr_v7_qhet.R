# v7: MVMR::qhet_mvmr (reference code, ref_MVMR/, GitHub commit 8be0d94) on the 318-variant primary MVMR set.
# Fixes the v6 export, which took only the first row of the returned data frame (qh[1, ]) and recycled it over the three exposures.
# The return object has one row per exposure ("Exposure i" = column i of BXGs, here LTL, BMI, LYM); its structure is saved and the
# export is checked row by row. Uncertainty: nonparametric bootstrap over variants (percentile 95% interval) at rho = 0, using a
# vectorised copy of the reference Qtemp objective (identical formula; checked to reproduce the reference point estimates). The
# reference BCa interval needs the 'boot' package, which is not installed here.
suppressMessages(library(MendelianRandomization))
for (f in list.files("ref_MVMR", "\\.R$", full.names = TRUE)) source(f)
d0 <- read.csv("mr_v4_inputs/mvmr_instruments.csv")
fb <- Reduce(`|`, lapply(c("LTL", "BMI", "LYM", "S_mvpEU", "D_mvpEU", "H_mvpEU"), function(k) d0[[paste0(k, "_p")]] <= 1e-300 | d0[[paste0(k, "_p")]] >= 0.999 | d0[[paste0(k, "_b")]] == 0))
d <- d0[!fb, ]; stopifnot(nrow(d) == 318)
ex <- c("LTL", "BMI", "LYM"); bx <- as.matrix(d[, paste0(ex, "_b")]); sx <- as.matrix(d[, paste0(ex, "_se")])
stopifnot(identical(colnames(bx), paste0(ex, "_b")))
# vectorised Qtemp (same objective as ref_MVMR/qhet_mvmr.R: w = seY^2 + b' Sigma_l b + tau2; tau2 from optimize on (-10, 10); optim Nelder-Mead from 0)
qfast <- function(gy, sy, px, sxm, pcor) {
  L <- nrow(px); k <- ncol(px)
  wq <- function(b, tau2) { sb <- sweep(sxm, 2, b, `*`); sy^2 + rowSums((sb %*% pcor) * sb) + tau2 }
  PL <- function(tau2) { bc <- stats::optim(rep(0, k), function(b) sum((gy - px %*% b)^2 / wq(b, tau2)))$par
                         (sum((gy - px %*% bc)^2 / wq(bc, tau2)) - (L - 2))^2 }
  tau <- stats::optimize(PL, interval = c(-10, 10))$minimum
  stats::optim(rep(0, k), function(b) sum((gy - px %*% b)^2 / wq(b, tau)))$par
}
Q <- list(); struct <- c(); chk <- list(); RET <- list()
for (oc in c("S_mvpEU", "D_mvpEU", "H_mvpEU")) {
  fm <- format_mvmr(BXGs = bx, BYG = d[[paste0(oc, "_b")]], seBXGs = sx, seBYG = d[[paste0(oc, "_se")]], RSID = d$varId)
  stopifnot(identical(names(fm)[4:6], c("betaX1", "betaX2", "betaX3")), all.equal(unname(as.matrix(fm[, 4:6])), unname(bx)) == TRUE)
  for (rho in c(0, 0.1, 0.2, 0.3)) {
    R <- matrix(rho, 3, 3); diag(R) <- 1
    qh <- suppressWarnings(qhet_mvmr(fm, R, CI = FALSE, iterations = 10)); RET[[sprintf('%s_rho%.1f', oc, rho)]] <- list(object = qh, exposure_columns = colnames(bx), fm_names = names(fm))
    struct <- c(struct, sprintf("outcome %s rho %.1f: class %s; dim %s; rownames %s; colnames %s; values %s", oc, rho, paste(class(qh), collapse = "/"),
                               paste(dim(qh), collapse = "x"), paste(rownames(qh), collapse = "|"), paste(names(qh), collapse = "|"), paste(sprintf("%.15g", qh[, 1]), collapse = "|")))
    stopifnot(nrow(qh) == 3, identical(rownames(qh), paste("Exposure", 1:3)))
    out <- data.frame(outcome = oc, rho = rho, exposure = ex, exposure_row = rownames(qh), qhet_estimate = qh[["Effect Estimates"]])
    stopifnot(all(out$qhet_estimate == qh[1:3, 1]), length(unique(out$qhet_estimate)) == 3)
    fast <- qfast(fm$betaYG, fm$sebetaYG, bx, sx, R)
    chk[[length(chk) + 1]] <- data.frame(outcome = oc, rho = rho, exposure = ex, reference = qh[, 1], vectorised = fast, abs_diff = abs(qh[, 1] - fast))
    Q[[length(Q) + 1]] <- out
  }
}
QQ <- do.call(rbind, Q); CK <- do.call(rbind, chk)
writeLines(struct, "step9_qhet_return_objects.txt")
write.csv(CK, "step9_qhet_vectorised_check.csv", row.names = FALSE)
saveRDS(RET, "qhet_bootstrap/qhet_return_objects.rds")
print(CK, digits = 6); stopifnot(max(CK$abs_diff) < 1e-4)
# bootstrap at rho = 0
B <- as.integer(Sys.getenv("QHET_B", "1000")); BT <- list(); dir.create("qhet_bootstrap", showWarnings = FALSE)
# v8: resampling indices are drawn in the parent process from a fixed seed (one seed per outcome) before the parallel fits, so the
# replicates do not depend on the number of cores or on parallel RNG streams; every replicate (including failed fits) is saved.
for (oc in c("S_mvpEU", "D_mvpEU", "H_mvpEU")) {
  gy <- d[[paste0(oc, "_b")]]; sy <- d[[paste0(oc, "_se")]]
  seed <- 20260930 + match(oc, c("S_mvpEU", "D_mvpEU", "H_mvpEU")); set.seed(seed)
  IX <- lapply(1:B, function(i) sample.int(nrow(d), replace = TRUE))
  bs <- parallel::mclapply(1:B, function(i) { ix <- IX[[i]]
    tryCatch(c(qfast(gy[ix], sy[ix], bx[ix, , drop = FALSE], sx[ix, , drop = FALSE], diag(3)), 1), error = function(e) c(NA, NA, NA, 0)) }, mc.cores = 2)
  M <- do.call(rbind, bs); colnames(M) <- c(ex, "fit_ok")
  write.csv(data.frame(replicate = 1:B, seed = seed, M, first_index = sapply(IX, `[`, 1)), sprintf("qhet_bootstrap/%s_rho0_replicates.csv", oc), row.names = FALSE)
  ok <- M[, "fit_ok"] == 1
  for (j in 1:3) BT[[length(BT) + 1]] <- data.frame(outcome = oc, rho = 0, exposure = ex[j], bootstrap_B = B, n_failed = sum(!ok), seed = seed, boot_se = sd(M[ok, j]),
                                                   pct_lo = unname(quantile(M[ok, j], 0.025)), pct_hi = unname(quantile(M[ok, j], 0.975)))
}
BT <- do.call(rbind, BT)
QQ <- merge(QQ, BT, by = c("outcome", "rho", "exposure"), all.x = TRUE, sort = FALSE)
QQ <- QQ[order(QQ$outcome, QQ$rho, match(QQ$exposure, ex)), ]
write.csv(QQ, "step9_mvmr_qhet.csv", row.names = FALSE); print(QQ, digits = 4)
