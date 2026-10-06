# v5: conditional F by multi-start BFGS (global minimum; v4 used a single start and reported a local minimum);
# sensitivity excluding variants whose SE for any exposure or outcome came from the N/allele-frequency fallback.
# Multivariable MR (MendelianRandomization 0.10.0): LTL adjusted for BMI and lymphocyte count; outcomes MVP European SBP, DBP, hypertension.
# Conditional F statistics follow Sanderson, Spiller and Bowden (Stat Med 2021) with zero covariance between exposure estimates
# (the exposures come from different GWAS); Q-statistic for heterogeneity reported by mr_mvivw.
suppressMessages(library(MendelianRandomization))
d0 <- read.csv("mr_v4_inputs/mvmr_instruments.csv")
fb <- Reduce(`|`, lapply(c("LTL", "BMI", "LYM", "S_mvpEU", "D_mvpEU", "H_mvpEU"), function(k) d0[[paste0(k, "_p")]] <= 1e-300 | d0[[paste0(k, "_p")]] >= 0.999 | d0[[paste0(k, "_b")]] == 0))
write.csv(data.frame(varId = d0$varId[fb]), "step7_mvmr_fallback_variants.csv", row.names = FALSE)
ALL <- list()
for (instset in c("all 330 variants", "excluding variants with any fallback SE")) {
d <- if (instset == "all 330 variants") d0 else d0[!fb, ]
condF <- function(bx, sx, j) {
  k <- ncol(bx); L <- nrow(bx)
  obj <- function(delta) sum((bx[, j] - bx[, -j, drop = FALSE] %*% delta)^2 / (sx[, j]^2 + (sx[, -j, drop = FALSE]^2) %*% (delta^2)))
  st <- list(rep(0, k - 1), unname(coef(lm(bx[, j] ~ 0 + bx[, -j]))), rep(1, k - 1), rep(-1, k - 1))
  if (k - 1 == 2) st <- c(st, list(c(1, -1), c(-1, 1)))
  min(sapply(st, function(s0) optim(s0, obj, method = "BFGS", control = list(maxit = 2000, reltol = 1e-12))$value)) / (L - (k - 1))
}
out <- list()
sets <- list(c("LTL", "BMI", "LYM"), c("LTL", "BMI"), c("LTL", "LYM"))
out <- list()
for (oc in c("S_mvpEU", "D_mvpEU", "H_mvpEU")) for (s in sets) {
  bx <- as.matrix(d[, paste0(s, "_b")]); sx <- as.matrix(d[, paste0(s, "_se")])
  inp <- mr_mvinput(bx = bx, bxse = sx, by = d[[paste0(oc, "_b")]], byse = d[[paste0(oc, "_se")]], exposure = s, outcome = oc)
  iv <- mr_mvivw(inp, model = "random", distribution = "t-dist")
  eg <- mr_mvegger(inp, distribution = "t-dist")
  for (j in seq_along(s)) {
    out[[length(out) + 1]] <- data.frame(outcome = oc, model = paste(s, collapse = " + "), exposure = s[j], method = "MV-IVW (random, t)",
      b = iv@Estimate[j], se = iv@StdError[j], lo = iv@CILower[j], hi = iv@CIUpper[j], p = iv@Pvalue[j], nsnp = nrow(d),
      Q = iv@Heter.Stat[1], Q_p = iv@Heter.Stat[2], condF = condF(bx, sx, j), egger_intercept = NA, egger_intercept_p = NA)
    out[[length(out) + 1]] <- data.frame(outcome = oc, model = paste(s, collapse = " + "), exposure = s[j], method = "MV-Egger (t)",
      b = eg@Estimate[j], se = eg@StdError.Est[j], lo = eg@CILower.Est[j], hi = eg@CIUpper.Est[j], p = eg@Pvalue.Est[j], nsnp = nrow(d),
      Q = NA, Q_p = NA, condF = NA, egger_intercept = eg@Intercept, egger_intercept_p = eg@Pvalue.Int)
  }
  # univariable IVW on the same instrument set, for contrast
  u <- mr_ivw(mr_input(bx = d$LTL_b, bxse = d$LTL_se, by = d[[paste0(oc, "_b")]], byse = d[[paste0(oc, "_se")]]), model = "random", distribution = "t-dist")
  if (identical(s, sets[[1]])) out[[length(out) + 1]] <- data.frame(outcome = oc, model = "LTL only (same variants)", exposure = "LTL", method = "IVW (random, t)",
      b = u@Estimate, se = u@StdError, lo = u@CILower, hi = u@CIUpper, p = u@Pvalue, nsnp = nrow(d), Q = NA, Q_p = NA, condF = NA, egger_intercept = NA, egger_intercept_p = NA)
}
R <- do.call(rbind, out); R$instrument_set <- instset; ALL[[instset]] <- R }
R <- do.call(rbind, ALL); write.csv(R, "step7_mvmr_results.csv", row.names = FALSE)
print(R[R$method != "MV-Egger (t)" & R$exposure == "LTL", c("instrument_set", "outcome", "model", "b", "lo", "hi", "p", "nsnp", "condF")], digits = 3)
