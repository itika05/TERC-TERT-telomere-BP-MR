# v9 (review M4 and M2): robust MR estimators that model a pleiotropic component, and the Nakao-weight sensitivity analysis.
# Inputs: the primary 116-instrument files (mr_v4_inputs/forward_selected_r2_0.001_*.csv) and the Nakao-weight files
# (mr_v9_inputs/, from step11_nakao_weights.py). Packages: MendelianRandomization 0.10.0 (mr_conmix, mr_cML, mr_ivw, mr_median),
# mr.raps from GitHub qingyuanzhao/mr.raps (commit dd79b5b), installed with its rsnps import removed (rsnps is used only by a
# plotting helper and is not installable here; the estimating functions are unchanged). RAPS: mr.raps.overdispersed.robust (Huber loss,
# k = 1.345, over-dispersion; sandwich SE), i.e. the estimator in Zhao et al. 2020 Ann Stat; the wrapper mr.raps() now also fits a
# shrinkage mixture model that failed on these data, so the documented robust over-dispersed function is called directly.
suppressMessages({library(MendelianRandomization); library(mr.raps)})
NOUT <- c(S_mvpEU = 425740, D_mvpEU = 425743, H_mvpEU = 318398, H_fg12 = 500264)   # outcome N used only by mr_cML (sample size)
NEXP <- 464716                                                                          # Codd et al. BOLT-LMM N
out <- list(); add <- function(...) out[[length(out) + 1]] <<- data.frame(...)
for (code in names(NOUT)) {
  d <- read.csv(sprintf("mr_v4_inputs/forward_selected_r2_0.001_%s.csv", code)); stopifnot(nrow(d) == 116)
  inp <- mr_input(bx = d$bx, bxse = d$sx, by = d$by, byse = d$sy, snps = d$varId)
  iv <- mr_ivw(inp, model = "random", distribution = "t-dist")
  add(code = code, weights = "Codd 2021", method = "IVW (MRE, t) [reference]", b = iv$Estimate, se = iv$StdError, lo = iv$CILower, hi = iv$CIUpper, p = iv$Pvalue, k = iv$SNPs, note = "")
  cm <- mr_conmix(inp, psi = 0, CIMin = -1, CIMax = 1, CIStep = 0.0005)
  ncm <- length(cm$CILower)
  add(code = code, weights = "Codd 2021", method = "Contamination mixture", b = cm$Estimate, se = NA, lo = min(cm$CILower), hi = max(cm$CIUpper), p = cm$Pvalue, k = cm$SNPs,
      note = sprintf("profile-likelihood CI (grid step 0.0005)%s; psi = %.4f (package default: 1.5 x SD of ratio estimates)", ifelse(ncm > 1, sprintf("; CI is a union of %d intervals", ncm), ""), cm$Psi))
  cml <- mr_cML(inp, MA = TRUE, DP = TRUE, n = min(NEXP, NOUT[[code]]), num_pert = 200, random_seed = 314)
  add(code = code, weights = "Codd 2021", method = "cML-MA-BIC-DP", b = cml$Estimate, se = cml$StdError, lo = cml$CILower, hi = cml$CIUpper, p = cml$Pvalue, k = cml$SNPs,
      note = sprintf("n = %d; 200 data perturbations, seed 314", min(NEXP, NOUT[[code]])))
  rp <- mr.raps.overdispersed.robust(d$bx, d$by, d$sx, d$sy, loss.function = "huber", suppress.warning = TRUE)
  add(code = code, weights = "Codd 2021", method = "MR-RAPS (Huber loss, over-dispersion)", b = rp$beta.hat, se = rp$beta.se, lo = rp$beta.hat - 1.959964 * rp$beta.se, hi = rp$beta.hat + 1.959964 * rp$beta.se,
      p = rp$beta.p.value, k = nrow(d), note = sprintf("tau2 = %.3g; normal CI", rp$tau2.hat))
  # Nakao-weight sensitivity (same instruments, same outcome rows)
  n <- read.csv(sprintf("mr_v9_inputs/forward_nakaoEU_weights_%s.csv", code)); stopifnot(identical(n$varId, d$varId), all.equal(n$by, d$by) == TRUE)
  ninp <- mr_input(bx = n$bx, bxse = n$sx, by = n$by, byse = n$sy, snps = n$varId)
  niv <- mr_ivw(ninp, model = "random", distribution = "t-dist")
  add(code = code, weights = "Nakao 2026 EU (portal)", method = "IVW (MRE, t)", b = niv$Estimate, se = niv$StdError, lo = niv$CILower, hi = niv$CIUpper, p = niv$Pvalue, k = niv$SNPs,
      note = sprintf("I2 %.0f%%", 100 * max(0, (niv$Heter.Stat[1] - (niv$SNPs - 1)) / niv$Heter.Stat[1])))
  nwm <- mr_median(ninp, weighting = "weighted", iterations = 1000, seed = 314159)
  add(code = code, weights = "Nakao 2026 EU (portal)", method = "Weighted median", b = nwm$Estimate, se = nwm$StdError, lo = nwm$CILower, hi = nwm$CIUpper, p = nwm$Pvalue, k = nwm$SNPs, note = "")
  neg <- mr_egger(ninp, distribution = "t-dist")
  add(code = code, weights = "Nakao 2026 EU (portal)", method = "MR-Egger (t)", b = neg$Estimate, se = neg$StdError.Est, lo = neg$CILower.Est, hi = neg$CIUpper.Est, p = neg$Pvalue.Est, k = neg$SNPs,
      note = sprintf("intercept %.4f, P %.3f", neg$Intercept, neg$Pvalue.Int))
}
O <- do.call(rbind, out); write.csv(O, "step11_robust_mr.csv", row.names = FALSE); print(O[, 1:8], digits = 4)
