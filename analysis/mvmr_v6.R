# v6 MVMR. Primary instrument set: the 330 jointly clumped variants minus the 12 with any SE that could not be taken from P (318 variants;
# this removes the 5 hypertension outcome records that had a quantitative-trait fallback SE). Estimators: MendelianRandomization mr_mvivw
# (random effects, t) and mr_mvegger; univariable IVW on the same variants. mr_mvivw uses first-order weights (outcome SE only) and ignores
# exposure measurement error and covariance between exposure estimates. Conditional F: reference MVMR::strength_mvmr (GitHub commit 8be0d94)
# with gencov = 0 and with gencov_l = D_l R D_l for an assumed common correlation rho of exposure estimation errors (sample overlap);
# for comparison, the custom multi-start minimisation used in v5. Robust estimation: MVMR::qhet_mvmr (Q-statistic minimisation, which accounts
# for weak instruments given pcor) under the same rho values (point estimates; bootstrap CIs not computed).
suppressMessages(library(MendelianRandomization))
for (f in list.files("ref_MVMR", "\\.R$", full.names = TRUE)) source(f)
d0 <- read.csv("mr_v4_inputs/mvmr_instruments.csv")
fb <- Reduce(`|`, lapply(c("LTL", "BMI", "LYM", "S_mvpEU", "D_mvpEU", "H_mvpEU"), function(k) d0[[paste0(k, "_p")]] <= 1e-300 | d0[[paste0(k, "_p")]] >= 0.999 | d0[[paste0(k, "_b")]] == 0))
ex <- c("LTL", "BMI", "LYM")
condF_custom <- function(bx, sx, j) { k <- ncol(bx); L <- nrow(bx)
  obj <- function(delta) sum((bx[, j] - bx[, -j, drop = FALSE] %*% delta)^2 / (sx[, j]^2 + (sx[, -j, drop = FALSE]^2) %*% (delta^2)))
  st <- list(rep(0, k - 1), unname(coef(lm(bx[, j] ~ 0 + bx[, -j]))), rep(1, k - 1), rep(-1, k - 1), c(1, -1), c(-1, 1))
  min(sapply(st, function(s0) optim(s0, obj, method = "BFGS", control = list(maxit = 2000, reltol = 1e-12))$value)) / (L - (k - 1)) }
est <- list(); Fs <- list(); Q <- list()
for (set in c("primary (318; fallback SEs excluded)", "all 330 variants")) {
  d <- if (grepl("^primary", set)) d0[!fb, ] else d0
  bx <- as.matrix(d[, paste0(ex, "_b")]); sx <- as.matrix(d[, paste0(ex, "_se")])
  for (oc in c("S_mvpEU", "D_mvpEU", "H_mvpEU")) {
    inp <- mr_mvinput(bx = bx, bxse = sx, by = d[[paste0(oc, "_b")]], byse = d[[paste0(oc, "_se")]], exposure = ex, outcome = oc)
    iv <- mr_mvivw(inp, model = "random", distribution = "t-dist"); eg <- mr_mvegger(inp, distribution = "t-dist")
    u <- mr_ivw(mr_input(bx = d$LTL_b, bxse = d$LTL_se, by = d[[paste0(oc, "_b")]], byse = d[[paste0(oc, "_se")]]), model = "random", distribution = "t-dist")
    for (j in 1:3) est[[length(est) + 1]] <- data.frame(set = set, outcome = oc, exposure = ex[j], method = c("MV-IVW (random, t)", "MV-Egger (t)"),
      b = c(iv@Estimate[j], eg@Estimate[j]), se = c(iv@StdError[j], eg@StdError.Est[j]), lo = c(iv@CILower[j], eg@CILower.Est[j]), hi = c(iv@CIUpper[j], eg@CIUpper.Est[j]),
      p = c(iv@Pvalue[j], eg@Pvalue.Est[j]), nsnp = nrow(d), Q = c(iv@Heter.Stat[1], NA), Q_p = c(iv@Heter.Stat[2], NA), egger_intercept = c(NA, eg@Intercept), egger_intercept_p = c(NA, eg@Pvalue.Int))
    est[[length(est) + 1]] <- data.frame(set = set, outcome = oc, exposure = "LTL", method = "IVW univariable, same variants", b = u@Estimate, se = u@StdError, lo = u@CILower, hi = u@CIUpper, p = u@Pvalue,
      nsnp = nrow(d), Q = NA, Q_p = NA, egger_intercept = NA, egger_intercept_p = NA)
    if (grepl("^primary", set)) {
      fm <- format_mvmr(BXGs = bx, BYG = d[[paste0(oc, "_b")]], seBXGs = sx, seBYG = d[[paste0(oc, "_se")]], RSID = d$varId)
      for (rho in c(0, 0.1, 0.2, 0.3)) {
        R <- matrix(rho, 3, 3); diag(R) <- 1
        qh <- tryCatch(suppressWarnings(qhet_mvmr(fm, R, CI = FALSE, iterations = 10)), error = function(e) NULL)
        if (!is.null(qh)) Q[[length(Q) + 1]] <- data.frame(outcome = oc, rho = rho, exposure = ex, qhet_estimate = as.numeric(unlist(qh[1, ])))
      }
    }
  }
  if (grepl("^primary", set)) {
    fm <- format_mvmr(BXGs = bx, BYG = d$S_mvpEU_b, seBXGs = sx, seBYG = d$S_mvpEU_se, RSID = d$varId)
    for (rho in c(-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3)) {
      R <- matrix(rho, 3, 3); diag(R) <- 1
      gc <- if (rho == 0) 0 else lapply(seq_len(nrow(d)), function(l) outer(sx[l, ], sx[l, ]) * R)
      Fr <- suppressWarnings(capture.output(ff <- strength_mvmr(fm, gencov = gc)))
      Fs[[length(Fs) + 1]] <- data.frame(rho = rho, exposure = ex, F_reference_strength_mvmr = as.numeric(ff[1, ]),
        F_custom_multistart = if (rho == 0) sapply(1:3, function(j) condF_custom(bx, sx, j)) else NA, nsnp = nrow(d))
    }
  }
}
E <- do.call(rbind, est); write.csv(E, "step8_mvmr_v6_estimates.csv", row.names = FALSE)
F <- do.call(rbind, Fs); write.csv(F, "step8_mvmr_v6_conditionalF.csv", row.names = FALSE)
QQ <- do.call(rbind, Q); write.csv(QQ, "step8_mvmr_v6_qhet.csv", row.names = FALSE)
print(E[E$exposure == "LTL" & E$method != "MV-Egger (t)", c("set", "outcome", "method", "b", "lo", "hi", "p", "nsnp")], digits = 3); print(F, digits = 3); print(QQ[QQ$exposure == "LTL", ], digits = 3)
