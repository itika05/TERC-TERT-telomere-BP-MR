# MR estimation with established packages on the v3 instrument sets.
# MendelianRandomization 0.10.0 (Burgess; github.com/cran/MendelianRandomization), MRPRESSO 1.0 (github.com/rondolab/MR-PRESSO).
# Settings: IVW multiplicative random effects (model="random": residual SE scaled up only when > 1) with t(k-1) inference;
# MR-Egger random effects with t(k-2); weighted median (1000 bootstrap, seed 314159); weighted MBE (phi = 1, delta SE, 1000 bootstrap, seed 314159);
# MR-PRESSO (NbDistribution = 2400 > k/alpha = 2320 for k = 116; forward primary outcomes only, run separately (mode "presso"); SignifThreshold = 0.05, seed 314159); covariance-aware IVW/Egger (correl) on the r2 < 0.1 set.
suppressMessages({library(MendelianRandomization); library(MRPRESSO)})
args <- commandArgs(TRUE)
run_one <- function(d, label, corr = NULL, extra = list(), presso = FALSE) {
  inp <- if (is.null(corr)) mr_input(bx = d$bx, bxse = d$sx, by = d$by, byse = d$sy, snps = d$varId)
         else mr_input(bx = d$bx, bxse = d$sx, by = d$by, byse = d$sy, snps = d$varId, correlation = corr)
  out <- list()
  add <- function(method, est, se, lo, hi, p, k, extra2 = list()) out[[length(out) + 1]] <<- c(list(analysis = label, method = method, b = est, se = se, lo = lo, hi = hi, p = p, nsnp = k), extra, extra2)
  k <- nrow(d)
  if (is.null(corr)) {
    iv <- mr_ivw(inp, model = "random", distribution = "t-dist")
    q <- iv@Heter.Stat
    add("IVW (MRE, t)", iv@Estimate, iv@StdError, iv@CILower, iv@CIUpper, iv@Pvalue, k, list(Q = q[1], Q_p = q[2], I2 = max(0, (q[1] - (k - 1)) / q[1]) * 100))
    eg <- mr_egger(inp, robust = FALSE, penalized = FALSE, distribution = "t-dist")
    add("MR-Egger (t)", eg@Estimate, eg@StdError.Est, eg@CILower.Est, eg@CIUpper.Est, eg@Pvalue.Est, k,
        list(intercept = eg@Intercept, intercept_se = eg@StdError.Int, intercept_p = eg@Pvalue.Int, I2GX = eg@I.sq))
    set.seed(314159); wm <- mr_median(inp, weighting = "weighted", iterations = 1000, seed = 314159)
    add("Weighted median", wm@Estimate, wm@StdError, wm@CILower, wm@CIUpper, wm@Pvalue, k)
    mb <- mr_mbe(inp, weighting = "weighted", stderror = "delta", phi = 1, iterations = 1000, seed = 314159)
    add("Weighted mode (MBE)", mb@Estimate, mb@StdError, mb@CILower, mb@CIUpper, mb@Pvalue, k)
    if (presso) { set.seed(314159)
    df <- data.frame(bx = d$bx, sx = d$sx, by = d$by, sy = d$sy)
    pr <- tryCatch(mr_presso(BetaOutcome = "by", BetaExposure = "bx", SdOutcome = "sy", SdExposure = "sx", OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
                             data = df, NbDistribution = as.integer(Sys.getenv("PRESSO_N", "2400")), SignifThreshold = 0.05), error = function(e) NULL)
    if (!is.null(pr)) {
      mr <- pr$`Main MR results`; gp <- pr$`MR-PRESSO results`$`Global Test`$Pvalue
      nout <- if (is.null(pr$`MR-PRESSO results`$`Distortion Test`)) 0 else length(pr$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`)
      row <- if (!is.na(mr[2, "Causal Estimate"])) 2 else 1
      tq <- qt(0.975, k - nout - 1)
      add(paste0("MR-PRESSO (", ifelse(row == 2, "outlier-corrected", "raw"), ")"), mr[row, "Causal Estimate"], mr[row, "Sd"],
          mr[row, "Causal Estimate"] - tq * mr[row, "Sd"], mr[row, "Causal Estimate"] + tq * mr[row, "Sd"], mr[row, "P-value"], k - nout,
          list(global_p = as.numeric(sub("<", "", gp)), n_outliers = nout,
               distortion_p = if (is.null(pr$`MR-PRESSO results`$`Distortion Test`)) NA else as.numeric(sub("<", "", pr$`MR-PRESSO results`$`Distortion Test`$Pvalue))))
    } }
  } else {
    iv <- mr_ivw(inp, model = "random", correl = TRUE, distribution = "t-dist")
    add("IVW, correlated instruments (GLS)", iv@Estimate, iv@StdError, iv@CILower, iv@CIUpper, iv@Pvalue, k)
    eg <- mr_egger(inp, correl = TRUE, distribution = "t-dist")
    add("MR-Egger, correlated instruments (GLS)", eg@Estimate, eg@StdError.Est, eg@CILower.Est, eg@CIUpper.Est, eg@Pvalue.Est, k,
        list(intercept = eg@Intercept, intercept_p = eg@Pvalue.Int))
  }
  rbind_fill(lapply(out, function(x) as.data.frame(lapply(x, function(v) if (length(v) == 0 || is.null(v)) NA else v[1]), stringsAsFactors = FALSE)))
}
rbind_fill <- function(L) { L <- L[!sapply(L, is.null)]; cols <- unique(unlist(lapply(L, names)))
  do.call(rbind, lapply(L, function(d) { for (c in setdiff(cols, names(d))) d[[c]] <- NA; d[, cols] })) }
res <- list(); bind <- function(x) res[[length(res) + 1]] <<- x
# v6: reverse MR with source LTL statistics (custom and PLINK-clumped instruments); forward PLINK-clumped sets; FinnGen without ambiguous palindromes
for (ex in c("SBP", "DBP", "Hypertension")) for (sel in c("custom", "plink")) {
  d <- read.csv(sprintf("mr_v6_inputs/reverse_%s_%s_L_src.csv", sel, ex))
  meta <- list(direction = "reverse", code = paste0(ex, "->L_src"), trait = ex, source = "LTL (Codd 2021 source statistics)", family = ifelse(sel == "custom", "primary", "secondary"))
  lab <- ifelse(sel == "custom", "r2<0.001 (primary)", "PLINK 1.9 clumping")
  bind(run_one(d, lab, extra = meta)); cat(ex, sel, "\n")
  if (sel == "custom" && ex != "Hypertension") { s <- d[d$steiger_z > 0, ]; bind(run_one(s, "r2<0.001, Steiger-filtered", extra = meta)[1, ]) }
}
for (code in c("S_mvpEU", "D_mvpEU", "H_mvpEU", "H_fg12")) {
  d <- read.csv(sprintf("mr_v6_inputs/forward_plink_%s.csv", code))
  meta <- list(direction = "forward", code = code, trait = d$trait[1], source = d$source[1], family = "secondary")
  bind(run_one(d, "PLINK 1.9 clumping", extra = meta)); cat("plink", code, nrow(d), "\n")
}
d <- read.csv("mr_v4_inputs/forward_selected_r2_0.001_H_fg12.csv"); d <- d[!as.logical(d$palindromic_ambiguous), ]
bind(run_one(d, "r2<0.001, excluding ambiguous palindromic", extra = list(direction = "forward", code = "H_fg12", trait = "Hypertension", source = "FinnGen R12", family = "secondary"))[1, ])
R <- rbind_fill(res); write.csv(R, "step8_mr_v6_R.csv", row.names = FALSE); cat("done", nrow(R), "\n")
