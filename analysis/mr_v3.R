# MR estimation with established packages on the v3 instrument sets.
# MendelianRandomization 0.10.0 (Burgess; github.com/cran/MendelianRandomization), MRPRESSO 1.0 (github.com/rondolab/MR-PRESSO).
# Settings: IVW multiplicative random effects (model="random": residual SE scaled up only when > 1) with t(k-1) inference;
# MR-Egger random effects with t(k-2); weighted median (1000 bootstrap, seed 314159); weighted MBE (phi = 1, delta SE, 1000 bootstrap, seed 314159);
# MR-PRESSO (NbDistribution = 2400 > k/alpha = 2340 for k = 117; forward primary outcomes only, run separately (mode "presso"); SignifThreshold = 0.05, seed 314159); covariance-aware IVW/Egger (correl) on the r2 < 0.1 set.
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
MODE <- if (length(args)) args[1] else "main"
res <- list()
bind <- function(x) res[[length(res) + 1]] <<- x
fwd_codes <- c("S_mvpEU", "D_mvpEU", "H_mvpEU", "H_fg", "S_gh", "D_gh", "S_bbj", "D_bbj", "S_mvpAA", "D_mvpAA", "H_mvpAA", "S_mvpHS", "D_mvpHS", "H_mvpHS")
Cm <- as.matrix(read.csv("mr_v3_inputs/forward_r2_0.1_correlation.csv", row.names = 1, check.names = FALSE))
for (code in fwd_codes) {
  d <- read.csv(sprintf("mr_v3_inputs/forward_selected_r2_0.001_%s.csv", code))
  meta <- list(direction = "forward", code = code, trait = d$trait[1], source = d$source[1], family = d$family[1])
  if (MODE == "presso") { if (code == args[2]) { x <- run_one(d, "r2<0.001 (primary)", extra = meta, presso = TRUE); write.csv(x[grepl("PRESSO", x$method), ], sprintf("step3v3_presso_%s.csv", code), row.names = FALSE) }; next }
  bind(run_one(d, "r2<0.001 (primary)", extra = meta, presso = FALSE)); cat(code, "\n")
  if (code %in% c("S_mvpEU", "D_mvpEU", "H_mvpEU")) {
    s <- d[d$steiger_z > 0, ]; bind(run_one(s, "r2<0.001, Steiger-filtered", extra = meta)[1, ])
    s <- d[!as.logical(d$terc_tert_region), ]; bind(run_one(s, "r2<0.001, excluding TERC/TERT regions", extra = meta)[1, ])
    s <- d[!as.logical(d$palindromic_ambiguous), ]; bind(run_one(s, "r2<0.001, excluding ambiguous palindromic", extra = meta)[1, ])
    for (src in c("sx", "sx_nAF")) { inp <- mr_input(bx = d$bx, bxse = d[[src]], by = d$by, byse = d$sy)
      iv <- mr_ivw(inp, model = "random", weights = "delta", distribution = "t-dist")
      bind(c(list(analysis = ifelse(src == "sx", "r2<0.001, second-order weights, exposure SE from P", "r2<0.001, second-order weights, exposure SE from N and allele frequency"),
                  method = "IVW (MRE, t, second-order weights)", b = iv@Estimate, se = iv@StdError, lo = iv@CILower, hi = iv@CIUpper, p = iv@Pvalue, nsnp = nrow(d)), meta) |> as.data.frame()) }
    d1 <- read.csv(sprintf("mr_v3_inputs/forward_selected_r2_0.1_%s.csv", code))
    cc <- Cm[d1$varId, d1$varId]
    bind(run_one(d1, "r2<0.1, covariance-aware", corr = cc, extra = meta))
    bind(run_one(d1, "r2<0.1, ignoring correlation (for contrast)", extra = meta)[1, ])
  }
}
if (MODE == "presso") quit(save = "no")
for (ex in c("SBP", "DBP", "Hypertension")) for (oc in c("L_codd", "L_nkSA")) {
  d <- read.csv(sprintf("mr_v3_inputs/reverse_%s_%s.csv", ex, oc))
  meta <- list(direction = "reverse", code = paste(ex, oc, sep = "->"), trait = ex, source = d$outcome[1], family = ifelse(oc == "L_codd", "primary", "secondary"))
  bind(run_one(d, "r2<0.001 (primary)", extra = meta, presso = FALSE)); cat(ex, oc, "\n")
  if (ex != "Hypertension") { s <- d[d$steiger_z > 0, ]; bind(run_one(s, "r2<0.001, Steiger-filtered", extra = meta)[1, ]) }
}
for (f in list.files(".", "^step3v3_presso_.*csv$")) res[[length(res) + 1]] <- read.csv(f)
R <- rbind_fill(res)
write.csv(R, "step3v3_mr_results_R.csv", row.names = FALSE)
cat("done", nrow(R), "\n")
