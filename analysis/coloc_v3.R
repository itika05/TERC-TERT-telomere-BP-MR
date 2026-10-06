# coloc (Wallace; v6.0.3) and susieR (0.12.35) on v3 regional inputs.
# coloc.abf: quantitative traits sdY = 1 (standardised units; prior SD 0.15), hypertension type "cc" (prior SD 0.2, log OR);
# priors p1 = p2 = 1e-4, p12 in {1e-6, 5e-6, 1e-5, 5e-5}; windows +/-50, 100, 250, 500 kb around the max-|z| LTL variant.
# SuSiE: susie_rss on z-scores with the 1000G EUR signed LD matrix (503 individuals), L = 10, n = median N; coloc.susie on
# every pair of credible sets; LD-mismatch diagnostic: kriging_rss (susieR) and estimate_s_rss.
suppressMessages({library(coloc); library(susieR)})
rd <- function(loc, w, t, secol = "se") { x <- read.csv(sprintf("coloc_v3_inputs/%s_%dkb_%s.csv", loc, w, t)); x$use_se <- x[[secol]]; x }
ds <- function(x, t) {
  d <- list(beta = setNames(x$beta, x$varId), varbeta = setNames(x$use_se^2, x$varId), snp = x$varId, position = x$pos, N = median(x$n, na.rm = TRUE))
  if (t == "HYPERTENSION") { d$type <- "cc" } else { d$type <- "quant"; d$sdY <- 1 }
  d }
abf <- list()
for (loc in c("TERC", "TERT")) for (w in c(50, 100, 250, 500)) for (secol in c("se", "se_reported")) {
  L <- rd(loc, w, "LTL", secol)
  for (t in c("SBP", "DBP", "HYPERTENSION")) { B <- rd(loc, w, t, secol)
    for (p12 in c(1e-6, 5e-6, 1e-5, 5e-5)) {
      r <- coloc.abf(ds(L, "LTL"), ds(B, t), p1 = 1e-4, p2 = 1e-4, p12 = p12)
      s <- r$summary
      abf[[length(abf) + 1]] <- data.frame(locus = loc, window_kb = w, se_source = ifelse(secol == "se", "P-derived unless floored (primary)", "portal-reported"),
        trait = t, p12 = p12, nsnps = s[["nsnps"]], PP0 = s[["PP.H0.abf"]], PP1 = s[["PP.H1.abf"]], PP2 = s[["PP.H2.abf"]], PP3 = s[["PP.H3.abf"]], PP4 = s[["PP.H4.abf"]],
        top_snp_H4 = r$results$snp[which.max(r$results$SNP.PP.H4)])
    } } }
A <- do.call(rbind, abf); write.csv(A, "step4v3_coloc_abf_R.csv", row.names = FALSE); cat("abf done\n")

sus <- list(); cs_out <- list(); diag <- list()
for (loc in c("TERC", "TERT")) for (w in c(100, 250)) {
  LD <- as.matrix(read.csv(sprintf("coloc_v3_inputs/%s_%dkb_LD.csv", loc, w), row.names = 1, check.names = FALSE))
  fits <- list()
  for (t in c("LTL", "SBP", "DBP", "HYPERTENSION")) {
    x <- rd(loc, w, t); z <- x$beta / x$use_se; n <- median(x$n, na.rm = TRUE)
    s_est <- estimate_s_rss(z, LD, n = n)
    kr <- kriging_rss(z, LD, n = n); nbad <- sum(kr$conditional_dist$logLR > 2 & abs(kr$conditional_dist$z) > 2)
    f <- tryCatch(susie_rss(z = z, R = LD, n = n, L = 10, max_iter = 500, estimate_residual_variance = FALSE), error = function(e) NULL)
    ncs <- if (is.null(f) || is.null(f$sets$cs)) 0 else length(f$sets$cs)
    diag[[length(diag) + 1]] <- data.frame(locus = loc, window_kb = w, trait = t, n_variants = length(z), N = n, s_estimate = s_est,
                                           kriging_flagged = nbad, n_credible_sets = ncs, converged = if (is.null(f)) NA else f$converged)
    if (ncs > 0) for (i in seq_along(f$sets$cs)) { v <- f$sets$cs[[i]]; lead <- v[which.max(f$pip[v])]
      cs_out[[length(cs_out) + 1]] <- data.frame(locus = loc, window_kb = w, trait = t, cs = names(f$sets$cs)[i], size = length(v),
        lead = x$varId[lead], lead_rsid = x$rsid[lead], lead_pip = f$pip[lead], purity_min_abs_corr = f$sets$purity[i, "min.abs.corr"],
        members = paste(x$rsid[v], collapse = ";")) }
    f$snp_names <- x$varId; if (!is.null(f)) { names(f$pip) <- x$varId; colnames(f$alpha) <- x$varId; colnames(f$lbf_variable) <- x$varId }
    fits[[t]] <- f
  }
  for (t in c("SBP", "DBP", "HYPERTENSION")) {
    a <- fits[["LTL"]]; b <- fits[[t]]
    if (is.null(a) || is.null(b) || is.null(a$sets$cs) || is.null(b$sets$cs)) {
      sus[[length(sus) + 1]] <- data.frame(locus = loc, window_kb = w, trait = t, hit1 = NA, hit2 = NA, nsnps = NA, PP3 = NA, PP4 = NA,
        note = "no credible set in one or both traits"); next }
    for (p12 in c(1e-6, 1e-5, 5e-5)) {
      r <- coloc.susie(a, b, p1 = 1e-4, p2 = 1e-4, p12 = p12)
      s <- r$summary
      sus[[length(sus) + 1]] <- data.frame(locus = loc, window_kb = w, trait = t, p12 = p12, hit1 = s$hit1, hit2 = s$hit2, nsnps = s$nsnps,
                                           PP3 = s$PP.H3.abf, PP4 = s$PP.H4.abf, note = "")
    } }
}
write.csv(do.call(rbind, diag), "step4v3_susie_diagnostics.csv", row.names = FALSE)
write.csv(do.call(rbind, cs_out), "step4v3_susie_credible_sets.csv", row.names = FALSE)
write.csv(do.call(rbind, lapply(sus, function(x) { if (!"p12" %in% names(x)) x$p12 <- NA; x[, c("locus","window_kb","trait","p12","hit1","hit2","nsnps","PP3","PP4","note")] })), "step4v3_coloc_susie.csv", row.names = FALSE)
cat("susie done\n")
