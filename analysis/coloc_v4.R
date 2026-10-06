# Colocalization v4 (coloc 6.0.3, susieR 0.12.35): LTL (Codd 2021 UK Biobank, source summary statistics) versus MVP European SBP, DBP and
# hypertension and FinnGen R12 hypertension. coloc.abf: quantitative traits sdY = 1 (standardized / inverse-normal units), hypertension type "cc";
# p1 = p2 = 1e-4, p12 in {1e-6, 5e-6, 1e-5, 5e-5}; windows +/-50, 100, 250 kb around the max-|z| LTL variant.
# SuSiE (+/-100 and 250 kb): susie_rss with 1000 Genomes EUR LD, L = 10; LD-mismatch diagnostics (estimate_s_rss, kriging_rss); coloc.susie.
suppressMessages({library(coloc); library(susieR)})
TR <- c("SBP", "DBP", "HYPERTENSION", "FG_HYPERTENSION")
rd <- function(loc, w, t) read.csv(sprintf("coloc_v4_inputs/%s_%dkb_%s.csv", loc, w, t))
mk <- function(b, se, snp, pos, n, cc) { d <- list(beta = setNames(b, snp), varbeta = setNames(se^2, snp), snp = snp, position = pos, N = n)
  if (cc) d$type <- "cc" else { d$type <- "quant"; d$sdY <- 1 }; d }
abf <- list()
for (loc in c("TERC", "TERT")) for (w in c(50, 100, 250)) for (t in TR) {
  x <- rd(loc, w, t)
  A <- mk(x$b_ltl, x$se_ltl, x$varId, x$pos, 464716, FALSE)
  B <- mk(x$b_out, x$se_out, x$varId, x$pos, median(x$n_out), grepl("HYPERTENSION", t))
  for (p12 in c(1e-6, 5e-6, 1e-5, 5e-5)) {
    r <- suppressMessages(coloc.abf(A, B, p1 = 1e-4, p2 = 1e-4, p12 = p12)); s <- r$summary
    top <- r$results[which.max(r$results$SNP.PP.H4), ]
    abf[[length(abf) + 1]] <- data.frame(locus = loc, window_kb = w, trait = t, p12 = p12, nsnps = s[["nsnps"]], PP0 = s[["PP.H0.abf"]], PP1 = s[["PP.H1.abf"]],
      PP2 = s[["PP.H2.abf"]], PP3 = s[["PP.H3.abf"]], PP4 = s[["PP.H4.abf"]], top_snp_H4 = top$snp, top_snp_H4_pp = top$SNP.PP.H4,
      top_rsid = x$rsid[match(top$snp, x$varId)])
  } }
A <- do.call(rbind, abf); write.csv(A, "step6_coloc_v4_abf.csv", row.names = FALSE); cat("abf done\n")

diag <- list(); cs_out <- list(); sus <- list()
fitone <- function(z, LD, n) tryCatch(susie_rss(z = z, R = LD, n = n, L = 10, max_iter = 500, estimate_residual_variance = FALSE), error = function(e) NULL)
for (loc in c("TERC", "TERT")) for (w in c(100, 250)) for (t in TR) {
  x <- rd(loc, w, t); x <- x[x$in_panel %in% c(TRUE, "True"), ]
  LD <- as.matrix(read.csv(sprintf("coloc_v4_inputs/%s_%dkb_%s_LD.csv", loc, w, t), row.names = 1, check.names = FALSE))
  fits <- list()
  for (role in c("LTL", t)) {
    if (role == "LTL") { z <- x$b_ltl / x$se_ltl; n <- 464716 } else { z <- x$b_out / x$se_out; n <- median(x$n_out) }
    s_est <- estimate_s_rss(z, LD, n = n)
    kr <- kriging_rss(z, LD, n = n); nbad <- sum(kr$conditional_dist$logLR > 2 & abs(kr$conditional_dist$z) > 2)
    f <- fitone(z, LD, n); ncs <- if (is.null(f) || is.null(f$sets$cs)) 0 else length(f$sets$cs)
    diag[[length(diag) + 1]] <- data.frame(locus = loc, window_kb = w, pair = t, trait = role, n_variants = length(z), N = n, s_estimate = s_est,
      kriging_flagged = nbad, n_credible_sets = ncs, converged = if (is.null(f)) NA else f$converged)
    if (ncs > 0) for (i in seq_along(f$sets$cs)) { v <- f$sets$cs[[i]]; lead <- v[which.max(f$pip[v])]
      cs_out[[length(cs_out) + 1]] <- data.frame(locus = loc, window_kb = w, pair = t, trait = role, cs = names(f$sets$cs)[i], size = length(v),
        lead = x$varId[lead], lead_rsid = x$rsid[lead], lead_pip = f$pip[lead], lead_abs_z = abs(z[lead]),
        purity_min_abs_corr = f$sets$purity[i, "min.abs.corr"], members = paste(x$rsid[v], collapse = ";")) }
    if (!is.null(f)) { names(f$pip) <- x$varId; colnames(f$alpha) <- x$varId; colnames(f$lbf_variable) <- x$varId }
    fits[[role]] <- f
  }
  a <- fits[["LTL"]]; b <- fits[[t]]
  if (is.null(a) || is.null(b) || is.null(a$sets$cs) || is.null(b$sets$cs)) {
    sus[[length(sus) + 1]] <- data.frame(locus = loc, window_kb = w, trait = t, p12 = NA, hit1 = NA, hit2 = NA, nsnps = NA, PP3 = NA, PP4 = NA,
      note = "no credible set in one or both traits"); next }
  for (p12 in c(1e-6, 1e-5, 5e-5)) {
    r <- suppressMessages(coloc.susie(a, b, p1 = 1e-4, p2 = 1e-4, p12 = p12)); s <- r$summary
    sus[[length(sus) + 1]] <- data.frame(locus = loc, window_kb = w, trait = t, p12 = p12, hit1 = s$hit1, hit2 = s$hit2, nsnps = s$nsnps, PP3 = s$PP.H3.abf, PP4 = s$PP.H4.abf, note = "")
  }
}
write.csv(do.call(rbind, diag), "step6_coloc_v4_susie_diagnostics.csv", row.names = FALSE)
write.csv(do.call(rbind, cs_out), "step6_coloc_v4_susie_credible_sets.csv", row.names = FALSE)
write.csv(do.call(rbind, sus), "step6_coloc_v4_susie.csv", row.names = FALSE)
cat("susie done\n")
