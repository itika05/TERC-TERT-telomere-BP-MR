# v6: SuSiE/coloc.susie sensitivity to the LD definition: signed phased-haplotype LD (as in v4/v5) versus genotype-dosage LD in the same
# 503 individuals; L = 10; all credible-set pairs retained with identifiers; plus GLS IVW with the dosage LD matrix.
suppressMessages({library(coloc); library(susieR); library(MendelianRandomization)})
TR <- c("SBP", "DBP", "HYPERTENSION", "FG_HYPERTENSION"); diag <- list(); cs_out <- list(); sus <- list()
for (loc in c("TERC", "TERT")) for (w in c(100, 250)) for (t in TR) for (ld in c("haplotype", "dosage")) {
  x <- read.csv(sprintf("coloc_v4_inputs/%s_%dkb_%s.csv", loc, w, t)); x <- x[x$in_panel %in% c(TRUE, "True"), ]
  f <- if (ld == "haplotype") sprintf("coloc_v4_inputs/%s_%dkb_%s_LD.csv", loc, w, t) else sprintf("coloc_v6_inputs/%s_%dkb_%s_LDdos.csv", loc, w, t)
  LD <- as.matrix(read.csv(f, row.names = 1, check.names = FALSE)); stopifnot(all(rownames(LD) == x$varId))
  fits <- list()
  for (role in c("LTL", t)) {
    if (role == "LTL") { z <- x$b_ltl / x$se_ltl; n <- 464716 } else { z <- x$b_out / x$se_out; n <- median(x$n_out) }
    s_est <- estimate_s_rss(z, LD, n = n); kr <- kriging_rss(z, LD, n = n)
    nbad <- sum(kr$conditional_dist$logLR > 2 & abs(kr$conditional_dist$z) > 2)
    fit <- tryCatch(susie_rss(z = z, R = LD, n = n, L = 10, max_iter = 500, estimate_residual_variance = FALSE), error = function(e) NULL)
    ncs <- if (is.null(fit) || is.null(fit$sets$cs)) 0 else length(fit$sets$cs)
    lowz <- 0
    if (ncs > 0) for (i in seq_along(fit$sets$cs)) { v <- fit$sets$cs[[i]]; lead <- v[which.max(fit$pip[v])]; if (abs(z[lead]) < 2) lowz <- lowz + 1
      cs_out[[length(cs_out) + 1]] <- data.frame(ld = ld, locus = loc, window_kb = w, pair = t, trait = role, cs = names(fit$sets$cs)[i], size = length(v), lead_varId = x$varId[lead], lead_rsid = x$rsid[lead],
        lead_pip = fit$pip[lead], lead_abs_z = abs(z[lead]), purity_min_abs_corr = fit$sets$purity[i, "min.abs.corr"], members = paste(x$rsid[v], collapse = ";")) }
    diag[[length(diag) + 1]] <- data.frame(ld = ld, locus = loc, window_kb = w, pair = t, trait = role, n_variants = length(z), s_estimate = s_est, kriging_flagged = nbad,
      n_credible_sets = ncs, cs_led_by_abs_z_lt_2 = lowz, converged = if (is.null(fit)) NA else fit$converged, saturated = ncs >= 10)
    if (!is.null(fit)) { names(fit$pip) <- x$varId; colnames(fit$alpha) <- x$varId; colnames(fit$lbf_variable) <- x$varId }
    fits[[role]] <- fit }
  a <- fits[["LTL"]]; b <- fits[[t]]
  if (is.null(a) || is.null(b) || is.null(a$sets$cs) || is.null(b$sets$cs)) { sus[[length(sus) + 1]] <- data.frame(ld = ld, locus = loc, window_kb = w, trait = t, hit1 = NA, hit2 = NA, PP3 = NA, PP4 = NA, note = "no credible set in one or both traits"); next }
  r <- suppressMessages(coloc.susie(a, b, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)); s <- r$summary
  sus[[length(sus) + 1]] <- data.frame(ld = ld, locus = loc, window_kb = w, trait = t, hit1 = s$hit1, hit2 = s$hit2, PP3 = s$PP.H3.abf, PP4 = s$PP.H4.abf, note = "")
  cat(loc, w, t, ld, "\n")
}
write.csv(do.call(rbind, diag), "step8_susie_ld_sensitivity_diagnostics.csv", row.names = FALSE)
write.csv(do.call(rbind, cs_out), "step8_susie_ld_sensitivity_credible_sets.csv", row.names = FALSE)
write.csv(do.call(rbind, sus), "step8_susie_ld_sensitivity_coloc.csv", row.names = FALSE)
# GLS with dosage LD
out <- list()
for (code in c("S_mvpEU", "D_mvpEU", "H_mvpEU")) for (ld in c("haplotype", "dosage")) {
  d1 <- read.csv(sprintf("mr_v4_inputs/forward_selected_r2_0.1_%s.csv", code))
  Cm <- as.matrix(read.csv(if (ld == "haplotype") "mr_v4_inputs/forward_r2_0.1_correlation.csv" else "mr_v6_inputs/forward_r2_0.1_correlation_dosage.csv", row.names = 1, check.names = FALSE))
  iv <- mr_ivw(mr_input(bx = d1$bx, bxse = d1$sx, by = d1$by, byse = d1$sy, correlation = Cm[d1$varId, d1$varId]), model = "random", correl = TRUE, distribution = "t-dist")
  out[[length(out) + 1]] <- data.frame(code = code, ld = ld, b = iv@Estimate, se = iv@StdError, lo = iv@CILower, hi = iv@CIUpper, p = iv@Pvalue, nsnp = nrow(d1)) }
G <- do.call(rbind, out); write.csv(G, "step8_gls_ld_sensitivity.csv", row.names = FALSE); print(G, digits = 3)
