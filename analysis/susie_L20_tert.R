# v5: SuSiE at TERT with L = 20 (L = 10 saturated for LTL in v4); otherwise identical to coloc_v4.R
suppressMessages({library(coloc); library(susieR)})
TR <- c("SBP", "DBP", "HYPERTENSION", "FG_HYPERTENSION")
rd <- function(loc, w, t) read.csv(sprintf("coloc_v4_inputs/%s_%dkb_%s.csv", loc, w, t))
fitone <- function(z, LD, n, L) tryCatch(susie_rss(z = z, R = LD, n = n, L = L, max_iter = 1000, estimate_residual_variance = FALSE), error = function(e) NULL)
diag <- list(); sus <- list(); cs_out <- list()
for (w in c(100, 250)) for (t in TR) {
  x <- rd("TERT", w, t); x <- x[x$in_panel %in% c(TRUE, "True"), ]
  LD <- as.matrix(read.csv(sprintf("coloc_v4_inputs/TERT_%dkb_%s_LD.csv", w, t), row.names = 1, check.names = FALSE))
  fits <- list()
  for (role in c("LTL", t)) {
    if (role == "LTL") { z <- x$b_ltl / x$se_ltl; n <- 464716 } else { z <- x$b_out / x$se_out; n <- median(x$n_out) }
    f <- fitone(z, LD, n, 20); ncs <- if (is.null(f) || is.null(f$sets$cs)) 0 else length(f$sets$cs)
    diag[[length(diag) + 1]] <- data.frame(locus = "TERT", window_kb = w, pair = t, trait = role, L = 20, n_variants = length(z), n_credible_sets = ncs,
      converged = if (is.null(f)) NA else f$converged, niter = if (is.null(f)) NA else f$niter)
    if (ncs > 0) for (i in seq_along(f$sets$cs)) { v <- f$sets$cs[[i]]; lead <- v[which.max(f$pip[v])]
      cs_out[[length(cs_out) + 1]] <- data.frame(locus = "TERT", window_kb = w, pair = t, trait = role, cs = names(f$sets$cs)[i], size = length(v),
        lead_rsid = x$rsid[lead], lead_pip = f$pip[lead], lead_abs_z = abs(z[lead]), purity_min_abs_corr = f$sets$purity[i, "min.abs.corr"]) }
    if (!is.null(f)) { names(f$pip) <- x$varId; colnames(f$alpha) <- x$varId; colnames(f$lbf_variable) <- x$varId }
    fits[[role]] <- f }
  a <- fits[["LTL"]]; b <- fits[[t]]
  if (is.null(a) || is.null(b) || is.null(a$sets$cs) || is.null(b$sets$cs)) {
    sus[[length(sus) + 1]] <- data.frame(locus = "TERT", window_kb = w, trait = t, p12 = 1e-5, hit1 = NA, hit2 = NA, PP3 = NA, PP4 = NA, maxPP4 = NA, note = "no credible set in one or both traits"); next }
  r <- suppressMessages(coloc.susie(a, b, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)); s <- r$summary
  k <- which.max(s$PP.H4.abf)
  sus[[length(sus) + 1]] <- data.frame(locus = "TERT", window_kb = w, trait = t, p12 = 1e-5, hit1 = s$hit1[k], hit2 = s$hit2[k], PP3 = s$PP.H3.abf[k], PP4 = s$PP.H4.abf[k],
     maxPP4 = max(s$PP.H4.abf), note = sprintf("%d pairs tested", nrow(s)))
  cat(w, t, "\n")
}
write.csv(do.call(rbind, diag), "step7_susie_L20_TERT_diagnostics.csv", row.names = FALSE)
write.csv(do.call(rbind, cs_out), "step7_susie_L20_TERT_credible_sets.csv", row.names = FALSE)
write.csv(do.call(rbind, sus), "step7_susie_L20_TERT_coloc.csv", row.names = FALSE)
