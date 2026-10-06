# v10 (review of v9, major 2): single-signal colocalization of LTL with MVP European hypertension repeated with an allele-frequency
# model SE for the hypertension records, because their P-derived SEs fall with |z| (Supplementary Table S79). For each window,
# s = SE_P * sqrt(2 f (1 - f)) is computed over the variants with 1000 Genomes EUR MAF >= 5%, c = median(s), and the model SE is
# c / sqrt(2 f (1 - f)) for every variant (variants with MAF < 1% dropped). Everything else is as in the primary coloc.abf analysis
# (p1 = p2 = 1e-4; p12 = 1e-6, 5e-6, 1e-5, 5e-5; type "cc" with log-odds beta and its variance; LTL quantitative, sdY = 1).
suppressMessages(library(coloc))
out <- list()
for (loc in c("TERC", "TERT")) for (w in c(50, 100, 250)) {
  x <- read.csv(sprintf("coloc_v4_inputs/%s_%dkb_HYPERTENSION.csv", loc, w)); x <- x[!is.na(x$eur_af_alt) & x$se_out > 0, ]
  maf <- pmin(x$eur_af_alt, 1 - x$eur_af_alt); x <- x[maf >= 0.01, ]; maf <- pmin(x$eur_af_alt, 1 - x$eur_af_alt)
  pq <- 2 * x$eur_af_alt * (1 - x$eur_af_alt); s <- x$se_out * sqrt(pq); c0 <- median(s[maf >= 0.05]); se_m <- c0 / sqrt(pq)
  for (se_src in c("P-derived (primary)", "allele-frequency model")) for (p12 in c(1e-6, 5e-6, 1e-5, 5e-5)) {
    v <- if (se_src == "P-derived (primary)") x$se_out^2 else se_m^2
    d1 <- list(beta = x$b_ltl, varbeta = x$se_ltl^2, snp = x$varId, type = "quant", sdY = 1, N = 464716)
    d2 <- list(beta = x$b_out, varbeta = v, snp = x$varId, type = "cc", N = median(x$n_out))
    r <- suppressMessages(coloc.abf(d1, d2, p1 = 1e-4, p2 = 1e-4, p12 = p12))$summary
    out[[length(out) + 1]] <- data.frame(locus = loc, window_kb = w, se_source = se_src, p12 = p12, nsnps = r[["nsnps"]], c_median_s = c0,
                                         PP0 = r[["PP.H0.abf"]], PP1 = r[["PP.H1.abf"]], PP2 = r[["PP.H2.abf"]], PP3 = r[["PP.H3.abf"]], PP4 = r[["PP.H4.abf"]])
  }
}
O <- do.call(rbind, out); write.csv(O, "step12_htn_coloc_model_se.csv", row.names = FALSE)
print(aggregate(PP4 ~ locus + window_kb + se_source, O, function(z) sprintf("%.2f-%.2f", min(z), max(z))))
