# v12 (review of v11): conditional colocalization of LTL with MVP European hypertension, previously not run because the case
# fraction was not known. The portal effective N (318,397.62) equals 4/(1/cases + 1/controls) for dbGaP phs002453 analysis
# pha005545.1 (Phe_401 Hypertension, European; 322,448 cases, 105,690 controls; Supplementary Table S86), so s = 322448 / 428138.
# Settings are those of coloc_v9_conditional.R (coloc.signals, method "cond", iterative and allbutone modes, p1 = p2 = 1e-4,
# p12 = 1e-5, pthr = 1e-6, maxhits = 5; 1000 Genomes EUR haplotype and dosage LD; reference-panel MAF). Because the MVP
# hypertension P-derived SEs fall with |z| (Supplementary Table S79), each analysis is repeated with the allele-frequency model SE
# of coloc_v10_htn_modelSE.R (c / sqrt(2f(1-f)), c = median of SE x sqrt(2f(1-f)) over panel variants with MAF >= 5%).
# Note: for case-control data coloc:::est_cond recomputes the conditional variances from MAF, N and s (VMAF.cc), so the supplied
# varbeta enters only for variants excluded from conditioning (r^2 > 0.8 with a signal, beta set to 0); the two SE versions therefore
# give the same posteriors, and these analyses do not depend on the MVP P-derived SEs.
suppressMessages(library(coloc))
S_CASE <- 322448 / 428138; N_TOT <- 428138; res <- list()   # total N with the case fraction, as for FinnGen (500,264 with s = 154,630/500,264)
for (mode in c("iterative", "allbutone")) for (loc in c("TERC", "TERT")) for (w in c(100, 250)) for (ld in c("haplotype", "dosage")) for (se_src in c("P-derived", "allele-frequency model")) {
  x <- read.csv(sprintf("coloc_v4_inputs/%s_%dkb_HYPERTENSION.csv", loc, w)); x <- x[x$in_panel %in% c(TRUE, "True"), ]
  f <- if (ld == "haplotype") sprintf("coloc_v4_inputs/%s_%dkb_HYPERTENSION_LD.csv", loc, w) else sprintf("coloc_v6_inputs/%s_%dkb_HYPERTENSION_LDdos.csv", loc, w)
  LD <- as.matrix(read.csv(f, row.names = 1, check.names = FALSE)); stopifnot(all(rownames(LD) == x$varId)); dimnames(LD) <- list(x$varId, x$varId)
  stopifnot(all(!is.na(x$eur_af_alt))); maf <- setNames(pmin(x$eur_af_alt, 1 - x$eur_af_alt), x$varId)
  pq <- 2 * x$eur_af_alt * (1 - x$eur_af_alt); c0 <- median((x$se_out * sqrt(pq))[maf >= 0.05])
  se <- if (se_src == "P-derived") x$se_out else c0 / sqrt(pq)
  d1 <- list(beta = setNames(x$b_ltl, x$varId), varbeta = setNames(x$se_ltl^2, x$varId), snp = x$varId, position = x$pos, type = "quant", sdY = 1, N = 464716, MAF = maf, LD = LD)
  d2 <- list(beta = setNames(x$b_out, x$varId), varbeta = setNames(se^2, x$varId), snp = x$varId, position = x$pos, type = "cc", s = S_CASE, N = N_TOT, MAF = maf, LD = LD)
  r <- tryCatch(suppressMessages(suppressWarnings(coloc.signals(d1, d2, method = "cond", mode = mode, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5, maxhits = 5, pthr = 1e-6))), error = function(e) e)
  if (inherits(r, "error")) { res[[length(res) + 1]] <- data.frame(mode = mode, locus = loc, window_kb = w, ld = ld, se_source = se_src, hit1 = NA, hit2 = NA, hit1_rsid = NA, hit2_rsid = NA, nsnps = nrow(x), PP.H3 = NA, PP.H4 = NA, note = conditionMessage(r)); next }
  s <- as.data.frame(r$summary); rs <- setNames(x$rsid, x$varId)
  res[[length(res) + 1]] <- data.frame(mode = mode, locus = loc, window_kb = w, ld = ld, se_source = se_src, hit1 = s$hit1, hit2 = s$hit2, hit1_rsid = rs[s$hit1], hit2_rsid = rs[s$hit2],
                                       nsnps = s$nsnps, PP.H3 = s$PP.H3.abf, PP.H4 = s$PP.H4.abf, note = "")
  cat(mode, loc, w, ld, se_src, nrow(s), "pairs\n")
}
O <- do.call(rbind, res); rownames(O) <- NULL; write.csv(O, "step14_htn_conditional_coloc.csv", row.names = FALSE); print(O[, c(1:5, 8:9, 11:12)], digits = 3)
