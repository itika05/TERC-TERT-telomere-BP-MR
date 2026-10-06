# v9 (review M6): conditional colocalization. coloc 6.0.3 coloc.signals(method = "cond", mode = mode) performs stepwise
# conditional analysis of each trait on the reference LD matrix (signals added while the conditional P < pthr), then runs coloc.abf on
# every pair of conditioned signals. Inputs and LD are those of the v4/v6 regional analyses (1000 Genomes phase 3, 503 European
# individuals), with phased-haplotype and genotype-dosage LD. Settings: p1 = p2 = 1e-4, p12 = 1e-5, pthr = 1e-6 (coloc default),
# maxhits = 5. LTL N = 464,716; outcome N = median per-variant N; quantitative traits sdY = 1; hypertension type "cc".
# A UK Biobank-derived LD reference (e.g., released UK Biobank LD matrices) is not held locally and could not be downloaded in the
# analysis environment, so these analyses use the same 1000 Genomes reference as the SuSiE analyses.
suppressMessages(library(coloc))
TR <- c("SBP", "DBP", "HYPERTENSION", "FG_HYPERTENSION"); res <- list(); hits <- list()
# MVP hypertension is not analysed: conditioning a case-control trait needs the case fraction s, and the MVP European HYPERTENSION record
# gives only an effective N (case and control counts were not available). FinnGen I9_HYPTENS: s = 154,630 / 500,264.
for (mode in c("iterative", "allbutone")) for (loc in c("TERC", "TERT")) for (w in c(100, 250)) for (t in setdiff(TR, "HYPERTENSION")) for (ld in c("haplotype", "dosage")) {
  x <- read.csv(sprintf("coloc_v4_inputs/%s_%dkb_%s.csv", loc, w, t)); x <- x[x$in_panel %in% c(TRUE, "True"), ]
  f <- if (ld == "haplotype") sprintf("coloc_v4_inputs/%s_%dkb_%s_LD.csv", loc, w, t) else sprintf("coloc_v6_inputs/%s_%dkb_%s_LDdos.csv", loc, w, t)
  LD <- as.matrix(read.csv(f, row.names = 1, check.names = FALSE)); stopifnot(all(rownames(LD) == x$varId)); dimnames(LD) <- list(x$varId, x$varId)
  stopifnot(all(!is.na(x$eur_af_alt))); maf <- setNames(pmin(x$eur_af_alt, 1 - x$eur_af_alt), x$varId)   # 1000 Genomes EUR (reference panel) MAF
  d1 <- list(beta = setNames(x$b_ltl, x$varId), varbeta = setNames(x$se_ltl^2, x$varId), snp = x$varId, position = x$pos, type = "quant", sdY = 1, N = 464716, MAF = maf, LD = LD)
  d2 <- if (t == "FG_HYPERTENSION") list(beta = setNames(x$b_out, x$varId), varbeta = setNames(x$se_out^2, x$varId), snp = x$varId, position = x$pos, type = "cc", s = 154630 / 500264, N = median(x$n_out), MAF = maf, LD = LD)
        else list(beta = setNames(x$b_out, x$varId), varbeta = setNames(x$se_out^2, x$varId), snp = x$varId, position = x$pos, type = "quant", sdY = 1, N = median(x$n_out), MAF = maf, LD = LD)
  r <- tryCatch(suppressMessages(suppressWarnings(coloc.signals(d1, d2, method = "cond", mode = mode, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5, maxhits = 5, pthr = 1e-6))),
                error = function(e) e)
  if (inherits(r, "error")) { res[[length(res) + 1]] <- data.frame(mode = mode, locus = loc, window_kb = w, trait = t, ld = ld, hit1 = NA, hit2 = NA, hit1_rsid = NA, hit2_rsid = NA, nsnps = nrow(x), PP.H3 = NA, PP.H4 = NA, note = conditionMessage(r)); next }
  s <- as.data.frame(r$summary)
  rs <- setNames(x$rsid, x$varId)
  res[[length(res) + 1]] <- data.frame(mode = mode, locus = loc, window_kb = w, trait = t, ld = ld, hit1 = s$hit1, hit2 = s$hit2, hit1_rsid = rs[s$hit1], hit2_rsid = rs[s$hit2], nsnps = s$nsnps,
                                       PP.H3 = s$PP.H3.abf, PP.H4 = s$PP.H4.abf, note = "")
  cat(mode, loc, w, t, ld, nrow(s), "pairs\n")
}
O <- do.call(rbind, res); rownames(O) <- NULL; write.csv(O, "step11_conditional_coloc.csv", row.names = FALSE); print(O[, c(1:5, 8:9, 11:12)], digits = 3)
