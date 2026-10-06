# v7: residual standard error of MR-Egger (MendelianRandomization 0.10.0) for the result rows whose stored interval follows the package's
# residual-SE < 1 branch (requested by step9_output_consistency.py). Refits mr_egger on the named input and records RSE and outputs.
suppressMessages(library(MendelianRandomization))
req <- read.csv("step9_egger_rse_request.csv", stringsAsFactors = FALSE); out <- list()
S13 <- read.csv("step9_reverse_inputs_S13.csv")
for (kk in req$analysis_key) {
  f <- trimws(strsplit(kk, "\\|")[[1]]); dirn <- f[1]; code <- f[2]; an <- f[3]
  if (dirn == "forward" && an == "r2<0.001 (primary)") d <- read.csv(sprintf("mr_v4_inputs/forward_selected_r2_0.001_%s.csv", code))
  else if (dirn == "reverse" && an == "r2<0.001 (primary)") d <- S13[S13$code == code & S13$in_custom_primary == "True" | S13$code == code & S13$in_custom_primary == TRUE, ]
  else stop(paste("no input mapping for", kk))
  d <- d[complete.cases(d[, c("bx", "sx", "by", "sy")]), ]
  e <- mr_egger(mr_input(bx = d$bx, bxse = d$sx, by = d$by, byse = d$sy), distribution = "t-dist")
  out[[length(out) + 1]] <- data.frame(analysis_key = kk, egger_k_refit = nrow(d), egger_RSE = e@RSE, egger_b_refit = e@Estimate, egger_lo_refit = e@CILower.Est, egger_hi_refit = e@CIUpper.Est)
}
write.csv(do.call(rbind, out), "step9_egger_rse.csv", row.names = FALSE)
