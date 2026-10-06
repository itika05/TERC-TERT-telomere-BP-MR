# v10 (review of v9, major 4): independent reconstruction of the contamination-mixture P value and CI from the profile likelihood of
# Burgess et al. 2020 (Nat Commun 11:376), written from the paper's definition, compared with MendelianRandomization::mr_conmix.
# For each variant j: ratio r_j = by/bx, s_j = sy/|bx|; valid component N(r_j; theta, s_j); invalid component N(r_j; 0, sqrt(psi^2 + s_j^2));
# the profile log-likelihood at theta takes, for each variant, the larger of the two densities. P = P(chi2_1 > 2[l(theta_hat) - l(0)]).
suppressMessages(library(MendelianRandomization))
out <- list()
for (code in c("S_mvpEU", "D_mvpEU", "H_mvpEU", "H_fg12")) {
  d <- read.csv(sprintf("mr_v4_inputs/forward_selected_r2_0.001_%s.csv", code)); r <- d$by / d$bx; s <- d$sy / abs(d$bx)
  psi <- 1.5 * sd(r)
  ll <- function(th) sum(log(pmax(dnorm(r, th, s), dnorm(r, 0, sqrt(psi^2 + s^2)))))
  grid <- seq(-1, 1, by = 0.0005); L <- sapply(grid, ll); th <- grid[which.max(L)]
  p_own <- pchisq(2 * (max(L) - ll(0)), 1, lower.tail = FALSE); inci <- grid[2 * (max(L) - L) < qchisq(0.95, 1)]
  cm <- mr_conmix(mr_input(bx = d$bx, bxse = d$sx, by = d$by, byse = d$sy), psi = 0, CIMin = -1, CIMax = 1, CIStep = 0.0005)
  out[[code]] <- data.frame(code = code, psi_own = psi, psi_pkg = cm$Psi, est_own = th, est_pkg = cm$Estimate, p_own = p_own, p_pkg = cm$Pvalue,
                            ci_own = sprintf("%.4f to %.4f", min(inci), max(inci)), ci_pkg = sprintf("%.4f to %.4f", min(cm$CILower), max(cm$CIUpper)),
                            loglik_ratio = 2 * (max(L) - ll(0)), wald_p_from_ci = 2 * pnorm(-abs(th) / ((max(inci) - min(inci)) / 3.92)))
}
O <- do.call(rbind, out); write.csv(O, "step12_conmix_check.csv", row.names = FALSE); print(O, digits = 4)
# Over-dispersion factor phi as computed inside mr_conmix (MendelianRandomization 0.10.0) for the valid set at the maximum, and the P value
# with the likelihood-ratio statistic divided by phi^2 (the direction used for the package CI) instead of multiplied.
ph <- list()
for (code in c("S_mvpEU", "D_mvpEU", "H_mvpEU", "H_fg12")) {
  d <- read.csv(sprintf("mr_v4_inputs/forward_selected_r2_0.001_%s.csv", code)); r <- d$by / d$bx; s <- d$sy / abs(d$bx); psi <- 1.5 * sd(r)
  grid <- seq(-1, 1, by = 0.0005); best <- -Inf; vb <- NULL
  for (th in grid) { li <- dnorm(r, th, s, log = TRUE); le <- dnorm(r, 0, sqrt(psi^2 + s^2), log = TRUE); v <- li > le; l <- sum(ifelse(v, li, le)); if (l > best) { best <- l; vb <- v } }
  w <- s[vb]^-2; phi <- if (sum(vb) < 1.5) 1 else max(sqrt(sum((r[vb] - weighted.mean(r[vb], w))^2 * w) / (sum(vb) - 1)), 1)
  li0 <- dnorm(r, 0, s, log = TRUE); le0 <- dnorm(r, 0, sqrt(psi^2 + s^2), log = TRUE); l0 <- sum(ifelse(li0 > le0, li0, le0)); LR <- 2 * (best - l0)
  ph[[code]] <- data.frame(code = code, n_valid = sum(vb), phi = phi, LR = LR, p_package_rule = pchisq(LR * phi^2, 1, lower.tail = FALSE), p_LR_divided_by_phi2 = pchisq(LR / phi^2, 1, lower.tail = FALSE), p_no_phi = pchisq(LR, 1, lower.tail = FALSE))
}
P <- do.call(rbind, ph); write.csv(P, "step12_conmix_phi.csv", row.names = FALSE); print(P, digits = 4)
