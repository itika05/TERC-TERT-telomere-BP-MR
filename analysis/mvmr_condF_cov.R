# v5: conditional F for LTL, BMI and lymphocyte count under an assumed correlation (rho) between the estimation errors of the exposure
# estimates (sample overlap: UK Biobank contributes to all three exposure GWAS). Sigma_l = D_l R D_l with R_jk = rho (j != k), D_l = diag(se).
# Sanderson, Spiller and Bowden (Stat Med 2021) eq. for Q_x: sum_l (bx_lj - bx_l,-j' delta)^2 / (s_lj^2 + delta' S_-j delta - 2 delta' s_j,-j).
d <- read.csv("mr_v4_inputs/mvmr_instruments.csv")
ex <- c("LTL", "BMI", "LYM"); bx <- as.matrix(d[, paste0(ex, "_b")]); sx <- as.matrix(d[, paste0(ex, "_se")])
condF <- function(j, rho) {
  L <- nrow(bx); k <- ncol(bx)
  obj <- function(delta) { sd <- sx[, -j, drop = FALSE] %*% diag(delta, nrow = k - 1)  # s_la * delta_a
    v <- sx[, j]^2 + rowSums(sd^2) + rho * (rowSums(sd)^2 - rowSums(sd^2)) - 2 * rho * sx[, j] * rowSums(sd)
    sum((bx[, j] - bx[, -j, drop = FALSE] %*% delta)^2 / v) }
  st <- list(rep(0, k - 1), coef(lm(bx[, j] ~ 0 + bx[, -j])), c(1, 1), c(-1, -1), c(1, -1), c(-1, 1))
  fits <- lapply(st, function(s0) optim(s0, obj, method = "BFGS", control = list(maxit = 2000, reltol = 1e-12)))
  best <- fits[[which.min(sapply(fits, `[[`, "value"))]]
  c(F = best$value / (L - (k - 1)), d1 = best$par[1], d2 = best$par[2], conv = best$convergence) }
out <- expand.grid(exposure = ex, rho = c(-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3), stringsAsFactors = FALSE)
res <- t(mapply(function(e, r) condF(match(e, ex), r), out$exposure, out$rho))
out$conditional_F <- res[, "F"]; out$delta1 <- res[, "d1"]; out$delta2 <- res[, "d2"]; out$optim_convergence <- res[, "conv"]
out$n_variants <- nrow(bx); write.csv(out, "step7_mvmr_condF_covariance_grid.csv", row.names = FALSE); print(out)
