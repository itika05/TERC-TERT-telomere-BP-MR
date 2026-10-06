# v8: pointwise test of the qhet objective used for the bootstrap. The reference objective PL2_MVMR (the Q statistic with
# w = seY^2 + b' Sigma_l b + tau2) is taken verbatim from the body of MVMR::qhet_mvmr (ref_MVMR/qhet_mvmr.R, commit 8be0d94): the
# Qtemp function is extracted from the qhet_mvmr body and its statements up to the definition of covlist/gammahat/segamma/pihat are
# evaluated, then the second PL2_MVMR definition (the one using tau_i) is evaluated in the same environment. It is compared with the
# vectorised objective in mvmr_v7_qhet.R at random coefficient vectors and tau2 values for every outcome and rho.
suppressMessages(library(MendelianRandomization))
for (f in list.files("ref_MVMR", "\\.R$", full.names = TRUE)) source(f)
d0 <- read.csv("mr_v4_inputs/mvmr_instruments.csv")
fb <- Reduce(`|`, lapply(c("LTL", "BMI", "LYM", "S_mvpEU", "D_mvpEU", "H_mvpEU"), function(k) d0[[paste0(k, "_p")]] <= 1e-300 | d0[[paste0(k, "_p")]] >= 0.999 | d0[[paste0(k, "_b")]] == 0))
d <- d0[!fb, ]; ex <- c("LTL", "BMI", "LYM"); bx <- as.matrix(d[, paste0(ex, "_b")]); sx <- as.matrix(d[, paste0(ex, "_se")])
qb <- as.list(body(qhet_mvmr)); i_q <- which(sapply(qb, function(e) grepl("^Qtemp <- function", paste(deparse(e), collapse = " "))))
Qtemp <- eval(qb[[i_q]][[3]])                        # the function object as written in the reference code
tb <- as.list(body(Qtemp))
i_pl2 <- which(sapply(tb, function(e) { t <- paste(deparse(e), collapse = " "); grepl("^PL2_MVMR = function", t) && grepl("tau_i", t) }))
i_setup <- which(sapply(tb, function(e) grepl("^pihat <-", paste(deparse(e), collapse = " "))))
stopifnot(length(i_pl2) == 1, length(i_setup) == 1)
wq <- function(b, tau2, sy, sxm, pcor) { sb <- sweep(sxm, 2, b, `*`); sy^2 + rowSums((sb %*% pcor) * sb) + tau2 }
fast_obj <- function(b, tau2, gy, sy, px, sxm, pcor) sum((gy - px %*% b)^2 / wq(b, tau2, sy, sxm, pcor))
set.seed(20260930); out <- list()
for (oc in c("S_mvpEU", "D_mvpEU", "H_mvpEU")) {
  fm <- format_mvmr(BXGs = bx, BYG = d[[paste0(oc, "_b")]], seBXGs = sx, seBYG = d[[paste0(oc, "_se")]], RSID = d$varId)
  for (rho in c(0, 0.1, 0.2, 0.3)) {
    R <- matrix(rho, 3, 3); diag(R) <- 1
    env <- new.env(); env$r_input <- fm; env$pcor <- R
    for (k in 2:i_setup) eval(tb[[k]], env)           # exp.number, stderr, correlation, covlist, gammahat, segamma, pihat
    eval(tb[[i_pl2]], env)                            # PL2_MVMR as written (uses tau_i)
    for (j in 1:200) {
      b <- rnorm(3, 0, 0.5); tau2 <- runif(1, 0, 2e-5); env$tau_i <- tau2
      ref <- env$PL2_MVMR(b); fst <- fast_obj(b, tau2, fm$betaYG, fm$sebetaYG, bx, sx, R)
      out[[length(out) + 1]] <- data.frame(outcome = oc, rho = rho, point = j, b1 = b[1], b2 = b[2], b3 = b[3], tau2 = tau2, reference = ref, vectorised = fst, rel_diff = abs(ref - fst) / abs(ref))
    }
  }
}
O <- do.call(rbind, out)
Of <- O; for (cn in c("b1", "b2", "b3", "tau2", "reference", "vectorised", "rel_diff")) Of[[cn]] <- sprintf("%.17g", O[[cn]])   # full double precision so rel_diff can be re-derived
write.csv(Of, "step10_qhet_objective_test.csv", row.names = FALSE)
cat("points", nrow(O), "max rel diff", max(O$rel_diff), "\n"); stopifnot(max(O$rel_diff) < 1e-10)
