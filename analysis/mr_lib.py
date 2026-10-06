"""Two-sample MR and colocalisation estimators (Python re-implementations following
TwoSampleMR / MendelianRandomization / MR-PRESSO / coloc conventions).
All functions take harmonised per-variant arrays: bx, sx (exposure), by, sy (outcome)."""
import numpy as np
from scipy import stats
from scipy.stats import norm


def ivw(bx, sx, by, sy):
    """IVW with multiplicative random effects (TwoSampleMR 'mr_ivw' default: se scaled by max(1, phi))."""
    w = bx ** 2 / sy ** 2
    b = np.sum(bx * by / sy ** 2) / np.sum(w)
    se_fe = 1 / np.sqrt(np.sum(w))
    resid = (by - b * bx) / sy
    k = len(bx)
    Q = np.sum(resid ** 2)
    phi = Q / (k - 1) if k > 1 else 1
    se = se_fe * max(1, np.sqrt(phi))
    p = 2 * stats.t.sf(abs(b / se), k - 1) if k > 1 else 2 * norm.sf(abs(b / se))
    return dict(method='IVW (MRE)', b=b, se=se, p=p, nsnp=k, Q=Q, Q_df=k - 1, Q_p=stats.chi2.sf(Q, k - 1),
                I2=max(0, (Q - (k - 1)) / Q) * 100 if Q > 0 else 0)


def egger(bx, sx, by, sy):
    """MR-Egger: orient to positive bx, weighted regression by ~ bx with intercept, residual SE scaling."""
    s = np.sign(bx); s[s == 0] = 1
    bx, by = bx * s, by * s
    w = 1 / sy ** 2
    X = np.column_stack([np.ones_like(bx), bx])
    W = np.diag(w)
    XtWX = X.T @ W @ X
    beta = np.linalg.solve(XtWX, X.T @ W @ by)
    resid = by - X @ beta
    k = len(bx)
    sigma2 = np.sum(w * resid ** 2) / (k - 2)
    cov = np.linalg.inv(XtWX) * max(1, sigma2)
    se = np.sqrt(np.diag(cov))
    t = beta / se
    p = 2 * stats.t.sf(np.abs(t), k - 2)
    return dict(method='MR-Egger', b=beta[1], se=se[1], p=p[1], nsnp=k,
                intercept=beta[0], intercept_se=se[0], intercept_p=p[0])


def _wmedian(b, w):
    o = np.argsort(b); b, w = b[o], w[o]
    ws = np.cumsum(w) - 0.5 * w
    ws = ws / np.sum(w)
    below = np.max(np.where(ws < 0.5)[0]) if np.any(ws < 0.5) else 0
    if below >= len(b) - 1:
        return b[-1]
    return b[below] + (b[below + 1] - b[below]) * (0.5 - ws[below]) / (ws[below + 1] - ws[below])


def weighted_median(bx, sx, by, sy, nboot=1000, seed=1):
    br = by / bx
    se_r = np.sqrt(sy ** 2 / bx ** 2 + by ** 2 * sx ** 2 / bx ** 4)  # second-order weights
    w = 1 / se_r ** 2
    b = _wmedian(br, w)
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(nboot):
        bxb = rng.normal(bx, sx); byb = rng.normal(by, sy)
        boots.append(_wmedian(byb / bxb, w))
    se = np.std(boots, ddof=1)
    return dict(method='Weighted median', b=b, se=se, p=2 * norm.sf(abs(b / se)), nsnp=len(bx))


def _mode(br, w, phi=1):
    # Hartwig 2017 MBE: normal kernel density, bandwidth via modified Silverman
    s = 0.9 * min(np.std(br, ddof=1), stats.median_abs_deviation(br, scale='normal')) / len(br) ** 0.2
    h = max(1e-8, s * phi)
    grid = np.linspace(br.min() - 3 * h, br.max() + 3 * h, 1501)
    dens = np.sum(w[None, :] * norm.pdf((grid[:, None] - br[None, :]) / h), axis=1)
    return grid[np.argmax(dens)]


def weighted_mode(bx, sx, by, sy, nboot=1000, seed=1):
    br = by / bx
    se_r = np.sqrt(sy ** 2 / bx ** 2 + by ** 2 * sx ** 2 / bx ** 4)
    w = se_r ** -2 / np.sum(se_r ** -2)
    b = _mode(br, w)
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(nboot):
        brb = rng.normal(br, se_r)
        boots.append(_mode(brb, w))
    se = np.std(boots, ddof=1)
    return dict(method='Weighted mode', b=b, se=se, p=2 * norm.sf(abs(b / se)), nsnp=len(bx))


def mr_presso(bx, sx, by, sy, nboot=10000, seed=1, alpha=0.05):
    """MR-PRESSO (Verbanck 2018), vectorised: global test on leave-one-out residuals, per-variant outlier
    test (Bonferroni x k), outlier-corrected IVW (weighted regression through origin, residual SE)."""
    rng = np.random.default_rng(seed)
    k = len(bx); w = 1 / sy ** 2

    def loo(bx, by):
        # bx, by can be (nboot,k)
        Sxy = np.sum(w * bx * by, axis=-1, keepdims=True); Sxx = np.sum(w * bx ** 2, axis=-1, keepdims=True)
        b_loo = (Sxy - w * bx * by) / (Sxx - w * bx ** 2)
        return b_loo

    b_loo = loo(bx[None, :], by[None, :])[0]
    rss_obs = w * (by - b_loo * bx) ** 2
    RSSobs = rss_obs.sum()
    exc = np.zeros(k); gl = 0
    chunk = 1000
    for c in range(0, nboot, chunk):
        n = min(chunk, nboot - c)
        bxs = rng.normal(bx, sx, size=(n, k))
        bys = rng.normal(b_loo * bx, sy, size=(n, k))
        bl = loo(bxs, bys)
        r = w * (bys - bl * bxs) ** 2
        gl += np.sum(r.sum(1) >= RSSobs)
        exc += np.sum(r >= rss_obs[None, :], axis=0)
    global_p = (gl + 1) / (nboot + 1)
    out_p = np.minimum(1, (exc + 1) / (nboot + 1) * k)
    outliers = np.where(out_p < alpha)[0] if global_p < alpha else np.array([], int)
    m = np.ones(k, bool); m[outliers] = False
    X = bx[m]; Y = by[m]; W = w[m]
    b = np.sum(W * X * Y) / np.sum(W * X ** 2)
    df = m.sum() - 1
    se = np.sqrt(np.sum(W * (Y - b * X) ** 2) / df / np.sum(W * X ** 2))
    return dict(method='MR-PRESSO (outlier-corrected)', b=b, se=se, p=2 * stats.t.sf(abs(b / se), df), nsnp=int(m.sum()),
                global_p=global_p, n_outliers=len(outliers), outlier_idx=outliers.tolist())


def fstat(bx, sx):
    return (bx / sx) ** 2


def steiger_r2_cont(b, se, n):
    """Variance explained by a variant for a standardised continuous trait from beta, se, n (TwoSampleMR get_r_from_bsen)."""
    t = b / se
    r2 = t ** 2 / (t ** 2 + n - 2)
    return r2


def steiger_r2_binary(logor, af, prev, ncase, ncontrol):
    """Variance explained on the logistic latent scale: 2f(1-f)b^2 / (2f(1-f)b^2 + pi^2/3). No prevalence or ascertainment
    correction is applied (prev, ncase and ncontrol are accepted for interface compatibility and ignored); this differs from
    TwoSampleMR get_r_from_lor, which uses prevalence."""
    lor = np.asarray(logor)
    af = np.asarray(af)
    val = (np.pi ** 2) / 3
    r2 = lor ** 2 * af * (1 - af) * 2 / (lor ** 2 * af * (1 - af) * 2 + val)
    return r2


def steiger_direction(r2x, nx, r2y, ny):
    """Per-variant Steiger z (Fisher z difference); positive => exposure->outcome direction."""
    zx = np.arctanh(np.sqrt(r2x)); zy = np.arctanh(np.sqrt(r2y))
    z = (zx - zy) / np.sqrt(1 / (nx - 3) + 1 / (ny - 3))
    return z, norm.sf(z)


# ---------------- coloc.abf (Giambartolomei 2014; Wakefield ABF) ----------------

def _labf(b, se, sd_prior):
    V = se ** 2
    r = sd_prior ** 2 / (sd_prior ** 2 + V)
    z = b / se
    return 0.5 * (np.log(1 - r) + r * z ** 2)


def _logsum(x):
    m = np.max(x)
    return m + np.log(np.sum(np.exp(x - m)))


def _logdiff(a, b):
    m = max(a, b)
    return m + np.log(np.exp(a - m) - np.exp(b - m))


def coloc_abf(b1, se1, b2, se2, sd1=0.15, sd2=0.15, p1=1e-4, p2=1e-4, p12=1e-5):
    """Default priors as coloc: sd prior 0.15 for quantitative (standardised) traits, 0.2 for case-control log OR."""
    l1 = _labf(b1, se1, sd1); l2 = _labf(b2, se2, sd2)
    lsum = l1 + l2
    lH0 = 0
    lH1 = np.log(p1) + _logsum(l1)
    lH2 = np.log(p2) + _logsum(l2)
    lH3 = np.log(p1) + np.log(p2) + _logdiff(_logsum(l1) + _logsum(l2), _logsum(lsum))
    lH4 = np.log(p12) + _logsum(lsum)
    all_ = np.array([lH0, lH1, lH2, lH3, lH4])
    pp = np.exp(all_ - _logsum(all_))
    snp_pp4 = np.exp(lsum - _logsum(lsum))
    return dict(nsnps=len(b1), PP0=pp[0], PP1=pp[1], PP2=pp[2], PP3=pp[3], PP4=pp[4]), snp_pp4
