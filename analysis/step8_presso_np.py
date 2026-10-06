"""v6-v7: vectorized re-implementation of the MR-PRESSO global and outlier tests (Verbanck et al. 2018; MRPRESSO 1.0 R code logic) to allow a large
number of simulations (B). Algorithm, as in MRPRESSO::mr_presso with weights 1/SdOutcome^2:
  observed: leave-one-out IVW (no intercept) estimate b_-i for each variant i; RSSobs = sum_i w_i (by_i - b_-i bx_i)^2;
  random data (per simulation): bx*_i ~ N(bx_i, sx_i); by*_i ~ N(b_-i bx_i, sy_i) (prediction of the LOO model at the observed bx_i);
  global test: RSS of the random data, recomputing leave-one-out estimates on the random data; P = #(RSSexp > RSSobs)/B;
  outlier test for variant i: Dif_i = by_i - b_-i bx_i; Exp_i = by*_i - b_-i bx*_i; P_i = #(Exp_i^2 > Dif_i^2)/B, multiplied by k (Bonferroni), capped at 1;
  outliers: P_i <= 0.05 (as in MRPRESSO), outlier test run only if the global P < 0.05; corrected estimate: weighted least squares without outliers (lm with weights), t(k_retained - 1) inference.
Equivalence with MRPRESSO 1.0 is tested on shared simulated arrays in step9_presso_equivalence.py (the distortion test is not re-implemented)."""
import numpy as np, pandas as pd, sys
from scipy import stats
def _stats(bx, sx, by, sy):
    w = 1 / sy ** 2; Sxy, Sxx = (w * bx * by).sum(), (w * bx ** 2).sum()
    bloo = (Sxy - w * bx * by) / (Sxx - w * bx ** 2)
    return w, bloo, (w * (by - bloo * bx) ** 2).sum(), (by - bloo * bx) ** 2
def _accumulate(X, Y, w, bloo, rss_obs, dif2):
    """X, Y: (m, k) random exposure and outcome data. Returns per-replicate RSSexp, global exceedances, per-variant outlier exceedances."""
    sxy = (w * X * Y).sum(1, keepdims=True); sxx = (w * X ** 2).sum(1, keepdims=True)
    bl = (sxy - w * X * Y) / (sxx - w * X ** 2)
    rss_exp = (w * (Y - bl * X) ** 2).sum(1)
    return rss_exp, int((rss_exp > rss_obs).sum()), ((Y - bloo * X) ** 2 > dif2).sum(0)
def _finish(bx, by, w, B, n_glob, n_out, alpha=0.05):
    k = len(bx); global_p = n_glob / B
    # as MRPRESSO::mr_presso: outlier test only when the global P < SignifThreshold; P_i = min(1, k * #/B); outliers P_i <= SignifThreshold
    p_out = np.minimum(1, n_out / B * k) if global_p < alpha else np.full(k, np.nan)
    out = (p_out <= alpha) if global_p < alpha else np.zeros(k, bool)
    keep = ~out; ww = w[keep]; b = (ww * bx[keep] * by[keep]).sum() / (ww * bx[keep] ** 2).sum()
    res = np.sqrt(ww) * (by[keep] - b * bx[keep]); df = keep.sum() - 1; s2 = (res ** 2).sum() / df
    se = np.sqrt(s2 / (ww * bx[keep] ** 2).sum()); p = 2 * stats.t.sf(abs(b / se), df)
    return dict(global_p=global_p, global_count=n_glob, n_outliers=int(out.sum()), b=b, se=se, p=p, k_retained=int(keep.sum())), p_out, out
def presso(bx, sx, by, sy, B, seed, chunk=20000):
    bx, sx, by, sy = map(np.asarray, (bx, sx, by, sy)); k = len(bx)
    w, bloo, rss_obs, dif2 = _stats(bx, sx, by, sy)
    rng = np.random.default_rng(seed); n_glob = 0; n_out = np.zeros(k); done = 0
    while done < B:
        m = min(chunk, B - done)
        X = rng.normal(bx, sx, size=(m, k)); Y = rng.normal(bloo * bx, sy, size=(m, k))
        _, g, o = _accumulate(X, Y, w, bloo, rss_obs, dif2); n_glob += g; n_out += o; done += m
    return _finish(bx, by, w, B, n_glob, n_out)
def presso_from_arrays(bx, sx, by, sy, X, Y):
    """Same statistics computed on supplied random data (used to test equivalence with MRPRESSO on shared simulated arrays)."""
    bx, sx, by, sy = map(np.asarray, (bx, sx, by, sy)); w, bloo, rss_obs, dif2 = _stats(bx, sx, by, sy)
    rss_exp, g, o = _accumulate(X, Y, w, bloo, rss_obs, dif2)
    r, p_out, out = _finish(bx, by, w, X.shape[0], g, o)
    return r, p_out, out, dict(bloo=bloo, rss_obs=rss_obs, dif2=dif2, rss_exp=rss_exp, n_out=o)
if __name__ == '__main__':
    code, B, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    d = pd.read_csv(f'mr_v4_inputs/forward_selected_r2_0.001_{code}.csv')
    r, p_out, out = presso(d.bx, d.sx, d.by, d.sy, B, seed)
    print(code, B, seed, r, ';'.join(d.rsid[out]))
