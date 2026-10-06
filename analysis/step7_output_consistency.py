"""v5: internal consistency of stored MR outputs: recompute P and 95% CI from b and SE with the documented reference distribution
(IVW and IVW-GLS: t(k-1); MR-Egger: t(k-2); weighted median / mode: normal; MR-PRESSO: as reported by MRPRESSO, t(k_retained-1)) and compare."""
import numpy as np, pandas as pd
from scipy import stats
R = pd.read_csv('step6v4_mr_results_R.csv'); out = []
for r in R.itertuples():
    k = r.nsnp
    if r.method.startswith('MR-Egger'): df = k - 2
    elif r.method.startswith('IVW'): df = k - 1
    elif r.method.startswith('MR-PRESSO'): df = k - 1
    else: df = None
    z = r.b / r.se
    p = 2 * (stats.t.sf(abs(z), df) if df else stats.norm.sf(abs(z))); q = stats.t.ppf(0.975, df) if df else 1.959964
    out.append(dict(code=r.code, analysis=r.analysis, method=r.method, k=k, ref_dist=f't({df})' if df else 'normal', p_reported=r.p, p_recomputed=p,
                    log10_ratio=abs(np.log10(max(r.p, 1e-300)) - np.log10(max(p, 1e-300))), lo_reported=r.lo, lo_recomputed=r.b - q * r.se, hi_reported=r.hi, hi_recomputed=r.b + q * r.se))
O = pd.DataFrame(out); O['ci_abs_diff'] = np.maximum(abs(O.lo_reported - O.lo_recomputed), abs(O.hi_reported - O.hi_recomputed))
O.to_csv('step7_mr_output_consistency.csv', index=False)
print(O.groupby('method').agg(n=('k', 'size'), max_log10_P_diff=('log10_ratio', 'max'), max_ci_diff=('ci_abs_diff', 'max')).to_string())
print(O.sort_values('log10_ratio', ascending=False).head(6)[['code', 'analysis', 'method', 'p_reported', 'p_recomputed', 'ci_abs_diff']].to_string())
