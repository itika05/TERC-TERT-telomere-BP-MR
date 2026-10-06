"""Compare the purpose-written Python MR estimators (mr_lib) with MendelianRandomization 0.10.0 / MRPRESSO 1.0 on identical v3 inputs.
Differences are expressed in units of the R standard error. Known sources of difference are recorded in the 'explanation' column."""
import pandas as pd, numpy as np
R = pd.read_csv('step3v3_mr_results_R.csv'); P = pd.read_csv('step3v3_mr_results_python.csv')
R = R[R.analysis == 'r2<0.001 (primary)']
mapR = {'IVW (MRE, t)': 'IVW (MRE)', 'MR-Egger (t)': 'MR-Egger', 'Weighted median': 'Weighted median', 'Weighted mode (MBE)': 'Weighted mode'}
rows = []
for r in R.itertuples():
    m = 'MR-PRESSO (outlier-corrected)' if r.method.startswith('MR-PRESSO') else mapR.get(r.method)
    q = P[(P.code == r.code) & (P.method == m)]
    if q.empty or m is None: continue
    q = q.iloc[0]
    expl = {'IVW (MRE)': 'same estimator and inference', 'MR-Egger': 'same estimator; R uses t(k-2) and residual SE floor of 1',
            'Weighted median': 'Python uses second-order (delta) weights, MendelianRandomization default uses first-order weights; bootstrap SE differs by random stream', 'Weighted mode': 'R MBE uses modified Silverman bandwidth with phi = 1 and delta SE; Python bandwidth rule differs',
            'MR-PRESSO (outlier-corrected)': 'simulation-based; R NbDistribution 2400, Python 10000; outlier sets can differ'}[m]
    rows.append(dict(direction=r.direction, code=r.code, method=m, nsnp_R=r.nsnp, nsnp_py=q.nsnp, b_R=r.b, b_py=q.b, se_R=r.se, se_py=q.se, p_R=r.p, p_py=q.p,
                     diff_in_R_SE=(q.b - r.b) / r.se, se_ratio=q.se / r.se, explanation=expl))
V = pd.DataFrame(rows); V.to_csv('validation/mr_python_vs_R.csv', index=False)
print(V.groupby('method').agg(n=('code', 'size'), max_abs_diff_SE=('diff_in_R_SE', lambda x: x.abs().max()), median_se_ratio=('se_ratio', 'median')).round(3))
