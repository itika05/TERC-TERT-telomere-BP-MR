"""v6 combined MR results: v4 results (step6v4_mr_results_R.csv) with (i) the primary reverse-MR rows replaced by the analysis using source LTL
statistics (codes kept as '<exposure>->L_codd'; source column updated), (ii) PLINK 1.9 clumped sets and the FinnGen ambiguous-palindrome sensitivity
appended, (iii) MR-PRESSO at B = 1,000,000 (validated re-implementation, seed 1; identical outlier sets across seeds 1-5) as the reported
MR-PRESSO row, with the MRPRESSO 1.0 package run at B = 2,400 retained under a separate label."""
import numpy as np, pandas as pd
from scipy import stats
R = pd.read_csv('step6v4_mr_results_R.csv'); N = pd.read_csv('step8_mr_v6_R.csv')
rev = R.direction.eq('reverse') & R.code.str.endswith('L_codd')
R = R[~rev].copy()
N['code'] = N.code.str.replace('->L_src', '->L_codd', regex=False)
R.loc[R.method.str.startswith('MR-PRESSO'), 'method'] = 'MR-PRESSO (MRPRESSO 1.0, B = 2,400)'
C = pd.read_csv('step8_presso_convergence_runs.csv'); C = C[(C.B == 1000000) & (C.seed == 1)]
meta = R[R.analysis.eq('r2<0.001 (primary)') & R.method.eq('IVW (MRE, t)')].set_index('code')
pr = []
for r in C.itertuples():
    m = meta.loc[r.code]; q = stats.t.ppf(0.975, r.k_retained - 1)
    pr.append(dict(analysis='r2<0.001 (primary)', method='MR-PRESSO (outlier-corrected)', b=r.b, se=r.se, lo=r.b - q * r.se, hi=r.b + q * r.se, p=r.p, nsnp=r.k_retained,
                   direction='forward', code=r.code, trait=m.trait, source=m.source, family=m.family, global_p=1e-6 if r.global_count == 0 else r.global_p,
                   n_outliers=r.n_outliers, presso_B=1000000, presso_global_count=r.global_count))
out = pd.concat([R, N, pd.DataFrame(pr)], ignore_index=True); out.to_csv('step8_mr_results_v6.csv', index=False)
print(out[out.method.str.startswith('MR-PRESSO')][['code', 'method', 'b', 'lo', 'hi', 'nsnp', 'n_outliers']].to_string())
print(out[out.direction.eq('reverse') & out.analysis.eq('r2<0.001 (primary)') & out.method.eq('IVW (MRE, t)')][['code', 'b', 'lo', 'hi', 'p', 'nsnp', 'source']].to_string())
