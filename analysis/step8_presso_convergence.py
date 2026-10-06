"""v6: MR-PRESSO outlier-set and corrected-estimate convergence as the number of simulations increases (vectorized implementation,
step8_presso_np.py). B = 2,400; 10,000; 50,000; 100,000; 1,000,000; five fixed seeds each (1-5). Stability criterion (prespecified here):
the outlier set is identical across all five seeds at the given B."""
import numpy as np, pandas as pd
from step8_presso_np import presso
rows, sets = [], {}
for code in ['S_mvpEU', 'D_mvpEU', 'H_mvpEU']:
    d = pd.read_csv(f'mr_v4_inputs/forward_selected_r2_0.001_{code}.csv')
    for B in [2400, 10000, 50000, 100000, 1000000]:
        for seed in range(1, 6):
            r, p_out, out = presso(d.bx.values, d.sx.values, d.by.values, d.sy.values, B, seed)
            s = frozenset(d.rsid[out]); sets[(code, B, seed)] = s
            rows.append(dict(code=code, B=B, seed=seed, outlier_P_resolution=len(d) / B, **r, outliers=';'.join(sorted(s))))
R = pd.DataFrame(rows)
summ = []
for (code, B), g in R.groupby(['code', 'B']):
    ss = [sets[(code, B, s)] for s in range(1, 6)]; inter, uni = frozenset.intersection(*ss), frozenset.union(*ss)
    summ.append(dict(code=code, B=B, outliers_min=g.n_outliers.min(), outliers_max=g.n_outliers.max(), identical_across_5_seeds=len(set(ss)) == 1,
                     in_all=len(inter), in_any=len(uni), b_min=g.b.min(), b_max=g.b.max(), se_mean=g.se.mean(), global_p_max=g.global_p.max()))
S = pd.DataFrame(summ); R.to_csv('step8_presso_convergence_runs.csv', index=False); S.to_csv('step8_presso_convergence_summary.csv', index=False)
print(S.to_string())
