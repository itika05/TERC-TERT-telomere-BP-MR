"""v5: MR-PRESSO outlier stability across three seeds (314159 = v4 primary, 20261, 20262), NbDistribution = 2400 each.
Also: fixed-set IVW (multiplicative random effects) after removing (i) variants flagged in all three seeds and (ii) variants flagged in any seed."""
import pandas as pd, numpy as np, glob
from scipy import stats
R = pd.concat([pd.read_csv(f) for f in glob.glob('presso_stab/*_2400_*.csv')]).sort_values(['code', 'seed'])
def ivw(d):
    w = 1 / d.sy**2; b = (w * d.bx * d.by).sum() / (w * d.bx**2).sum(); k = len(d)
    res = (d.by - b * d.bx) / d.sy; phi = max(1, (res**2).sum() / (k - 1)); se = np.sqrt(phi / (w * d.bx**2).sum())
    t = stats.t.ppf(0.975, k - 1); return b, se, b - t * se, b + t * se, 2 * stats.t.sf(abs(b / se), k - 1), k
rows = []; seeds = []
for code, g in R.groupby('code'):
    d = pd.read_csv(f'mr_v4_inputs/forward_selected_r2_0.001_{code}.csv')
    sets = [set(s.split(';')) for s in g.outliers]
    allk, anyk = set.intersection(*sets), set.union(*sets)
    for _, r in g.iterrows():
        seeds.append(dict(code=code, seed=r.seed, n_outliers=r.n_outliers, b=r.b, se=r.se, p=r.p, global_p=r.global_p, minutes=round(r.minutes, 1),
                          jaccard_vs_v4=len(set(r.outliers.split(';')) & sets[list(g.seed).index(314159)]) / len(set(r.outliers.split(';')) | sets[list(g.seed).index(314159)])))
    for lab, ex in [('all instruments', set()), ('excl. outliers flagged in all 3 seeds', allk), ('excl. outliers flagged in any seed', anyk)]:
        b, se, lo, hi, p, k = ivw(d[~d.rsid.isin(ex)])
        rows.append(dict(code=code, set=lab, n_removed=len(ex), k=k, b=b, se=se, lo=lo, hi=hi, p=p, removed=';'.join(sorted(ex))))
pd.DataFrame(seeds).to_csv('step7_presso_seed_runs.csv', index=False)
pd.DataFrame(rows).to_csv('step7_presso_consensus_ivw.csv', index=False)
print(pd.DataFrame(seeds).to_string()); print(pd.DataFrame(rows).drop(columns='removed').to_string())
