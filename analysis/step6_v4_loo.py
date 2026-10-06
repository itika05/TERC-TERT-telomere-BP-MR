"""Leave-one-out (IVW, multiplicative random effects; identical to MendelianRandomization mr_ivw, see validation/mr_python_vs_R.csv),
per-instrument data for figures, and minimum detectable effects for reverse MR, on the v4 inputs."""
import pandas as pd, numpy as np
from scipy.stats import t as tdist, norm
from mr_lib import ivw
loo, mdes = [], []
for code in ['S_mvpEU', 'D_mvpEU', 'H_mvpEU', 'H_fg12']:
    d = pd.read_csv(f'mr_v4_inputs/forward_selected_r2_0.001_{code}.csv'); a = [d[c].values for c in ['bx', 'sx', 'by', 'sy']]
    for i in range(len(d)):
        m = np.ones(len(d), bool); m[i] = False; e = ivw(*[x[m] for x in a]); loo.append(dict(code=code, left_out=d.rsid.iloc[i], varId=d.varId.iloc[i], b=e['b'], se=e['se'], p=e['p']))
    pd.DataFrame(dict(code=code, varId=d.varId, rsid=d.rsid, gene=d.gene, F=(d.bx / d.sx) ** 2, bx=d.bx, sx=d.sx, by=d.by, sy=d.sy, steiger_z=d.steiger_z,
                      terc_tert=d.terc_tert_region)).to_csv(f'mr_v4_inputs/snp_level_{code}.csv', index=False)
pd.DataFrame(loo).to_csv('step6_v4_forward_leave_one_out.csv', index=False)
R = pd.read_csv('step6v4_mr_results_R.csv')
for r in R[(R.direction == 'reverse') & (R.analysis == 'r2<0.001 (primary)') & (R.method == 'IVW (MRE, t)')].itertuples():
    for al in [0.05, 0.05 / 6]:
        mdes.append(dict(code=r.code, alpha=al, k=r.nsnp, se=r.se, mde=(tdist.ppf(1 - al / 2, r.nsnp - 1) + norm.ppf(0.8)) * r.se))
pd.DataFrame(mdes).to_csv('step6_v4_reverse_mde.csv', index=False)
L = pd.DataFrame(loo); print(L.groupby('code').b.agg(['min', 'max']).round(4)); print(pd.DataFrame(mdes).round(4))
