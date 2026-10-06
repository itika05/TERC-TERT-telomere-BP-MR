"""Python estimators (mr_lib) on the same v3 inputs, for comparison with the R packages, plus leave-one-out, F statistics and
minimum detectable effects at alpha = 0.05 and alpha = 0.05/6."""
import pandas as pd, numpy as np
from scipy.stats import t as tdist, norm
from mr_lib import ivw, egger, weighted_median, weighted_mode, mr_presso
rows = []; loo = []
def mde(se, k, alpha): return (tdist.ppf(1 - alpha / 2, k - 1) + norm.ppf(0.8)) * se
for code in ['S_mvpEU', 'D_mvpEU', 'H_mvpEU', 'H_fg', 'S_gh', 'D_gh', 'S_bbj', 'D_bbj', 'S_mvpAA', 'D_mvpAA', 'H_mvpAA', 'S_mvpHS', 'D_mvpHS', 'H_mvpHS']:
    d = pd.read_csv(f'mr_v3_inputs/forward_selected_r2_0.001_{code}.csv'); a = [d[c].values for c in ['bx', 'sx', 'by', 'sy']]
    ests = [ivw(*a), egger(*a), weighted_median(*a), weighted_mode(*a)] + ([mr_presso(*a)] if code in ('S_mvpEU', 'D_mvpEU', 'H_mvpEU') else [])
    for e in ests: rows.append(dict(direction='forward', code=code, **{k: v for k, v in e.items() if k != 'outlier_idx'}))
    if code in ('S_mvpEU', 'D_mvpEU', 'H_mvpEU'):
        for i in range(len(d)):
            m = np.ones(len(d), bool); m[i] = False; e = ivw(*[x[m] for x in a]); loo.append(dict(code=code, left_out=d.rsid.iloc[i], varId=d.varId.iloc[i], b=e['b'], se=e['se'], p=e['p']))
        pd.DataFrame(dict(code=code, varId=d.varId, rsid=d.rsid, gene=d.gene, F=(d.bx / d.sx) ** 2, bx=d.bx, sx=d.sx, by=d.by, sy=d.sy, steiger_z=d.steiger_z)).to_csv(f'mr_v3_inputs/snp_level_{code}.csv', index=False)
for ex in ['SBP', 'DBP', 'Hypertension']:
    for oc in ['L_codd', 'L_nkSA']:
        d = pd.read_csv(f'mr_v3_inputs/reverse_{ex}_{oc}.csv'); a = [d[c].values for c in ['bx', 'sx', 'by', 'sy']]
        e = ivw(*a); k = len(d)
        rows.append(dict(direction='reverse', code=f'{ex}->{oc}', **{kk: v for kk, v in e.items() if kk != 'outlier_idx'}, mde_alpha05=mde(e['se'], k, 0.05), mde_alpha_bonf=mde(e['se'], k, 0.05 / 6)))
        for f in [egger, weighted_median, weighted_mode] + ([mr_presso] if oc == 'L_codd' else []):
            e = f(*a); rows.append(dict(direction='reverse', code=f'{ex}->{oc}', **{kk: v for kk, v in e.items() if kk != 'outlier_idx'}))
pd.DataFrame(rows).to_csv('step3v3_mr_results_python.csv', index=False); pd.DataFrame(loo).to_csv('step3v3_forward_leave_one_out.csv', index=False)
print(pd.DataFrame(rows)[['direction', 'code', 'method', 'nsnp', 'b', 'se', 'p', 'I2', 'intercept_p', 'n_outliers', 'mde_alpha05', 'mde_alpha_bonf']].round(4).to_string())
