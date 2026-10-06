"""v7 (Supplementary Table S46): internal consistency of EVERY row of the final MR result table (step8_mr_results_v6.csv).
Each row is identified by an explicit analysis key (direction | code | analysis | method | source); duplicate keys fail the build.
P and 95% CI are recomputed from the stored b and SE under the documented reference distribution and compared before rounding:
  IVW, IVW-GLS, MR-PRESSO: t(k - 1) (k = variants in the fit; for outlier-corrected MR-PRESSO, retained variants);
  MR-Egger (MendelianRandomization, distribution = "t-dist"): t(k - 2), EXCEPT that when the residual standard error is < 1 the package
  reports the wider of the normal interval and the t interval with SE x RSE, and P = max(normal P, t P with z / RSE) (egger.bounds and
  mr_egger source). Egger rows not reproduced by plain t(k - 2) are refitted in step9_egger_rse.R to obtain the RSE, and the package rule is applied;
  weighted median and mode: normal.
Tolerance 1e-6 absolute on CI bounds and 1e-6 relative on P. A row reproduced by no documented rule FAILS the build."""
import sys, subprocess
import numpy as np, pandas as pd
from scipy import stats
R = pd.read_csv('step8_mr_results_v6.csv')
key = R.direction.astype(str) + ' | ' + R.code + ' | ' + R.analysis + ' | ' + R.method + ' | ' + R.source.astype(str)
assert not key.duplicated().any(), 'duplicate analysis keys'
TOL = 1e-6; out = []
for (kk, r) in zip(key, R.itertuples()):
    k = int(r.nsnp); m = r.method
    if m.startswith('MR-Egger'): df = k - 2
    elif m.startswith('IVW') or m.startswith('MR-PRESSO'): df = k - 1
    else: df = None
    z = r.b / r.se
    cand = {}
    if df:
        q = stats.t.ppf(0.975, df); cand[f't({df})'] = (2 * stats.t.sf(abs(z), df), r.b - q * r.se, r.b + q * r.se)
    cand['normal'] = (2 * stats.norm.sf(abs(z)), r.b - 1.959963984540054 * r.se, r.b + 1.959963984540054 * r.se)
    expected = f't({df})' if df else 'normal'
    p, lo, hi = cand[expected]
    ok = abs(lo - r.lo) <= TOL and abs(hi - r.hi) <= TOL and abs(p - r.p) <= TOL * max(r.p, 1e-300) + 1e-300
    out.append(dict(analysis_key=kk, k=k, rule_expected=expected, rule_reproducing=expected if ok else None, p_stored=r.p, lo_stored=r.lo, hi_stored=r.hi,
                    p_recomputed=p, lo_recomputed=lo, hi_recomputed=hi, status='PASS' if ok else ('CHECK_EGGER_RSE' if m.startswith('MR-Egger') else 'FAIL'), b=r.b, se=r.se, df=df))
O = pd.DataFrame(out)
eg = O[O.status == 'CHECK_EGGER_RSE']
if len(eg):
    eg[['analysis_key']].to_csv('step9_egger_rse_request.csv', index=False)
    subprocess.run(['Rscript', 'step9_egger_rse.R'], check=True)
    E = pd.read_csv('step9_egger_rse.csv'); O = O.merge(E, on='analysis_key', how='left')
    for i in O.index[O.status == 'CHECK_EGGER_RSE']:
        r = O.loc[i]; qn, qt = 1.959963984540054, stats.t.ppf(0.975, r.df); rse = r.egger_RSE
        lo = min(r.b - qn * r.se, r.b - qt * r.se * rse); hi = max(r.b + qn * r.se, r.b + qt * r.se * rse)
        p = max(2 * stats.norm.sf(abs(r.b / r.se)), 2 * stats.t.sf(abs(r.b / r.se / rse), r.df))
        ok = rse < 1 and abs(lo - r.lo_stored) <= TOL and abs(hi - r.hi_stored) <= TOL and abs(p - r.p_stored) <= TOL * r.p_stored
        O.loc[i, ['rule_reproducing', 'p_recomputed', 'lo_recomputed', 'hi_recomputed', 'status']] = [f'package branch for residual SE < 1 (RSE {rse:.4f})', p, lo, hi, 'PASS' if ok else 'FAIL']
O = O.drop(columns=['b', 'se', 'df'])
O.to_csv('step9_mr_output_consistency_S46.csv', index=False)
print(O.status.value_counts().to_dict()); print(O[O.status != 'PASS'].to_string() if (O.status != 'PASS').any() else 'all rows reproduced')
print(O[O.rule_reproducing.astype(str).str.startswith('package branch')][['analysis_key', 'rule_reproducing', 'status']].to_string())
if (O.status != 'PASS').any(): sys.exit('output consistency FAILED')
