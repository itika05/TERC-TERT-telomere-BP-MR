"""v5: validate standard errors reconstructed from P (SE = |beta| / Phi^-1(1 - P/2)) against reported SEs.
(1) Knowledge Portal top-1,500 records per dataset (dataset-associations endpoint, which reports stdErr) for MVP EU and 4 other BP GWAS;
(2) Codd 2021 LTL: portal P-derived SE vs SE in the original figshare summary statistics for the 188 forward candidates;
(3) ledger of every instrument record in the primary analyses whose SE could not be taken from P (P floor, P>=0.999 or beta = 0)."""
import json, numpy as np, pandas as pd
from scipy.stats import norm
d = json.load(open('source_data/kp_v4_round2.json'))['cal']
rows = []
for k, recs in d.items():
    ds, tr = k.split('|'); X = pd.DataFrame(recs).iloc[:, :5]; X.columns = ['varId', 'beta', 'se_rep', 'p', 'n']
    X = X[(X.p > 1e-300) & (X.beta != 0) & (X.se_rep > 0)]
    X['se_P'] = X.beta.abs() / norm.isf(X.p / 2); r = X.se_P / X.se_rep
    rows.append(dict(source=f'Knowledge Portal {ds}', trait=tr, n_records=len(X), p_range=f'{X.p.min():.1e} to {X.p.max():.1e}', ratio_median=r.median(),
                     ratio_p2_5=r.quantile(.025), ratio_p97_5=r.quantile(.975), max_abs_pct_diff=100 * (r - 1).abs().max(),
                     pct_within_1pct=100 * ((r - 1).abs() < 0.01).mean(), pct_within_5pct=100 * ((r - 1).abs() < 0.05).mean()))
C = pd.read_csv('step6_codd_source_vs_portal.csv').dropna(subset=['se_source'])
r = C.se_portal_from_P / C.se_source
rows.append(dict(source='Codd 2021 LTL: portal P-derived vs original summary statistics', trait='LTL', n_records=len(C), p_range=f'{C.p_portal.min():.1e} to {C.p_portal.max():.1e}',
                 ratio_median=r.median(), ratio_p2_5=r.quantile(.025), ratio_p97_5=r.quantile(.975), max_abs_pct_diff=100 * (r - 1).abs().max(),
                 pct_within_1pct=100 * ((r - 1).abs() < 0.01).mean(), pct_within_5pct=100 * ((r - 1).abs() < 0.05).mean()))
R = pd.DataFrame(rows); R.to_csv('step7_se_validation.csv', index=False); print(R.to_string())
# fallback ledger. Forward: outcome SE from P; reverse: flags written by step6 (outcome_p_floor / outcome_p_one)
F = pd.read_csv('step6_v4_harmonised_inputs_all.csv')
Fl = F[(F.p_outcome >= 0.999) | (F.by == 0) | (F.p_outcome <= 1e-300) | (F.p_exposure <= 1e-300)].copy()
Fl['reason'] = np.select([Fl.p_outcome >= 0.999, Fl.by == 0, Fl.p_outcome <= 1e-300], ['outcome P >= 0.999', 'outcome beta = 0', 'outcome P at 1e-300 floor'], 'exposure P at 1e-300 floor')
Fl['direction'] = 'forward'
Rv = pd.read_csv('step6_v4_reverse_inputs_all.csv')
Rl = Rv[Rv.outcome_p_floor.astype(bool) | Rv.outcome_p_one.astype(bool)].copy()
Rl['reason'] = np.where(Rl.outcome_p_floor.astype(bool), 'outcome (LTL) P at 1e-300 floor', 'outcome (LTL) P >= 0.999 or beta = 0'); Rl['direction'] = 'reverse'
L = pd.concat([Fl, Rl], ignore_index=True)
src = L.reason.eq('exposure P at 1e-300 floor')
L['formula_used'] = np.where(src, 'none: exposure SE taken from the source summary statistics (only P is floored in the published table); listed for completeness', 'SE = 1/sqrt(2 N f (1-f)) (standardised trait; N = GWAS N; f = 1000 Genomes EUR allele frequency)')
L['fallback_used'] = ~src
L.to_csv('step7_se_fallback_ledger.csv', index=False)
print('fallback records:', int(L.fallback_used.sum()), 'listed:', len(L)); print(L.groupby(['direction', 'reason']).size().to_string() if len(L) else '')
print(L.groupby(['direction','selection' if 'selection' in L else 'direction']).size() if len(L) else '')
