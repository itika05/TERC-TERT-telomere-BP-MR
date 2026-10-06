"""v7 SE-fallback ledger (Supplementary Table S38): every record whose SE could not be recovered from P, with its status in the
current results. Reverse-MR portal LTL records are marked superseded (current reverse MR uses the source SEs; S13). MVMR records are
listed per trait (P <= 1e-300, P >= 0.999 or beta = 0) and are excluded from the only MVMR analysis now reported (318 variants)."""
import pandas as pd, numpy as np
L = pd.read_csv('step7_se_fallback_ledger.csv')
L['analysis_scope'] = np.where(L.direction == 'reverse', 'reverse MR (portal Codd LTL records, v4-v5)', np.where(L.file.str.contains('codd_published'), 'forward MR, published Codd instrument', 'forward MR, MVP African American transportability'))
S13 = pd.read_csv('step9_reverse_inputs_S13.csv')
def status(r):
    if r.direction == 'reverse':
        used = S13[(S13.code == f'{r.exposure}->L_codd') & (S13.rsid == r.rsid)]
        src = f'; current source SE {used.sy.iloc[0]:.6g}' if len(used) else ''
        return 'superseded: not used in any current result (reverse MR uses Codd source beta/SE, S13' + src + ')'
    if not r.fallback_used: return 'no fallback applied (exposure SE from source statistics); listed for completeness'
    return 'current: used in the secondary MVP African American SBP transportability analysis (quantitative trait)'
L['status_current_release'] = L.apply(status, axis=1)
M = pd.read_csv('mr_v4_inputs/mvmr_instruments.csv'); rows = []
for tr in ['LTL', 'BMI', 'LYM', 'S_mvpEU', 'D_mvpEU', 'H_mvpEU']:
    p, b = M[tr + '_p'], M[tr + '_b']; fb = (p <= 1e-300) | (p >= 0.999) | (b == 0)
    for r in M[fb].itertuples():
        rows.append(dict(direction='MVMR', varId=r.varId, trait=tr, reason='P <= 1e-300' if getattr(r, tr + '_p') <= 1e-300 else ('P >= 0.999' if getattr(r, tr + '_p') >= 0.999 else 'beta = 0'),
                         formula_used='SE = 1/sqrt(2 N f (1-f)) (valid only for standardised quantitative traits)', fallback_used=True,
                         analysis_scope='MVMR (330 jointly clumped variants)',
                         status_current_release='excluded: variant removed from the 318-variant MVMR, the only MVMR analysis reported'))
MV = pd.DataFrame(rows)
assert MV.varId.nunique() == 12, MV.varId.nunique()
assert set(MV.varId) == set(pd.read_csv('step7_mvmr_fallback_variants.csv').varId)
out = pd.concat([L, MV], ignore_index=True)
front = ['direction', 'analysis_scope', 'status_current_release', 'code', 'exposure', 'outcome', 'trait', 'varId', 'rsid', 'reason', 'formula_used', 'fallback_used']
out = out[front + [c for c in out.columns if c not in front]]
out.to_csv('step9_se_fallback_ledger_S38.csv', index=False)
print(out.groupby(['direction', 'status_current_release']).size().to_string())
print('MVMR fallback records by trait:', MV.groupby('trait').size().to_dict())
