"""Multivariable MR inputs (audit major concerns 1 and 4): LTL with body-mass index and lymphocyte count.
Exposures: LTL, Codd et al. 2021 (Knowledge Portal record, identical to source); BMI, Pulit et al. 2019 GIANT+UK Biobank (GWAS_UKBiobankGIANT_eu);
lymphocyte count, Chen et al. 2020 European (Chen2020_BCX_eu). Outcomes: MVP European SBP, DBP, hypertension (no UK Biobank).
Candidate instruments: the 116 LTL instruments, plus the 300 strongest European Knowledge Portal lead variants for BMI and for lymphocyte count
after 1-Mb distance pruning (autosomal, MAF >= 1%, MHC excluded); a candidate is kept only if P < 5e-8 for its trait in the dataset used.
The pooled set is clumped at r2 < 0.001 within 10 Mb using measured 1000 Genomes EUR LD, ordered by the smallest P across the three exposures.
SE: |beta|/z(P); where P <= 1e-300, P >= 0.999 or beta = 0, 1/sqrt(2 N f (1 - f)) with f from 1000 Genomes EUR."""
import json, numpy as np, pandas as pd
from scipy.stats import norm
from ld_ref import load, match, hap_matrix, corr

ref = pd.concat([load('ldref/eur_instruments.json'), load('ldref/eur_mvmr_candidates.json')]).drop_duplicates(['chr', 'pos', 'ref', 'alt']).reset_index(drop=True)
r1 = json.load(open('source_data/kp_v4_instrument_phenome_and_mvp_regional.json'))['inst']
r2 = json.load(open('source_data/kp_v4_round2.json'))
recs = {**{k: v for k, v in r2['r2'].items()}, **r1}
meta = r2['mvmeta']
EXP = {'LTL': ('Codd2021_LTL_Mixed', 'LTL'), 'BMI': ('GWAS_UKBiobankGIANT_eu', 'BMI'), 'LYM': ('Chen2020_BCX_eu', 'LymphoCount')}
OUT = {'S_mvpEU': ('Verma2024_MVPTraits_EU', 'SBP'), 'D_mvpEU': ('Verma2024_MVPTraits_EU', 'DBP'), 'H_mvpEU': ('Verma2024_MVPTraits_EU', 'HYPERTENSION')}
ltl = pd.read_csv('mr_v4_inputs/forward_selected_r2_0.001_S_mvpEU.csv')
cand = pd.DataFrame({'varId': list(ltl.varId) + [m[0] for m in meta['BMI']] + [m[0] for m in meta['LymphoCount']],
                     'source_trait': ['LTL'] * len(ltl) + ['BMI'] * len(meta['BMI']) + ['LYM'] * len(meta['LymphoCount'])}).drop_duplicates('varId')
rows = []
for r in cand.itertuples():
    rr = recs.get(r.varId, []); row = dict(varId=r.varId, source_trait=r.source_trait)
    for k, (ds, ph) in {**EXP, **OUT}.items():
        hit = [x for x in rr if x[0] == ds and x[1] == ph]
        if hit: row[k + '_b'], row[k + '_p'], row[k + '_n'] = hit[0][2], hit[0][3], hit[0][4]
    rows.append(row)
D = pd.DataFrame(rows)
m = match(D.varId, ref); D = D.merge(m, on='varId')
f = D.eur_af_alt.clip(0.005, 0.995)
for k in list(EXP) + list(OUT):
    p = D[k + '_p']; b = D[k + '_b']
    nf = 1 / np.sqrt(2 * D[k + '_n'] * f * (1 - f))  # used where P is floored (<= 1e-300), P ~ 1 or beta = 0 (SE not recoverable from P)
    D[k + '_se'] = np.where((p > 1e-300) & (p < 0.999) & (b != 0), np.abs(b) / norm.isf(np.clip(p, 1e-300, 0.999) / 2), nf)
D['own_p'] = [row[{'LTL': 'LTL', 'BMI': 'BMI', 'LYM': 'LYM'}[row.source_trait] + '_p'] for _, row in D.iterrows()]
D['usable'] = D.panel_status.str.startswith('matched') & (D.own_p < 5e-8) & D[[k + '_b' for k in list(EXP) + list(OUT)]].notna().all(axis=1)
D['min_p'] = D[[k + '_p' for k in EXP]].min(axis=1)
U = D[D.usable].sort_values('min_p').reset_index(drop=True)
t = U.varId.str.split(':', expand=True); U['chr'] = t[0]; U['pos'] = t[1].astype(int)
R = corr(hap_matrix(U, ref))
kept = []
for i in range(len(U)):
    if any(U.chr[j] == U.chr[i] and abs(U.pos[j] - U.pos[i]) < 10_000_000 and R[i, j] ** 2 >= 0.001 for j in kept): continue
    kept.append(i)
K = U.iloc[kept].copy()
D.to_csv('step6_mvmr_candidate_ledger.csv', index=False); K.to_csv('mr_v4_inputs/mvmr_instruments.csv', index=False)
print('candidates', len(D), 'usable', len(U), 'kept', len(K), K.source_trait.value_counts().to_dict())
