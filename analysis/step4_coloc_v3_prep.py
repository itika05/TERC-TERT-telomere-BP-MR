"""Step 4 v3 inputs: regional data with SEs that avoid the P floor (P-derived where valid, portal-reported where P is floored), a declared window-centre rule that does not depend on P-value
ties (centre = variant with the largest |beta/SE| for LTL in the European regional data), a P-floor audit, and signed LD matrices
from 1000 Genomes phase 3 EUR haplotypes aligned to the KP alternate (effect) allele, for coloc (R package) and SuSiE."""
import os, numpy as np, pandas as pd
from ld_ref import load, match, hap_matrix, corr

os.makedirs('coloc_v3_inputs', exist_ok=True)
from scipy.stats import norm
r = pd.read_csv('kp_regional_TERC_TERT.csv'); r['pos'] = r.varId.str.split(':').str[1].astype(int)
# SE used: P-derived where P is above the double-precision floor (full precision; the portal's reported SE is rounded to ~2 significant
# figures), the portal-reported SE where P sits at the floor (P-derived SE would be inflated there). Sensitivity: reported SE throughout.
r['se_reported'] = r.se; r['floor'] = r.p <= 1e-300
r['se_p'] = np.abs(r.beta) / norm.isf(r.p.clip(lower=1e-300) / 2)
r['se'] = np.where(r.floor, r.se_reported, r.se_p); r['z'] = r.beta / r.se
audit = []
for (loc, tr, an), g in r.groupby(['locus', 'trait', 'anc']):
    minp = g.p.min(); zs = g.z.abs().sort_values(ascending=False)
    audit.append(dict(locus=loc, trait=tr, ancestry=an, n_variants=len(g), n_at_P_floor=int((g.p <= 1e-300).sum()), n_tied_at_min_P=int((g.p == minp).sum()),
                      min_P=minp, top_by_z=g.loc[zs.index[0], 'rsid'], top_pos=int(g.loc[zs.index[0], 'pos']), top_abs_z=zs.iloc[0], second_abs_z=zs.iloc[1],
                      old_centre_rule_min_P_first_listed=g.loc[g.p.idxmin(), 'rsid']))
A = pd.DataFrame(audit); A.to_csv('step4v3_pfloor_audit.csv', index=False); print(A.to_string())

rows = []
for loc in ['TERC', 'TERT']:
    ref = load(f'ldref/eur_region_{loc}.json')
    L = r[(r.locus == loc) & (r.trait == 'LTL') & (r.anc == 'EU')]
    centre = L.loc[L.z.abs().idxmax()]
    for w in [50, 100, 250, 500]:
        W = r[(r.locus == loc) & (r.anc == 'EU') & ((r.pos - centre.pos).abs() <= w * 1000)]
        ok = W[(W.beta != 0) & (W.se > 0) & np.isfinite(W.se) & (W.se_reported > 0)]
        common = set.intersection(*[set(ok[ok.trait == t].varId) for t in ['LTL', 'SBP', 'DBP', 'HYPERTENSION']])
        m = match(sorted(common), ref)
        m = m[m.panel_status.str.startswith('matched')]
        m = m.assign(pos=m.varId.str.split(':').str[1].astype(int)).sort_values('pos')
        for t in ['LTL', 'SBP', 'DBP', 'HYPERTENSION']:
            x = W[W.trait == t].set_index('varId').loc[m.varId]
            x[['rsid', 'beta', 'se', 'se_reported', 'se_p', 'floor', 'p', 'n', 'afEU']].assign(pos=m.pos.values, eur_af_alt=m.eur_af_alt.values).to_csv(f'coloc_v3_inputs/{loc}_{w}kb_{t}.csv')
        if w in (250, 100):
            H = hap_matrix(m, ref); R = corr(H)
            pd.DataFrame(R, index=m.varId, columns=m.varId).to_csv(f'coloc_v3_inputs/{loc}_{w}kb_LD.csv')
        rows.append(dict(locus=loc, window_kb=w, centre=centre.rsid, centre_pos=int(centre.pos), variants_all_traits=len(common), in_panel=len(m)))
pd.DataFrame(rows).to_csv('step4v3_window_ledger.csv', index=False); print(pd.DataFrame(rows))
