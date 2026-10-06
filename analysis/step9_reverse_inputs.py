"""v7: one reverse-MR input table (Supplementary Table S13) holding exactly the instrument rows used by every reverse-MR result in
step8_mr_results_v6.csv, with raw source records, aligned values, allele-match status and set membership, followed by a reconstruction
test: IVW (multiplicative random effects, residual scale floored at 1, t(k-1)) is recomputed from S13 for every reverse IVW result and
compared with the stored result. The build fails if any difference exceeds the declared tolerance."""
import sys
import numpy as np, pandas as pd
from scipy.stats import t as tdist

C = pd.concat([pd.read_csv('source_data/codd2021_ukb_ltl_reverse_candidates.tsv', sep='\t'),
               pd.read_csv('source_data/codd2021_ukb_ltl_subset.tsv', sep='\t')]).drop_duplicates()
C.columns = ['src_rsid', 'src_p', 'chr', 'pos', 'src_effect_allele', 'src_other_allele', 'src_eaf', 'src_beta', 'src_se']
C['chr'] = C.chr.astype(str)
rows = []
for ex in ['SBP', 'DBP', 'Hypertension']:
    cu = pd.read_csv(f'mr_v6_inputs/reverse_custom_{ex}_L_src.csv'); pl = pd.read_csv(f'mr_v6_inputs/reverse_plink_{ex}_L_src.csv')
    u = pd.concat([cu, pl]).drop_duplicates('varId').drop(columns=['steiger_z'])
    # Steiger z is recomputed per set in step8_reverse_source.py from the same values; take it from the custom file where present
    u = u.merge(pd.concat([cu, pl]).drop_duplicates('varId')[['varId', 'steiger_z']], on='varId', how='left')
    u['in_custom_primary'] = u.varId.isin(cu.varId)
    u['in_plink'] = u.varId.isin(pl.varId)
    u['in_steiger_filtered'] = u.in_custom_primary & (u.steiger_z > 0) if ex != 'Hypertension' else False
    t = u.varId.str.split(':', expand=True); u['chr'] = t[0]; u['pos'] = t[1].astype(int); u['ref'] = t[2]; u['alt_effect'] = t[3]
    u = u.merge(C.drop_duplicates(['chr', 'pos', 'src_effect_allele', 'src_other_allele']), on=['chr', 'pos'], how='left')
    ok = ((u.src_effect_allele == u.alt_effect) & (u.src_other_allele == u.ref)) | ((u.src_effect_allele == u.ref) & (u.src_other_allele == u.alt_effect))
    u = u[ok | u.src_beta.isna()].drop_duplicates('varId')
    u['outcome_source'] = 'Codd 2021 UK Biobank LTL source summary statistics (figshare)'
    u['exposure_source'] = 'Knowledge Portal Verma2024_MVPTraits_EU; SE = |beta| / norm.isf(P/2)'
    u['code'] = f'{ex}->L_codd'
    rows.append(u)
for ex in ['SBP', 'DBP', 'Hypertension']:
    o = pd.read_csv(f'mr_v4_inputs/reverse_{ex}_L_nkSA.csv')
    o['code'] = f'{ex}->L_nkSA'; o['in_custom_primary'] = True; o['in_plink'] = False
    o['in_steiger_filtered'] = (o.steiger_z > 0) if ex != 'Hypertension' else False
    o['outcome_source'] = 'Knowledge Portal Nakao2026_LTL_SA (exploratory; SE from P)'
    o['exposure_source'] = 'Knowledge Portal Verma2024_MVPTraits_EU; SE = |beta| / norm.isf(P/2)'; o['match'] = 'portal (aligned to portal alternate allele)'
    rows.append(o)
S = pd.concat(rows, ignore_index=True)
cols = ['code', 'exposure', 'outcome', 'varId', 'rsid', 'in_custom_primary', 'in_plink', 'in_steiger_filtered', 'bx', 'sx', 'by', 'sy', 'steiger_z', 'match',
        'src_rsid', 'src_effect_allele', 'src_other_allele', 'src_eaf', 'src_beta', 'src_se', 'src_p', 'outcome_source', 'exposure_source']
S = S[[c for c in cols if c in S.columns]]
S.to_csv('step9_reverse_inputs_S13.csv', index=False)

# reconstruction test against stored results
def ivw(d):
    w = 1 / d.sy ** 2; b = (w * d.bx * d.by).sum() / (w * d.bx ** 2).sum(); k = len(d)
    rss = (w * (d.by - b * d.bx) ** 2).sum(); phi = max(1.0, rss / (k - 1)); se = np.sqrt(phi / (w * d.bx ** 2).sum())
    q = tdist.ppf(0.975, k - 1); p = 2 * tdist.sf(abs(b / se), k - 1)
    return b, se, b - q * se, b + q * se, p, k
R = pd.read_csv('step8_mr_results_v6.csv'); TOL = 1e-8; out = []
for (code, analysis), flag in [((c, 'r2<0.001 (primary)'), 'in_custom_primary') for c in S.code.unique()] + \
        [((c, 'PLINK 1.9 clumping'), 'in_plink') for c in S.code.unique() if c.endswith('L_codd')] + \
        [((c, 'r2<0.001, Steiger-filtered'), 'in_steiger_filtered') for c in S.code.unique() if not c.startswith('Hypertension')]:
    st = R[(R.code == code) & (R.analysis == analysis) & (R.method == 'IVW (MRE, t)')]
    if st.empty: continue
    st = st.iloc[0]; d = S[(S.code == code) & (S[flag] == True)]
    b, se, lo, hi, p, k = ivw(d)
    rel = max(abs(b - st.b) / max(abs(st.b), 1e-12), abs(se - st.se) / st.se, abs(p - st.p) / max(st.p, 1e-300))
    out.append(dict(code=code, analysis=analysis, k_input=k, k_stored=int(st.nsnp), b_input=b, b_stored=st.b, se_input=se, se_stored=st.se,
                    p_input=p, p_stored=st.p, max_rel_diff=rel, tolerance=TOL, status='PASS' if (k == st.nsnp and rel <= TOL) else 'FAIL'))
T = pd.DataFrame(out); T.to_csv('step9_reverse_reconstruction_test.csv', index=False)
print(T[['code', 'analysis', 'k_input', 'k_stored', 'b_input', 'b_stored', 'max_rel_diff', 'status']].to_string())
if (T.status != 'PASS').any(): sys.exit('reverse-MR reconstruction FAILED')
