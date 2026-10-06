"""v5: instrument-flow ledger (forward, reverse, published Codd instrument) reconstructed from the saved candidate files and ledgers."""
import numpy as np, pandas as pd
rows = []
add = lambda d, st, n, note='': rows.append(dict(direction=d, step=st, n=int(n), note=note))
K = pd.read_csv('kp_ltl_leads_wide.csv')
add('forward', 'KP European LTL lead variants (candidate pool)', len(K), 'Knowledge Portal LTL leads; ' + ', '.join(f'{k}: {v}' for k, v in K.src.value_counts().items()))
a = K[K.L_codd_p < 5e-8]; add('forward', 'Codd 2021 (UK Biobank) P < 5e-8', len(a))
b = a[a.rsid.notna()]; add('forward', 'with rsID', len(b))
c = b[np.minimum(b.maf, 1 - b.maf) >= 0.01]; add('forward', 'MAF >= 1%', len(c))
V = pd.read_csv('step6_v4_forward_variant_ledger.csv')
add('forward', 'candidates entering clumping', len(V))
add('forward', 'excluded: chrX (not assessable in 1000G haploid males)', (V.chr.astype(str) == 'X').sum())
add('forward', 'excluded: MHC chr6:25-34 Mb', V.mhc.sum())
add('forward', 'usable (matched to 1000G EUR, autosomal, non-MHC)', (V.usable & ~V.mhc).sum())
add('forward', 'selected after clumping r2 < 0.001 within 10 Mb (primary)', V['selected_r2_0.001'].sum(), f'{(V.usable & ~V.mhc).sum() - V["selected_r2_0.001"].sum()} removed for r2 >= 0.001 with a stronger instrument')
add('forward', 'selected after clumping r2 < 0.1 (correlated-instrument sensitivity)', V['selected_r2_0.1'].sum())
import glob
for f in sorted(glob.glob('mr_v4_inputs/forward_selected_r2_0.001_*.csv')):
    g = pd.read_csv(f); add('forward', f'primary instruments with outcome data: {g.code.iloc[0]} ({g.source.iloc[0]})', len(g))
R3 = pd.read_csv('step3v3_reverse_ledger.csv'); Rv = pd.read_csv('step6_v4_reverse_inputs_all.csv')
for r in R3.itertuples():
    add('reverse', f'{r.exposure}: MVP European leads, P < 5e-8, rsID, MAF >= 1% (candidates)', r.candidates)
    add('reverse', f'{r.exposure}: usable in 1000G EUR', r.usable)
    add('reverse', f'{r.exposure}: selected after clumping r2 < 0.001 within 10 Mb', r.selected)
    for o, g in Rv[Rv.exposure == r.exposure].groupby('outcome'): add('reverse', f'{r.exposure}: with outcome data, {o}', len(g))
P = pd.read_csv('step6_codd_published_instrument_lookup.csv')
add('published', 'Codd 2021 Supplementary Table 1 sentinels flagged MR = Y', len(P))
add('published', 'excluded: MHC', P.mhc.sum())
add('published', 'not found in Knowledge Portal (no portal varId)', P.portal_varId.isna().sum())
for t in ['SBP', 'DBP', 'HYPERTENSION']: add('published', f'with MVP European {t} data', P[f'{t}_b'].notna().sum())
for t in ['SBP', 'DBP', 'HYPERTENSION']:
    import glob as _g; f = f'mr_v4_inputs/forward_codd_published_{t[0]}_mvpEU.csv'; add('published', f'analysed for {t}', len(pd.read_csv(f)), 'beta = 0 records cannot yield an SE from P and are dropped' if t == 'HYPERTENSION' else '')
F = pd.DataFrame(rows); F.to_csv('step7_instrument_flow_ledger.csv', index=False); print(F.to_string())
# per-variant reasons for the published instrument
P['reason_not_analysed'] = np.select([P.portal_varId.isna(), P.SBP_b.isna()], ['no Knowledge Portal record at this position/allele pair', 'no MVP European SBP/DBP/hypertension record returned by the Knowledge Portal for this variant'], '')
P[['rsid', 'chr', 'pos', 'gene', 'effect_allele', 'other_allele', 'mhc', 'portal_varId', 'reason_not_analysed', 'HYPERTENSION_b']].to_csv('step7_published_instrument_reasons.csv', index=False)
print(P.reason_not_analysed.value_counts())
