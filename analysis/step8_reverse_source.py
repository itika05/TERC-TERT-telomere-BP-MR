"""v6: reverse MR with LTL outcome beta and SE taken from the original Codd 2021 summary statistics (figshare) for every reverse
instrument, matched on GRCh37 position and both alleles and aligned to the Knowledge Portal alternate allele, replacing the portal
LTL records (whose SE was derived from P, with six N/allele-frequency fallbacks). Also builds reverse inputs for the PLINK 1.9 clumped sets
and the forward PLINK-set FinnGen input including the eight variants retrieved from FinnGen R12."""
import numpy as np, pandas as pd
from scipy.stats import norm
from mr_lib import steiger_r2_cont, steiger_direction
N_CODD = 464716
C = pd.concat([pd.read_csv('source_data/codd2021_ukb_ltl_reverse_candidates.tsv', sep='\t'), pd.read_csv('source_data/codd2021_ukb_ltl_subset.tsv', sep='\t')]).drop_duplicates()
C.columns = ['rsid_c', 'p_c', 'chr', 'pos', 'ea', 'oa', 'eaf', 'b_c', 'se_c']; C['chr'] = C.chr.astype(str)
def src(varIds):
    t = pd.Series(varIds).str.split(':', expand=True); d = pd.DataFrame(dict(varId=varIds, chr=t[0], pos=t[1].astype(int), ref=t[2], alt=t[3]))
    m = d.merge(C, on=['chr', 'pos'], how='left')
    al = np.where((m.ea == m.alt) & (m.oa == m.ref), 1, np.where((m.ea == m.ref) & (m.oa == m.alt), -1, 0))
    m['by_src'] = np.where(al != 0, m.b_c * al, np.nan); m['sy_src'] = np.where(al != 0, m.se_c, np.nan); m['p_src'] = np.where(al != 0, m.p_c, np.nan)
    m['match'] = np.select([al == 1, al == -1, m.ea.isna()], ['same', 'swapped', 'absent from source'], 'allele mismatch')
    m['_ok'] = al != 0
    return m.sort_values('_ok', ascending=False).drop_duplicates('varId')[['varId', 'by_src', 'sy_src', 'p_src', 'match']]
RV = pd.read_csv('step3v3_reverse_variant_ledgers_all.csv'); led = []
for ex, pre in [('SBP', 'S_mvp'), ('DBP', 'D_mvp'), ('Hypertension', 'H_mvp')]:
    old = pd.read_csv(f'mr_v4_inputs/reverse_{ex}_L_codd.csv')
    cand = RV[(RV.exposure == ex) & RV.usable].drop_duplicates('rsid')
    plink = set(open(f'plink_clump/{ex}_plink_kept.txt').read().split())
    for sel, ids in [('custom', set(cand.loc[cand['selected_r2_0.001'], 'varId'])), ('plink', set(cand.loc[cand.rsid.isin(plink), 'varId']))]:
        c = cand[cand.varId.isin(ids)].copy(); s = src(c.varId.values); c = c.merge(s, on='varId', how='left')
        nx = c[pre + '_n'].median()
        c['bx'] = c[pre + '_b']; c['sx'] = c[pre + '_b'].abs() / norm.isf(c[pre + '_p'] / 2)
        c['by'] = c.by_src; c['sy'] = c.sy_src
        led.append(dict(exposure=ex, selection=sel, instruments=len(c), matched_same=int((c.match == 'same').sum()), matched_swapped=int((c.match == 'swapped').sum()),
                        absent=int((c.match == 'absent from source').sum()), allele_mismatch=int((c.match == 'allele mismatch').sum())))
        c = c[c.by.notna() & c.sx.notna() & np.isfinite(c.sx)]
        c['steiger_z'] = steiger_direction(steiger_r2_cont(c.bx, c.sx, nx), nx, steiger_r2_cont(c.by, c.sy, N_CODD), N_CODD)[0] if ex != 'Hypertension' else np.nan
        c.assign(exposure=ex, outcome='LTL (Codd 2021 source statistics)')[['exposure', 'outcome', 'varId', 'rsid', 'bx', 'sx', 'by', 'sy', 'p_src', 'match', 'steiger_z']].to_csv(
            f'mr_v6_inputs/reverse_{sel}_{ex}_L_src.csv', index=False)
        if sel == 'custom':
            cmp = old.merge(c[['varId', 'by', 'sy']], on='varId', suffixes=('_portal', '_src'))
            led[-1].update(dict(b_corr_portal_vs_source=np.corrcoef(cmp.by_portal, cmp.by_src)[0, 1], se_ratio_median=float((cmp.sy_portal / cmp.sy_src).median()),
                                se_ratio_min=float((cmp.sy_portal / cmp.sy_src).min()), se_ratio_max=float((cmp.sy_portal / cmp.sy_src).max())))
L = pd.DataFrame(led); L.to_csv('step8_reverse_source_ledger.csv', index=False); print(L.to_string())
# forward PLINK set inputs incl. FinnGen with the 8 retrieved records
k = set(open('plink_clump/forward_plink_kept.txt').read().split())
for code in ['S_mvpEU', 'D_mvpEU', 'H_mvpEU', 'H_fg12']:
    sup = pd.read_csv(f'mr_v4_inputs/forward_selected_r2_0.1_{code}.csv'); d = sup[sup.rsid.isin(k)].copy()
    if code == 'H_fg12':
        F8 = pd.read_csv('source_data/finngen_R12_I9_HYPTENS_plink_only_8.tsv', sep='\t'); V = pd.read_csv('step6_v4_forward_variant_ledger.csv')
        add = []
        for r in F8.itertuples():
            v = V[V.rsid == r.rsids].iloc[0]; ref, alt = v.varId.split(':')[2:4]
            s = 1 if (r.ref, r.alt) == (ref, alt) else (-1 if (r.ref, r.alt) == (alt, ref) else 0)
            add.append(dict(code='H_fg12', trait='Hypertension', source='FinnGen R12', family='secondary', selection='plink', varId=v.varId, rsid=v.rsid, gene=v.gene,
                            bx=v.bx_src, sx=v.sx_src, by=r.beta * s, sy=r.sebeta, p_exposure=v.p_src, p_outcome=r.pval, allele_check=s))
        A = pd.DataFrame(add); print('FinnGen 8 allele alignment', A.allele_check.tolist()); d = pd.concat([d, A[A.allele_check != 0]], ignore_index=True)
    d.to_csv(f'mr_v6_inputs/forward_plink_{code}.csv', index=False); print(code, len(d))
