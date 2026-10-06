"""Colocalization inputs, version 4 (audit major concern 5).
LTL: original UK Biobank summary statistics (Codd et al. 2021; beta and SE as published, so no P-value floor problem).
Blood pressure, the same samples as the primary MR: MVP European SBP, DBP and hypertension (Verma et al. 2024), retrieved per variant from the
Knowledge Portal (beta, P, N; SE = |beta|/z(P)). Independent hypertension sample: FinnGen R12 I9_HYPTENS (beta, SE; GRCh38, joined to GRCh37 by rsID
and alleles). None of these outcome samples includes UK Biobank.
Window centre: the variant with the largest |beta/SE| for LTL in the source data. Windows +/-50, 100 and 250 kb (MVP regional data were retrieved for
the +/-250 kb windows). Variants are kept per trait pair when present in both datasets with non-zero beta and finite SE; for SuSiE they must also be
matched to the 1000 Genomes EUR panel (signed LD aligned to the Knowledge Portal alternate allele)."""
import os, json, numpy as np, pandas as pd
from scipy.stats import norm
from ld_ref import load, match, hap_matrix, corr

os.makedirs('coloc_v4_inputs', exist_ok=True)
C = pd.read_csv('source_data/codd2021_ukb_ltl_subset.tsv', sep='\t'); C.columns = ['rsid', 'p', 'chr', 'pos', 'ea', 'oa', 'eaf', 'beta', 'se']; C['chr'] = C.chr.astype(str)
reg = json.load(open('source_data/kp_v4_instrument_phenome_and_mvp_regional.json'))['reg']
F = pd.read_csv('source_data/finngen_R12_I9_HYPTENS_subset_hg38.tsv', sep='\t'); F.columns = ['chr38', 'pos38', 'ref38', 'alt38', 'rsid', 'p', 'beta', 'se', 'af']
LOCI = {'TERC': ('3', 168_900_000, 170_100_000), 'TERT': ('5', 700_000, 1_900_000)}
N = {'SBP': 425740, 'DBP': 425743, 'HYPERTENSION': 318398}

rows = []
for v, recs in reg.items():
    c, p, a0, a1 = v.split(':')
    for ds, ph, b, pv, n in recs:
        if ds == 'Verma2024_MVPTraits_EU': rows.append(dict(varId=v, chr=c, pos=int(p), ref=a0, alt=a1, trait=ph, beta=b, p=pv, n=n))
M = pd.DataFrame(rows)
M['se'] = np.abs(M.beta) / norm.isf(M.p.clip(lower=1e-300, upper=0.999) / 2)
M = M[(M.beta != 0) & (M.p < 0.999)]

ledger, summ = [], []
for loc, (ch, a, b) in LOCI.items():
    ref = load(f'ldref/eur_region_{loc}.json')
    L = C[(C.chr == ch) & C.pos.between(a, b)].copy()
    # express LTL on the KP varId convention (ref:alt with alt as effect) using the MVP/1000G variant list, else ea as alt
    L['z'] = L.beta / L.se
    centre = L.loc[L.z.abs().idxmax()]
    for w in [50, 100, 250]:
        lo, hi = centre.pos - w * 1000, centre.pos + w * 1000
        Lw = L[L.pos.between(lo, hi)]
        for tr in ['SBP', 'DBP', 'HYPERTENSION', 'FG_HYPERTENSION']:
            if tr == 'FG_HYPERTENSION':
                O = F.merge(Lw[['rsid', 'chr', 'pos']], on='rsid')
                O = O.assign(varId=O.chr + ':' + O.pos.astype(str) + ':' + O.ref38 + ':' + O.alt38, ref=O.ref38, alt=O.alt38, n=500264)
            else:
                O = M[(M.trait == tr) & (M.chr == ch) & M.pos.between(lo, hi)].copy()
            J = O.merge(Lw[['chr', 'pos', 'ea', 'oa', 'beta', 'se', 'p', 'rsid']].rename(columns={'beta': 'b_ltl', 'se': 'se_ltl', 'p': 'p_ltl', 'rsid': 'rsid_ltl'}), on=['chr', 'pos'])
            s = np.where((J.ea == J.alt) & (J.oa == J.ref), 1, np.where((J.ea == J.ref) & (J.oa == J.alt), -1, 0))
            J = J[s != 0].assign(b_ltl=lambda x: x.b_ltl * s[s != 0]).drop_duplicates('varId')
            J = J[np.isfinite(J.se) & (J.se > 0) & (J.se_ltl > 0)].sort_values('pos')
            m = match(J.varId, ref); J = J.merge(m[['varId', 'panel_status', 'eur_af_alt']], on='varId')
            J['in_panel'] = J.panel_status.str.startswith('matched')
            J[['varId', 'rsid_ltl', 'pos', 'b_ltl', 'se_ltl', 'p_ltl', 'beta', 'se', 'p', 'n', 'in_panel', 'eur_af_alt']].rename(
                columns={'rsid_ltl': 'rsid', 'beta': 'b_out', 'se': 'se_out', 'p': 'p_out', 'n': 'n_out'}).to_csv(f'coloc_v4_inputs/{loc}_{w}kb_{tr}.csv', index=False)
            P = J[J.in_panel]
            if w in (100, 250) and len(P):
                mm = m[m.varId.isin(P.varId)].set_index('varId').loc[P.varId].reset_index()
                pd.DataFrame(corr(hap_matrix(mm, ref)), index=P.varId, columns=P.varId).to_csv(f'coloc_v4_inputs/{loc}_{w}kb_{tr}_LD.csv')
            summ.append(dict(locus=loc, window_kb=w, centre=centre.rsid, centre_pos=int(centre.pos), centre_abs_z=abs(centre.z), trait=tr,
                             n_variants=len(J), n_in_1000G_panel=int(J.in_panel.sum()), min_p_ltl=J.p_ltl.min() if len(J) else np.nan, min_p_trait=J.p.min() if len(J) else np.nan))
S = pd.DataFrame(summ); S.to_csv('step6_coloc_v4_window_ledger.csv', index=False); print(S.to_string())
