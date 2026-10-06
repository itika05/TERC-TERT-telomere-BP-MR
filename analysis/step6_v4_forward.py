"""Version 4 forward MR inputs (reviewer/editor audit of v3).

Changes from v3
  * The MHC (chr6:25-34 Mb, GRCh37) is excluded from the forward candidates before clumping, the same rule as for the reverse
    instruments (audit major concern 4).
  * Exposure estimates are taken from the original UK Biobank LTL summary statistics of Codd et al. (2021; figshare
    https://figshare.com/s/caa99dc0f76d62990195, BOLT-LMM beta and standard error), aligned to the Knowledge Portal alternate allele.
    Portal-derived values were identical (r = 1.0000; SE ratio 0.999-1.001), which is recorded in step6_codd_source_vs_portal.csv.
  * FinnGen: the Knowledge Portal FinnGen hypertension record is an early 2018 release (N = 52,313). It is replaced by FinnGen R12
    I9_HYPTENS (154,630 cases, 345,634 controls; public summary statistics, GRCh38) matched by rsID and alleles.
  * A second instrument set is the published LD-pruned MR instrument of Codd et al. (Supplementary Table 1, column 'MR' = Y).
Outputs: mr_v4_inputs/forward_{set}_{code}.csv, step6_v4_forward_variant_ledger.csv, step6_v4_forward_clump_log.csv,
         step6_v4_forward_pair_ledger.csv, step6_codd_source_vs_portal.csv, mr_v4_inputs/forward_r2_0.1_correlation.csv"""
import os, json, numpy as np, pandas as pd
from scipy.stats import norm
import step3_mr_v3 as s3
from mr_lib import steiger_r2_cont, steiger_r2_binary, steiger_direction

os.makedirs('mr_v4_inputs', exist_ok=True)
N_CODD = 464716
se_p = lambda b, p: np.abs(b) / norm.isf(np.clip(p, 1e-300, 1) / 2)
MHC = lambda chr_, pos: (chr_ == '6') & (pos >= 25_000_000) & (pos <= 34_000_000)
TT = lambda x: ((x.chr == '3') & x.pos.between(168_982_000, 169_983_000)) | ((x.chr == '5') & x.pos.between(753_000, 1_795_000))

# ---------------- source summary statistics
C = pd.read_csv('source_data/codd2021_ukb_ltl_subset.tsv', sep='\t')
C.columns = ['rsid_c', 'p_c', 'chr', 'pos', 'ea', 'oa', 'eaf', 'beta_c', 'se_c']
C['chr'] = C.chr.astype(str)
F = pd.read_csv('source_data/finngen_R12_I9_HYPTENS_subset_hg38.tsv', sep='\t')
F.columns = ['chr38', 'pos38', 'ref38', 'alt38', 'rsid', 'p_fg', 'b_fg', 'se_fg', 'af_fg']


def add_source(d):
    """d has varId (chr:pos:ref:alt, GRCh37, KP alt = effect). Adds bx_src, sx_src, p_src aligned to alt."""
    t = d.varId.str.split(':', expand=True)
    d = d.assign(chr=t[0], pos=t[1].astype(int), ref=t[2], alt=t[3])
    m = d.merge(C, on=['chr', 'pos'], how='left')
    al = np.where((m.ea == m.alt) & (m.oa == m.ref), 1, np.where((m.ea == m.ref) & (m.oa == m.alt), -1, 0))
    m['bx_src'] = np.where(al != 0, m.beta_c * al, np.nan); m['sx_src'] = np.where(al != 0, m.se_c, np.nan); m['p_src'] = np.where(al != 0, m.p_c, np.nan)
    m['_ok'] = al != 0
    return m.sort_values('_ok', ascending=False).drop_duplicates('varId').drop(columns=['rsid_c', 'p_c', 'ea', 'oa', 'eaf', 'beta_c', 'se_c', '_ok'])


def add_finngen(d):
    """Match FinnGen R12 by rsID; align FinnGen alt to KP alt."""
    m = d.merge(F, on='rsid', how='left')
    al = np.where((m.alt38 == m.alt) & (m.ref38 == m.ref), 1, np.where((m.alt38 == m.ref) & (m.ref38 == m.alt), -1, 0))
    m['H_fg12_b'] = np.where(al != 0, m.b_fg * al, np.nan); m['H_fg12_se'] = np.where(al != 0, m.se_fg, np.nan); m['H_fg12_p'] = np.where(al != 0, m.p_fg, np.nan)
    m['_ok'] = al != 0
    return m.sort_values('_ok', ascending=False).drop_duplicates('varId').drop(columns=['chr38', 'pos38', 'ref38', 'alt38', 'p_fg', 'b_fg', 'se_fg', 'af_fg', '_ok'])


if __name__ == '__main__':
    led3 = pd.read_csv('step3v3_forward_variant_ledger.csv', dtype={'chr': str})
    base = led3.drop(columns=['panel_status', 'sign', 'palindromic', 'eur_af_alt', 'usable', 'selected_r2_0.001', 'selected_r2_0.1', 'chr', 'pos'])
    fw = s3.prep(base)
    fw['mhc'] = MHC(fw.chr, fw.pos)
    fw.loc[fw.mhc, 'panel_status'] = 'excluded: MHC (chr6:25-34 Mb)'
    fw.loc[fw.mhc, 'usable'] = False
    kept, drops, pairs, Rm, posmap = s3.clump(fw, 'L_codd_p')
    fw['selected_r2_0.001'] = fw.varId.isin(kept)
    k01, _, _, _, _ = s3.clump(fw, 'L_codd_p', r2max=0.1)
    fw['selected_r2_0.1'] = fw.varId.isin(k01)
    fw = add_source(fw.drop(columns=['chr', 'pos'])); fw = add_finngen(fw)
    fw.drop(columns=['idx']).to_csv('step6_v4_forward_variant_ledger.csv', index=False)
    drops.to_csv('step6_v4_forward_clump_log.csv', index=False); pairs.to_csv('step6_v4_forward_pair_ledger.csv', index=False)
    print('candidates', len(fw), 'MHC excluded', int(fw.mhc.sum()), 'kept r2<0.001', len(kept), 'kept r2<0.1', len(k01))
    S01 = fw[fw['selected_r2_0.1']].sort_values('L_codd_p'); idx = [posmap[v] for v in S01.varId]
    pd.DataFrame(Rm[np.ix_(idx, idx)], index=S01.varId, columns=S01.varId).to_csv('mr_v4_inputs/forward_r2_0.1_correlation.csv')

    # source vs portal check
    sel = fw[fw['selected_r2_0.001']]
    sp = pd.DataFrame(dict(varId=sel.varId, rsid=sel.rsid, b_portal=sel.L_codd_b, se_portal_from_P=se_p(sel.L_codd_b, sel.L_codd_p),
                           b_source=sel.bx_src, se_source=sel.sx_src, p_portal=sel.L_codd_p, p_source=sel.p_src))
    sp.to_csv('step6_codd_source_vs_portal.csv', index=False)
    print('source vs portal: r =', np.corrcoef(sp.b_portal, sp.b_source)[0, 1], 'b ratio', (sp.b_source / sp.b_portal).agg(['min', 'max']).round(4).tolist(),
          'SE ratio', (sp.se_source / sp.se_portal_from_P).agg(['min', 'max']).round(4).tolist(), 'missing source', sp.b_source.isna().sum())

    outcomes = [('S_mvpEU', 'SBP', 'MVP European', 'primary', 425740), ('D_mvpEU', 'DBP', 'MVP European', 'primary', 425743),
                ('H_mvpEU', 'Hypertension', 'MVP European', 'primary', 318398), ('H_fg12', 'Hypertension', 'FinnGen R12', 'secondary', 500264),
                ('S_gh', 'SBP', 'Genes & Health', 'secondary', 18536), ('D_gh', 'DBP', 'Genes & Health', 'secondary', 18536),
                ('S_bbj', 'SBP', 'Biobank Japan', 'secondary', 136597), ('D_bbj', 'DBP', 'Biobank Japan', 'secondary', 136615),
                ('S_mvpAA', 'SBP', 'MVP African American', 'secondary', 119331), ('D_mvpAA', 'DBP', 'MVP African American', 'secondary', 119332),
                ('H_mvpAA', 'Hypertension', 'MVP African American', 'secondary', 73034), ('S_mvpHS', 'SBP', 'MVP Hispanic', 'secondary', 57988),
                ('D_mvpHS', 'DBP', 'MVP Hispanic', 'secondary', 57990), ('H_mvpHS', 'Hypertension', 'MVP Hispanic', 'secondary', 52787)]
    for selcol in ['selected_r2_0.001', 'selected_r2_0.1']:
        iv = fw[fw[selcol]].copy(); f = iv.afEU.fillna(iv.maf).clip(0.001, 0.999)
        iv['bx'] = iv.bx_src; iv['sx'] = iv.sx_src; iv['sx_nAF'] = 1 / np.sqrt(2 * N_CODD * f * (1 - f))
        for code, trait, src, fam, nout in outcomes:
            m = iv.dropna(subset=[code + '_b', code + '_p']); m = m[m[code + '_b'] != 0].copy()
            m['by'] = m[code + '_b']
            m['sy'] = m[code + '_se'] if code + '_se' in m else se_p(m.by, m[code + '_p'])
            r2x = steiger_r2_cont(m.bx, m.sx, N_CODD)
            r2y = steiger_r2_binary(m.by, m.afEU.fillna(m.maf).values, None, None, None) if trait == 'Hypertension' else steiger_r2_cont(m.by, m.sy, nout)
            m['steiger_z'] = steiger_direction(r2x, N_CODD, r2y, nout)[0]
            m['terc_tert_region'] = TT(m).values
            m['palindromic_ambiguous'] = m.palindromic & (np.minimum(m.eur_af_alt, 1 - m.eur_af_alt) > 0.42)
            m.assign(code=code, trait=trait, source=src, family=fam, selection=selcol)[
                ['code', 'trait', 'source', 'family', 'selection', 'varId', 'rsid', 'gene', 'bx', 'sx', 'sx_nAF', 'by', 'sy', 'p_src', code + '_p', 'steiger_z',
                 'terc_tert_region', 'palindromic', 'palindromic_ambiguous', 'panel_status']].rename(columns={code + '_p': 'p_outcome', 'p_src': 'p_exposure'}).to_csv(
                f'mr_v4_inputs/forward_{selcol}_{code}.csv', index=False)

    # ---------------- published Codd MR instrument (Supplementary Table 1, MR = Y)
    st = pd.read_csv('source_data/codd2021_sentinels_parsed.csv', dtype={'Chr': str})
    st = st[(st.MR == 'Y') & (st.Chr != '23')].copy()
    st['pos'] = st.BP_GRCh37.astype(int); st['mhc'] = MHC(st.Chr, st.pos)
    r2 = json.load(open('source_data/kp_v4_round2.json'))['r2']
    rows = []
    for r in st.itertuples():
        rec = None
        for vid, sign in [(r.varA, 1), (r.varB, -1)]:  # varA = chr:pos:A2:A1 (KP alt = A1, the Codd effect allele)
            if r2.get(vid):
                rec = (vid, sign, r2[vid]); break
        src = C[(C.chr == r.Chr) & (C.pos == r.pos)]
        src = src[((src.ea == r.A1) & (src.oa == r.A2)) | ((src.ea == r.A2) & (src.oa == r.A1))]
        if src.empty: continue
        s = src.iloc[0]; bx = s.beta_c * (1 if s.ea == r.A1 else -1)
        out = dict(rsid=r._3, chr=r.Chr, pos=r.pos, gene=r.GENE, effect_allele=r.A1, other_allele=r.A2, bx=bx, sx=s.se_c, p_exposure=s.p_c, mhc=r.mhc,
                   portal_varId=rec[0] if rec else None)
        if rec:
            for ds, ph, b, p, n in rec[2]:
                if ds == 'Verma2024_MVPTraits_EU' and ph in ('SBP', 'DBP', 'HYPERTENSION'):
                    out[f'{ph}_b'] = b * rec[1]; out[f'{ph}_p'] = p; out[f'{ph}_n'] = n
        rows.append(out)
    P = pd.DataFrame(rows); P.to_csv('step6_codd_published_instrument_lookup.csv', index=False)
    print('published MR instrument: autosomal', len(st), 'with source data', len(P), 'MHC', int(P.mhc.sum()), 'with MVP', P.SBP_b.notna().sum())
    for code, ph, trait in [('S_mvpEU', 'SBP', 'SBP'), ('D_mvpEU', 'DBP', 'DBP'), ('H_mvpEU', 'HYPERTENSION', 'Hypertension')]:
        m = P[~P.mhc].dropna(subset=[f'{ph}_b']).copy(); m = m[m[f'{ph}_b'] != 0]
        m['by'] = m[f'{ph}_b']; m['sy'] = se_p(m.by, m[f'{ph}_p']); m['varId'] = m.portal_varId
        m['code'] = code; m['trait'] = trait; m['source'] = 'MVP European'; m['family'] = 'sensitivity'; m['selection'] = 'Codd 2021 published MR instrument'
        m['sx_nAF'] = m.sx; m['steiger_z'] = np.nan; m['terc_tert_region'] = TT(m.assign(chr=m.chr.astype(str))).values
        m['palindromic'] = m.apply(lambda x: {x.effect_allele, x.other_allele} in ({'A', 'T'}, {'C', 'G'}), axis=1); m['palindromic_ambiguous'] = False; m['panel_status'] = 'not clumped (published set)'
        m['p_outcome'] = m[f'{ph}_p']
        m[['code', 'trait', 'source', 'family', 'selection', 'varId', 'rsid', 'gene', 'bx', 'sx', 'sx_nAF', 'by', 'sy', 'p_exposure', 'p_outcome', 'steiger_z',
           'terc_tert_region', 'palindromic', 'palindromic_ambiguous', 'panel_status']].to_csv(f'mr_v4_inputs/forward_codd_published_{code}.csv', index=False)
        print(code, 'published-set instruments with MVP data:', len(m))
