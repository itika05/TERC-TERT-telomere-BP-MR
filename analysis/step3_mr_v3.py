"""Step 3, version 3: instrument selection with LD computed from 1000 Genomes phase 3 EUR phased haplotypes (503 individuals),
PLINK-style clumping at r2 < 0.001 within 10 Mb in both MR directions, full variant and pair ledgers, and export of harmonised inputs
for estimation with established R packages (mr_v3.R).

Status rules (reviewer request: never code missing LD as zero):
  * a candidate is usable only if it is matched to a panel record by position and alleles, polymorphic in EUR, and autosomal;
  * unusable candidates (absent from panel, allele mismatch, monomorphic, chrX) are EXCLUDED from the primary instrument set
    because their independence cannot be established; they are listed in the variant ledger;
  * every candidate pair within 10 Mb is listed in the pair ledger with r2 (known) or the reason it is unknown.
SE: dataset-level portal records give beta, P and N only (no SE). SE = |beta|/z(P). No forward instrument reaches the P floor
(min Codd P = 1.3e-293); for the 2 reverse-outcome Codd records at the floor, SE is taken from N and allele frequency
(1/sqrt(2 N f(1-f)), standardised trait) and flagged. Sensitivity: Codd SE from N and allele frequency for every instrument."""
import numpy as np, pandas as pd
from scipy.stats import norm
from ld_ref import load, match, hap_matrix, corr
from mr_lib import steiger_r2_cont, steiger_r2_binary, steiger_direction

REF = load('ldref/eur_instruments.json')
N_CODD = 464716
se_p = lambda b, p: np.abs(b) / norm.isf(np.clip(p, 1e-300, 1) / 2)
se_n = lambda f, n: 1 / np.sqrt(2 * n * f * (1 - f))


def clump(cand, pcol, r2max=0.001, win=10_000_000):
    """cand: DataFrame with varId, chr, pos, usable flag, idx/sign (panel), sorted by pcol ascending. Returns kept varIds, drop log, pair ledger."""
    cand = cand.sort_values(pcol).reset_index(drop=True)
    U = cand[cand.usable]
    H = hap_matrix(U, REF); Rm = corr(H); pos = {v: i for i, v in enumerate(U.varId)}
    pairs = []
    for c, g in cand.groupby('chr'):
        g = g.sort_values('pos')
        v, p, u, st = g.varId.values, g.pos.values, g.usable.values, g.panel_status.values
        for i in range(len(v)):
            for j in range(i + 1, len(v)):
                if p[j] - p[i] >= win: break
                if u[i] and u[j]:
                    r = Rm[pos[v[i]], pos[v[j]]]; pairs.append(dict(a=v[i], b=v[j], dist=p[j] - p[i], r=r, r2=r * r, ld_status='computed (both variants usable)'))
                else:
                    why = '; '.join(sorted(set(s for s, uu in [(st[i], u[i]), (st[j], u[j])] if not uu)))
                    pairs.append(dict(a=v[i], b=v[j], dist=p[j] - p[i], r=np.nan, r2=np.nan, ld_status='unknown: ' + why))
    kept, drops = [], []
    for row in cand.itertuples():
        if not row.usable:
            drops.append(dict(varId=row.varId, reason='excluded: ' + row.panel_status, by='')); continue
        hit = None
        for k in kept:
            kr = cand.loc[cand.varId == k].iloc[0]
            if kr.chr == row.chr and abs(kr.pos - row.pos) < win:
                r2 = Rm[pos[k], pos[row.varId]] ** 2
                if r2 >= r2max: hit = (k, r2); break
        if hit: drops.append(dict(varId=row.varId, reason=f'r2 = {hit[1]:.4f} >= {r2max} with stronger instrument', by=hit[0]))
        else: kept.append(row.varId)
    return kept, pd.DataFrame(drops), pd.DataFrame(pairs), Rm, pos


def prep(cand):
    t = cand.varId.str.split(':', expand=True); cand = cand.assign(chr=t[0], pos=t[1].astype(int))
    m = match(cand.varId, REF); cand = cand.merge(m, on='varId')
    cand.loc[cand.chr == 'X', 'panel_status'] = 'chrX (haploid males; not assessed)'
    cand['usable'] = cand.panel_status.str.startswith('matched') & (cand.chr != 'X')
    return cand


if __name__ == '__main__':
    # ------------------------------ forward: LTL -> BP
    d = pd.read_csv('kp_ltl_leads_wide.csv'); d = d[d.varId.str.count(':') == 3]
    fw = d[(d.L_codd_p < 5e-8) & d.rsid.notna() & (np.minimum(d.maf, 1 - d.maf) >= 0.01)].copy()
    fw = prep(fw)
    kept, drops, pairs, Rm, posmap = clump(fw, 'L_codd_p')
    fw['selected_r2_0.001'] = fw.varId.isin(kept)
    k01, _, _, _, _ = clump(fw, 'L_codd_p', r2max=0.1)
    fw['selected_r2_0.1'] = fw.varId.isin(k01)
    fw.drop(columns=['idx']).to_csv('step3v3_forward_variant_ledger.csv', index=False)
    drops.to_csv('step3v3_forward_clump_log.csv', index=False); pairs.to_csv('step3v3_forward_pair_ledger.csv', index=False)
    print('forward candidates', len(fw), fw.panel_status.value_counts().to_dict(), 'kept r2<0.001:', len(kept), 'kept r2<0.1:', len(k01))
    print('pairs', len(pairs), pairs.ld_status.value_counts().to_dict(), 'r2>=0.001 among computed:', (pairs.r2 >= 0.001).sum())
    # correlation matrix for the r2<0.1 set (covariance-aware IVW)
    S01 = fw[fw['selected_r2_0.1']].sort_values('L_codd_p')
    idx = [posmap[v] for v in S01.varId]; pd.DataFrame(Rm[np.ix_(idx, idx)], index=S01.varId, columns=S01.varId).to_csv('mr_v3_inputs/forward_r2_0.1_correlation.csv') if False else None
    import os; os.makedirs('mr_v3_inputs', exist_ok=True)
    pd.DataFrame(Rm[np.ix_(idx, idx)], index=S01.varId, columns=S01.varId).to_csv('mr_v3_inputs/forward_r2_0.1_correlation.csv')

    outcomes = [('S_mvpEU', 'SBP', 'MVP European', 'primary', 425740), ('D_mvpEU', 'DBP', 'MVP European', 'primary', 425743),
                ('H_mvpEU', 'Hypertension', 'MVP European', 'primary', 318398), ('H_fg', 'Hypertension', 'FinnGen', 'secondary', 52313),
                ('S_gh', 'SBP', 'Genes & Health', 'secondary', 18536), ('D_gh', 'DBP', 'Genes & Health', 'secondary', 18536),
                ('S_bbj', 'SBP', 'Biobank Japan', 'secondary', 136597), ('D_bbj', 'DBP', 'Biobank Japan', 'secondary', 136615),
                ('S_mvpAA', 'SBP', 'MVP African American', 'secondary', 119331), ('D_mvpAA', 'DBP', 'MVP African American', 'secondary', 119332),
                ('H_mvpAA', 'Hypertension', 'MVP African American', 'secondary', 73034), ('S_mvpHS', 'SBP', 'MVP Hispanic', 'secondary', 57988),
                ('D_mvpHS', 'DBP', 'MVP Hispanic', 'secondary', 57990), ('H_mvpHS', 'Hypertension', 'MVP Hispanic', 'secondary', 52787)]
    TT = lambda x: ((x.chr == '3') & x.pos.between(168_982_000, 169_983_000)) | ((x.chr == '5') & x.pos.between(753_000, 1_795_000))
    rows = []
    for sel in ['selected_r2_0.001', 'selected_r2_0.1']:
        iv = fw[fw[sel]].copy()
        f = iv.afEU.fillna(iv.maf).clip(0.001, 0.999)
        iv['bx'] = iv.L_codd_b; iv['sx'] = se_p(iv.L_codd_b, iv.L_codd_p); iv['sx_nAF'] = se_n(f, N_CODD)
        for code, trait, src, fam, nout in outcomes:
            m = iv.dropna(subset=[code + '_b', code + '_p']); m = m[m[code + '_b'] != 0].copy()
            m['by'] = m[code + '_b']; m['sy'] = se_p(m.by, m[code + '_p'])
            r2x = steiger_r2_cont(m.bx, m.sx, N_CODD)
            r2y = steiger_r2_binary(m.by, m.afEU.fillna(m.maf).values, None, None, None) if trait == 'Hypertension' else steiger_r2_cont(m.by, m.sy, nout)
            m['steiger_z'] = steiger_direction(r2x, N_CODD, r2y, nout)[0]
            m['terc_tert_region'] = TT(m).values
            m['palindromic_ambiguous'] = m.palindromic & (np.minimum(m.eur_af_alt, 1 - m.eur_af_alt) > 0.42)
            m.assign(code=code, trait=trait, source=src, family=fam, selection=sel)[
                ['code', 'trait', 'source', 'family', 'selection', 'varId', 'rsid', 'gene', 'bx', 'sx', 'sx_nAF', 'by', 'sy', 'L_codd_p', code + '_p', 'steiger_z',
                 'terc_tert_region', 'palindromic', 'palindromic_ambiguous', 'panel_status']].rename(columns={code + '_p': 'p_outcome', 'L_codd_p': 'p_exposure'}).to_csv(
                f'mr_v3_inputs/forward_{sel}_{code}.csv', index=False)
            rows.append(dict(direction='forward', selection=sel, code=code, n=len(m)))

    # ------------------------------ reverse: BP -> LTL
    r = pd.read_csv('kp_bp_leads_reverse_wide.csv').merge(pd.read_csv('kp_bp_leads_rsid.csv'), on='varId', how='left'); r = r[r.varId.str.count(':') == 3]
    t = r.varId.str.split(':', expand=True); r = r.assign(chr=t[0], pos=t[1].astype(int))
    r = r[~((r.chr == '6') & r.pos.between(25_000_000, 34_000_000))]
    r = r[r.rsid.notna() & (np.minimum(r.maf, 1 - r.maf) >= 0.01)].drop(columns=['chr', 'pos'])
    led = []
    for ex, exname in [('S_mvp', 'SBP'), ('D_mvp', 'DBP'), ('H_mvp', 'Hypertension')]:
        c = prep(r[r[ex + '_p'] < 5e-8].copy())
        kept, drops, pairs, _, _ = clump(c, ex + '_p')
        c['selected_r2_0.001'] = c.varId.isin(kept)
        c.drop(columns=['idx']).assign(exposure=exname).to_csv(f'step3v3_reverse_{exname}_variant_ledger.csv', index=False)
        drops.to_csv(f'step3v3_reverse_{exname}_clump_log.csv', index=False); pairs.to_csv(f'step3v3_reverse_{exname}_pair_ledger.csv', index=False)
        led.append(dict(exposure=exname, candidates=len(c), usable=int(c.usable.sum()), **{f'status: {k}': v for k, v in c.panel_status.value_counts().items()},
                        pairs_within_10Mb=len(pairs), pairs_computed=int(pairs.r2.notna().sum()), pairs_unknown=int(pairs.r2.isna().sum()),
                        pairs_r2_ge_0_001=int((pairs.r2 >= 0.001).sum()), selected=len(kept)))
        iv = c[c['selected_r2_0.001']].copy()
        for oc, ocname, nout in [('L_codd', 'LTL (Codd 2021)', N_CODD), ('L_nkSA', 'LTL South Asian (Nakao 2026)', 11277)]:
            m = iv.dropna(subset=[oc + '_b']); m = m[(m[oc + '_b'] != 0)].copy()
            m['bx'] = m[ex + '_b']; m['sx'] = se_p(m.bx, m[ex + '_p']); m['by'] = m[oc + '_b']; m['sy'] = se_p(m.by, m[oc + '_p'])
            f = m.eur_af_alt.clip(0.001, 0.999)
            floor = m[oc + '_p'] <= 1e-300; m['outcome_p_floor'] = floor
            m.loc[floor, 'sy'] = se_n(f[floor], nout)
            p1 = ~np.isfinite(m.sy) | (m.sy <= 0); m['outcome_p_one'] = p1   # P recorded as 1 with non-zero beta: SE not recoverable from P
            m.loc[p1, 'sy'] = se_n(f[p1], nout)
            nx = m[ex + '_n'].median()
            if exname != 'Hypertension':
                m['steiger_z'] = steiger_direction(steiger_r2_cont(m.bx, m.sx, nx), nx, steiger_r2_cont(m.by, m.sy, nout), nout)[0]
            else:
                m['steiger_z'] = np.nan
            m.assign(exposure=exname, outcome=ocname)[['exposure', 'outcome', 'varId', 'rsid', 'bx', 'sx', 'by', 'sy', 'outcome_p_floor', 'outcome_p_one', 'steiger_z', 'panel_status']].to_csv(
                f'mr_v3_inputs/reverse_{exname}_{oc}.csv', index=False)
    L = pd.DataFrame(led); L.to_csv('step3v3_reverse_ledger.csv', index=False); print(L.T)
    fl = dict(direction='forward', candidates=len(fw), usable=int(fw.usable.sum()), **{f'status: {k}': v for k, v in fw.panel_status.value_counts().items()},
              pairs_within_10Mb=len(pairs := pd.read_csv('step3v3_forward_pair_ledger.csv')), pairs_computed=int(pairs.r2.notna().sum()),
              pairs_unknown=int(pairs.r2.isna().sum()), pairs_r2_ge_0_001=int((pairs.r2 >= 0.001).sum()), selected=int(fw['selected_r2_0.001'].sum()),
              selected_r2_0_1=int(fw['selected_r2_0.1'].sum()))
    pd.DataFrame([fl]).to_csv('step3v3_forward_ledger.csv', index=False); print(pd.DataFrame([fl]).T)
