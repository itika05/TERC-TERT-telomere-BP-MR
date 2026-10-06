"""Step 3: two-sample MR, genetically predicted leukocyte telomere length (LTL) -> BP / hypertension.
Exposure: Codd et al. 2021 LTL GWAS (UK Biobank, N=464,716 (GWAS sample size as reported)), per-allele betas as deposited in the CMD Knowledge Portal.
Instrument selection: KP European LTL lead variants with Codd 2021 P < 5e-8; greedy LD pruning by Codd P,
  drop if r2 >= 0.05 with a stronger kept instrument within 10 Mb (Ensembl REST, 1000G phase 3 EUR; Ensembl does not report r2 < 0.05)
  or if within 1 Mb of a kept instrument. Stricter r2 < 0.001 clumping would need a local PLINK reference panel.
Outcomes (no UK Biobank samples): MVP (Verma 2024) European SBP, DBP, hypertension; FinnGen hypertension;
  Genes & Health (Huang 2021/2022, British Pakistani/Bangladeshi) SBP, DBP; Biobank Japan SBP, DBP;
  MVP African-American and Hispanic SBP/DBP/HTN.
Units: BP in the KP's standardised (SD) units per SD of LTL; hypertension as log OR per SD LTL.
Alleles: all betas are for the KP 'alt' allele of the same varId (internal KP harmonisation), so no strand flipping is needed.
SE for every dataset-level estimate is derived from beta and P (|beta| / z)."""
import pandas as pd, numpy as np, json
from scipy.stats import norm
from mr_lib import ivw, egger, weighted_median, weighted_mode, mr_presso, fstat, steiger_r2_cont, steiger_r2_binary, steiger_direction

d = pd.read_csv('kp_ltl_leads_wide.csv')
d = d[d.varId.str.count(':') == 3].copy()
d[['chr', 'pos', 'ref', 'alt']] = d.varId.str.split(':', expand=True)
d['pos'] = d.pos.astype(int)

def se_from_p(b, p):
    z = norm.isf(np.clip(p, 1e-300, 1) / 2)
    return np.abs(b) / z

d['bx'] = d.L_codd_b
d['sx'] = se_from_p(d.L_codd_b, d.L_codd_p)
iv = d[(d.L_codd_p < 5e-8) & d.rsid.notna() & (np.minimum(d.maf, 1 - d.maf) >= 0.01)].sort_values('L_codd_p').copy()

# ---------- LD pruning
ld = pd.read_csv('ld_eur_pairs_v2.csv')  # all 327 candidate pairs within 10 Mb queried; r2 blank = below Ensembl reporting floor (r2 < 0.05)
LD = {}
for a_, b_, r in ld[['a', 'b', 'r2']].itertuples(index=False):
    LD[(a_, b_)] = r; LD[(b_, a_)] = r
kept = []
dropped = []
for row in iv.itertuples():
    ok = True
    for k in kept:
        if k.chr != row.chr or abs(k.pos - row.pos) >= 10_000_000:
            continue
        r2 = LD.get((k.rsid, row.rsid), np.nan)
        if not np.isnan(r2) and r2 >= 0.05:
            ok = False; dropped.append((row.rsid, k.rsid, f'r2={r2:.3f} (1000G EUR)')); break
        if abs(k.pos - row.pos) < 1_000_000:
            ok = False; dropped.append((row.rsid, k.rsid, f'within 1 Mb (r2 {"<0.05" if np.isnan(r2) else f"={r2:.3f}"})')); break
    if ok:
        kept.append(row)
iv = pd.DataFrame(kept)
iv['F'] = fstat(iv.bx, iv.sx)
iv.to_csv('step3_instruments_LTL.csv', index=False)
print('candidates', len(pd.read_csv('kp_ltl_leads_wide.csv')), 'dropped', len(dropped))
pd.DataFrame(dropped, columns=['dropped', 'correlated_with', 'reason']).to_csv('step3_instruments_pruned_out.csv', index=False)
print('instruments after pruning:', len(iv), ' mean F', round(iv.F.mean(), 1), ' min F', round(iv.F.min(), 1))

# ---------- sample sizes (dataset-level n from the KP records for the two study SNPs)
nmap = {'S_mvpEU': 425740, 'D_mvpEU': 425743, 'H_mvpEU': 318398, 'H_fg': 52313, 'S_gh': 18536, 'D_gh': 18536,
        'S_bbj': 136597, 'D_bbj': 136615, 'S_mvpAA': 119331, 'D_mvpAA': 119332, 'H_mvpAA': 73034,
        'S_mvpHS': 57988, 'D_mvpHS': 57990, 'H_mvpHS': 52787}
outcomes = [('S_mvpEU', 'SBP', 'MVP European (Verma 2024)', 'EUR'), ('D_mvpEU', 'DBP', 'MVP European (Verma 2024)', 'EUR'),
            ('H_mvpEU', 'Hypertension', 'MVP European (Verma 2024)', 'EUR'), ('H_fg', 'Hypertension', 'FinnGen', 'EUR'),
            ('S_gh', 'SBP', 'Genes & Health (Huang 2021)', 'SAS'), ('D_gh', 'DBP', 'Genes & Health (Huang 2021)', 'SAS'),
            ('S_bbj', 'SBP', 'Biobank Japan', 'EAS'), ('D_bbj', 'DBP', 'Biobank Japan', 'EAS'),
            ('S_mvpAA', 'SBP', 'MVP African American', 'AFR'), ('D_mvpAA', 'DBP', 'MVP African American', 'AFR'),
            ('H_mvpAA', 'Hypertension', 'MVP African American', 'AFR'),
            ('S_mvpHS', 'SBP', 'MVP Hispanic', 'AMR'), ('D_mvpHS', 'DBP', 'MVP Hispanic', 'AMR'), ('H_mvpHS', 'Hypertension', 'MVP Hispanic', 'AMR')]
TERC_TERT = lambda x: ((x.chr == '3') & (x.pos.between(168_982_000, 169_983_000))) | ((x.chr == '5') & (x.pos.between(753_000, 1_795_000)))

res = []; loo_all = []; snp_level = []
for code, trait, src, anc in outcomes:
    m = iv.dropna(subset=[code + '_b', code + '_p']).copy()
    m = m[m[code + '_b'] != 0]
    m['by'] = m[code + '_b']; m['sy'] = se_from_p(m.by, m[code + '_p'])
    bx, sx, by, sy = m.bx.values, m.sx.values, m.by.values, m.sy.values
    base = dict(outcome=trait, source=src, ancestry=anc, code=code)
    r_ivw = ivw(bx, sx, by, sy)
    ests = [r_ivw, egger(bx, sx, by, sy), weighted_median(bx, sx, by, sy), weighted_mode(bx, sx, by, sy)]
    pr = mr_presso(bx, sx, by, sy)
    ests.append(pr)
    # Steiger
    r2x = steiger_r2_cont(bx, sx, 464716)
    if trait == 'Hypertension':
        af = m.afEU.fillna(m.maf).values
        r2y = steiger_r2_binary(by, af, None, None, None)
    else:
        r2y = steiger_r2_cont(by, sy, nmap[code])
    z, pz = steiger_direction(r2x, 464716, r2y, nmap[code])
    keep = z > 0
    r_st = ivw(bx[keep], sx[keep], by[keep], sy[keep]); r_st['method'] = 'IVW, Steiger-filtered'
    ests.append(r_st)
    nt = ~TERC_TERT(m).values
    r_nt = ivw(bx[nt], sx[nt], by[nt], sy[nt]); r_nt['method'] = 'IVW, excluding TERC & TERT loci'
    ests.append(r_nt)
    for e in ests:
        x = dict(base); x.update({k: v for k, v in e.items() if k != 'outlier_idx'})
        x['steiger_prop_correct'] = keep.mean()
        res.append(x)
    # leave-one-out (IVW)
    for i in range(len(bx)):
        mm = np.ones(len(bx), bool); mm[i] = False
        e = ivw(bx[mm], sx[mm], by[mm], sy[mm])
        loo_all.append(dict(code=code, left_out=m.rsid.iloc[i], b=e['b'], se=e['se'], p=e['p']))
    for i in range(len(bx)):
        snp_level.append(dict(code=code, rsid=m.rsid.iloc[i], varId=m.varId.iloc[i], gene=m.gene.iloc[i], bx=bx[i], sx=sx[i], by=by[i], sy=sy[i],
                              steiger_z=z[i], presso_outlier=i in pr['outlier_idx']))

R = pd.DataFrame(res)
R['OR_or_b'] = np.where(R.outcome == 'Hypertension', np.exp(R.b), R.b)
R['lo'] = np.where(R.outcome == 'Hypertension', np.exp(R.b - 1.96 * R.se), R.b - 1.96 * R.se)
R['hi'] = np.where(R.outcome == 'Hypertension', np.exp(R.b + 1.96 * R.se), R.b + 1.96 * R.se)
R.to_csv('step3_mr_forward_results.csv', index=False)
pd.DataFrame(loo_all).to_csv('step3_mr_forward_leave_one_out.csv', index=False)
pd.DataFrame(snp_level).to_csv('step3_mr_forward_snp_level.csv', index=False)
pd.set_option('display.width', 250)
print(R[['outcome', 'source', 'method', 'nsnp', 'b', 'se', 'p', 'OR_or_b', 'lo', 'hi', 'Q_p', 'I2', 'intercept', 'intercept_p', 'global_p', 'n_outliers']].round(4).to_string())
