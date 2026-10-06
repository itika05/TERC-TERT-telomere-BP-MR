"""v9 (review minor 3 follow-up): pattern of the P-derived SE in MVP European hypertension records.
For a variant with allele frequency f, the SE of a per-allele effect is approximately c / sqrt(2 f (1 - f)) when all variants share
one sample size and imputation quality, so s = SE * sqrt(2 f (1 - f)) should not depend on |z|. This is checked for the MVP SBP, DBP and
hypertension records in the TERT +/-250 kb window and for the 116 primary instruments, and, as a reference, for FinnGen R12 I9_HYPTENS,
whose files report both SE and P. A sensitivity analysis then repeats the MVP hypertension IVW with an allele-frequency model SE
(c = median of s over the 116 instruments). Output: step11_hypertension_se_pattern.csv (S79), step11_hypertension_se_ivw.csv."""
import numpy as np, pandas as pd
from scipy.stats import norm, spearmanr, t as tdist
rows = []
def summ(label, s, z):
    rho, p = spearmanr(s, z); qb = pd.cut(z, [0, 1, 2, 3, np.inf])
    med = pd.Series(np.asarray(s)).groupby(np.asarray(qb.astype(str))).median()
    rows.append(dict(dataset=label, n=len(s), spearman_rho_s_vs_absz=rho, spearman_p=p, **{f'median_s_absz_{k}': v for k, v in med.items()}))
for tr in ['SBP', 'DBP', 'HYPERTENSION']:
    T = pd.read_csv(f'coloc_v4_inputs/TERT_250kb_{tr}.csv'); T = T[(T.se_out > 0) & T.eur_af_alt.between(0.05, 0.95)]
    summ(f'MVP EU {tr}, TERT +/-250 kb window (MAF >= 5%)', T.se_out * np.sqrt(2 * T.eur_af_alt * (1 - T.eur_af_alt)), norm.isf(T.p_out / 2))
F = pd.read_csv('source_data/finngen_R12_I9_HYPTENS_subset_hg38.tsv', sep='\t'); F = F[(F.pval > 1e-300) & (F.pval < 0.999) & (F.beta != 0) & F.af_alt.between(0.05, 0.95)]
summ('FinnGen R12 I9_HYPTENS, all variants in the held extract, P-derived SE (MAF >= 5%)', F.beta.abs() / norm.isf(F.pval / 2) * np.sqrt(2 * F.af_alt * (1 - F.af_alt)), norm.isf(F.pval / 2))
rows.append(dict(dataset='FinnGen R12 I9_HYPTENS: P-derived SE / reported SE, median (min to max)', n=len(F),
                 spearman_rho_s_vs_absz=np.nan, note=f'{np.median(F.beta.abs() / norm.isf(F.pval / 2) / F.sebeta):.6f} ({(F.beta.abs() / norm.isf(F.pval / 2) / F.sebeta).min():.4f} to {(F.beta.abs() / norm.isf(F.pval / 2) / F.sebeta).max():.4f})'))
AF = pd.read_csv('step8_finngen_harmonisation_audit_v6.csv')[['varId', 'eur_af_kp_alt']]
def ivw(b, se, by, sy):
    r = by / b; w = (b / sy) ** 2; bh = (w * r).sum() / w.sum(); k = len(b); Q = (w * (r - bh) ** 2).sum()
    s = np.sqrt(max(1, Q / (k - 1)) / w.sum()); q = tdist.ppf(0.975, k - 1); return bh, s, bh - q * s, bh + q * s, 2 * tdist.sf(abs(bh / s), k - 1), k
IV = []
for code in ['S_mvpEU', 'D_mvpEU', 'H_mvpEU']:
    d = pd.read_csv(f'mr_v4_inputs/forward_selected_r2_0.001_{code}.csv').merge(AF, on='varId', how='left'); assert d.eur_af_kp_alt.notna().all() and len(d) == 116
    pq = 2 * d.eur_af_kp_alt * (1 - d.eur_af_kp_alt); s = d.sy * np.sqrt(pq); z = (d.by / d.sy).abs()
    summ(f'MVP EU {code[0]} instruments (k = 116)', s, z)
    base = ivw(d.bx, d.sx, d.by, d.sy); IV.append(dict(code=code, outcome_se='P-derived (primary)', b=base[0], se=base[1], lo=base[2], hi=base[3], p=base[4], k=base[5]))
    sm = np.median(s) / np.sqrt(pq); alt = ivw(d.bx, d.sx, d.by, sm)
    IV.append(dict(code=code, outcome_se='allele-frequency model: median(SE*sqrt(2pq)) / sqrt(2pq)', b=alt[0], se=alt[1], lo=alt[2], hi=alt[3], p=alt[4], k=alt[5]))
O = pd.DataFrame(rows); O.to_csv('step11_hypertension_se_pattern.csv', index=False); print(O.to_string())
I = pd.DataFrame(IV); I.to_csv('step11_hypertension_se_ivw.csv', index=False); print(I.to_string())
ref = pd.read_csv('step8_mr_results_v6.csv'); ref = ref[(ref.analysis == 'r2<0.001 (primary)') & (ref.method == 'IVW (MRE, t)')].set_index('code')
for r in I[I.outcome_se.str.startswith('P-derived')].itertuples(): assert abs(r.b - ref.loc[r.code].b) < 1e-9 and abs(r.se - ref.loc[r.code].se) < 1e-9, r
