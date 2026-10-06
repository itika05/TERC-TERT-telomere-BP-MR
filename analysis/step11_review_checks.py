"""v9: data checks requested by the pre-submission review (all from files already in the release; nothing new retrieved).
Outputs
  step11_steiger_reverse_removed.csv   (S76) reverse-MR instruments removed by Steiger filtering, with locus annotation from Codd et al.
                                        Supplementary Table 1 and the effect of removing each one (review M7)
  step11_rs2736100_hypertension_se.csv (S77) the MVP European hypertension record of rs2736100 and its P-derived SE in the context of
                                        the other TERT-window records (review minor 3)
  step11_review_value_checks.csv       (S78) exact values behind review queries (M4 Egger vs weighted median; minor 1, 2, 8, 11)
"""
import json, numpy as np, pandas as pd
from scipy import stats
from scipy.stats import norm
rows = []
def chk(item, quantity, value, source): rows.append(dict(review_item=item, quantity=quantity, value=value, source=source))

# ---- M7: Steiger-removed reverse instruments ------------------------------------------------------------------------------------
S13 = pd.read_csv('step9_reverse_inputs_S13.csv')
raw = pd.read_csv('source_data/codd2021_supp_table1_sentinels.tsv', sep='\t', skiprows=3, header=None)
# Codd Supplementary Table 1: col0 locus number (blank for further signals at the same locus), col1 chr, col2 sentinel, col3 position, col4 locus name
raw = raw[[0, 1, 2, 3, 4]].copy(); raw.columns = ['locus_no', 'chr', 'sentinel', 'pos', 'locus_name']
raw = raw[raw.chr.notna() & raw.pos.notna()].copy()
raw['locus_no'] = raw.locus_no.ffill(); raw['locus_name'] = raw.groupby('locus_no').locus_name.transform(lambda s: s.dropna().iloc[0] if s.notna().any() else np.nan)
raw['chr'] = raw.chr.astype(str).str.replace(r'\.0$', '', regex=True); raw['pos'] = raw.pos.astype(float).astype(int)
def ivw(x):
    b = x.by / x.bx; w = (x.bx / x.sy) ** 2; bh = (w * b).sum() / w.sum(); k = len(x)
    Q = (w * (b - bh) ** 2).sum(); se = np.sqrt(max(1, Q / (k - 1)) / w.sum()); return bh, se, 2 * stats.t.sf(abs(bh / se), k - 1), k
R = pd.read_csv('step8_mr_results_v6.csv')
out = []
for code, ex in [('SBP->L_codd', 'SBP'), ('DBP->L_codd', 'DBP')]:
    x = S13[(S13.code == code) & S13.in_custom_primary.astype(bool)].copy(); full = ivw(x)
    ref = R[(R.code == code) & (R.analysis == 'r2<0.001 (primary)') & (R.method == 'IVW (MRE, t)')].iloc[0]
    assert abs(full[0] - ref.b) < 1e-10 and full[3] == ref.nsnp, (code, full, ref.b)
    kept = x[x.in_steiger_filtered.astype(bool)]; st = ivw(kept)
    refs = R[(R.code == code) & (R.analysis == 'r2<0.001, Steiger-filtered') & (R.method == 'IVW (MRE, t)')].iloc[0]
    assert abs(st[0] - refs.b) < 1e-10 and st[3] == refs.nsnp
    wsum = ((x.bx / x.sy) ** 2).sum()
    for i, r in x[~x.in_steiger_filtered.astype(bool)].iterrows():
        ch, pos = r.varId.split(':')[:2]; pos = int(pos)
        c = raw[raw.chr == ch].copy(); c['d'] = (c.pos - pos).abs(); n = c.sort_values('d').iloc[0]
        d1 = ivw(x.drop(i))
        out.append(dict(reverse_exposure=ex, rsid=r.rsid, varId=r.varId, exposure_beta=r.bx, exposure_se=r.sx, ltl_beta=r.by, ltl_se=r.sy, ltl_source_p=r.src_p,
                        steiger_z=r.steiger_z, ratio_estimate=r.by / r.bx, ivw_weight_share=(r.bx / r.sy) ** 2 / wsum,
                        nearest_codd_sentinel=n.sentinel, codd_locus=n.locus_name, distance_bp=int(n.d),
                        ivw_full=full[0], ivw_without_this_variant=d1[0], p_without_this_variant=d1[2], ivw_without_all_removed=st[0], p_without_all_removed=st[2], k_full=full[3]))
O = pd.DataFrame(out); O.to_csv('step11_steiger_reverse_removed.csv', index=False); print(O[['reverse_exposure', 'rsid', 'codd_locus', 'distance_bp', 'steiger_z', 'ratio_estimate', 'ivw_without_this_variant']].to_string())

# ---- minor 3: rs2736100 MVP hypertension SE --------------------------------------------------------------------------------------
T = pd.read_csv('coloc_v4_inputs/TERT_250kb_HYPERTENSION.csv'); T = T[(T.se_out > 0) & T.eur_af_alt.between(0.05, 0.95)].copy()
T['se_x_sqrt_2pq'] = T.se_out * np.sqrt(2 * T.eur_af_alt * (1 - T.eur_af_alt))     # ~ 1/sqrt(N_eff) if records shared one N and imputation quality
med = T.se_x_sqrt_2pq.median(); T['se_expected_from_window'] = med / np.sqrt(2 * T.eur_af_alt * (1 - T.eur_af_alt))
T['se_ratio_observed_to_expected'] = T.se_out / T.se_expected_from_window; T['percentile_of_ratio'] = T.se_ratio_observed_to_expected.rank(pct=True) * 100
r = T[T.rsid == 'rs2736100'].iloc[0]
# precision: the P-derived SE uses beta and P as printed in the portal record; perturb each by half a unit in its last printed digit
b, p = -0.004988, 0.4809
se0 = abs(b) / norm.isf(p / 2); se_hi = abs(b - 0.0000005) / norm.isf((p - 0.00005) / 2); se_lo = abs(b + 0.0000005) / norm.isf((p + 0.00005) / 2)
sel = T[T.rsid.isin(['rs7726159', 'rs7705526', 'rs2736100', 'rs2853677', 'rs2736098'])][['rsid', 'varId', 'eur_af_alt', 'b_out', 'se_out', 'p_out', 'n_out', 'se_expected_from_window', 'se_ratio_observed_to_expected', 'percentile_of_ratio']]
sel = sel.assign(note='')
sel.loc[sel.rsid == 'rs2736100', 'note'] = f'portal record beta {b}, P {p}, N {r.n_out:.2f}; SE from P {se0:.6f} (rounding range {min(se_lo, se_hi):.6f} to {max(se_lo, se_hi):.6f})'
sel.to_csv('step11_rs2736100_hypertension_se.csv', index=False); print(sel.to_string())
chk('Minor 3', 'rs2736100 MVP EU hypertension: P-derived SE / SE expected from window median and allele frequency', round(float(r.se_ratio_observed_to_expected), 3), 'S77')
chk('Minor 3', 'rs2736100: percentile of that ratio among TERT-window records (MAF >= 5%)', round(float(r.percentile_of_ratio), 1), 'S77')
chk('Minor 3', 'rs2736100: SE range from rounding of printed beta and P', f'{min(se_lo, se_hi):.6f} to {max(se_lo, se_hi):.6f}', 'S77')
chk('Minor 3', 'TERT-window hypertension records with ratio >= rs2736100 ratio', int((T.se_ratio_observed_to_expected >= r.se_ratio_observed_to_expected).sum()), 'S77')
chk('Minor 3', 'TERT-window hypertension records (MAF >= 5%)', len(T), 'S77')

# ---- M4: exact estimator values ---------------------------------------------------------------------------------------------------
for m in ['IVW (MRE, t)', 'MR-Egger (t)', 'Weighted median', 'Weighted mode (MBE)']:
    x = R[(R.code == 'S_mvpEU') & (R.analysis == 'r2<0.001 (primary)') & (R.method == m)].iloc[0]
    chk('Major 4', f'SBP {m}: estimate (SE)', f'{x.b:.6f} ({x.se:.6f})', 'S11')

# ---- minor 1: published instrument accounting --------------------------------------------------------------------------------------
PUB = pd.read_csv('step8_published_vs_portal_pool.csv')
chk('Minor 1', 'published Codd instrument variants', len(PUB), 'S65')
chk('Minor 1', 'identical to one of our instruments', int(PUB.same_variant.sum()), 'S65')
chk('Minor 1', 'within 500 kb of an instrument (includes the identical ones)', int(PUB.within_500kb.sum()), 'S65')
chk('Minor 1', 'within 500 kb but not identical', int((PUB.within_500kb & ~PUB.same_variant).sum()), 'S65')
chk('Minor 1', 'farther than 500 kb', int((~PUB.within_500kb).sum()), 'S65')

# ---- minor 2: Figure 2c k label ---------------------------------------------------------------------------------------------------
for c in ['S_mvpEU', 'D_mvpEU', 'H_mvpEU']:
    chk('Minor 2', f'{c} primary IVW k (value written in Figure 2c)', int(R[(R.code == c) & (R.analysis == 'r2<0.001 (primary)') & (R.method == 'IVW (MRE, t)')].iloc[0].nsnp), 'S11; fig2_mr_v7.py')

# ---- minor 8: Hispanic rs10936599 SBP record ----------------------------------------------------------------------------------------
A = pd.read_csv('kp_snp_lookup_by_ancestry_checked.csv'); h = A[(A.rs == 'rs10936599') & (A.trait == 'SBP') & (A.ancestry == 'HS')].iloc[0]
chk('Minor 8', 'Hispanic rs10936599 SBP: beta, SE, P as reported', f'{h.beta_alt}, {h.se}, {h.p}', 'kp_snp_lookup_by_ancestry_checked.csv')
chk('Minor 8', 'z from beta/SE vs z from P', f'{h.z_se:.2f} vs {abs(h.z_p):.2f}', 'same')
chk('Minor 8', 'two-sided P implied by beta/SE', f'{2 * norm.sf(abs(h.z_se)):.2f}', 'computed')
chk('Minor 8', 'SE implied by beta and P vs reported SE (exclusion rule: implied < 0.5 x reported)', f'{h.se_p:.5f} vs {h.se} (ratio {h.se_p / h.se:.2f})', 'step2_lookup.py')

# ---- minor 11 / M5: effective N of binary MVP records ----------------------------------------------------------------------------
KD = pd.read_csv('kp_snp_lookup_by_dataset.csv'); e = KD[(KD.varId_b37 == '5:1286516:C:A') & (KD.trait == 'EssentialHYPERTENSION') & (KD.dataset == 'Verma2024_MVPTraits_EU')].iloc[0]
neff = 4 / (1 / 320429 + 1 / 107275)
chk('Minor 11', 'MVP EU EssentialHYPERTENSION portal N (rs2736100 record)', f'{e.n:.2f}', 'kp_snp_lookup_by_dataset.csv')
chk('Minor 11', '4/(1/cases + 1/controls) for dbGaP pha005550 Phe_401_1.EUR (320,429 cases, 107,275 controls)', f'{neff:.2f}', 'dbGaP pha005550 description')
h_ = KD[(KD.varId_b37 == '5:1286516:C:A') & (KD.trait == 'HYPERTENSION') & (KD.dataset == 'Verma2024_MVPTraits_EU')].iloc[0]
chk('Minor 11', 'MVP EU HYPERTENSION portal N used in this study (rs2736100 record)', f'{h_.n:.2f}', 'kp_snp_lookup_by_dataset.csv')
pd.DataFrame(rows).to_csv('step11_review_value_checks.csv', index=False); print(pd.DataFrame(rows).to_string())
