"""Supplementary Figure S6 (v6): verification analyses. (a) MR-PRESSO convergence with the number of simulations (five seeds per B);
(b) custom (phased-haplotype r2) versus PLINK 1.9 clumping in both MR directions, with IVW estimates; (c) P-derived versus reported SE;
(d) MVMR conditional F: reference MVMR::strength_mvmr versus assumed exposure-error correlation, with the custom v5 value; (e) haplotype versus
genotype-dosage LD in the same 503 individuals; (f) FinnGen R12 versus 1000 Genomes FIN allele frequency for the 116 instruments."""
from figstyle import *
import json, numpy as np, pandas as pd
from scipy.stats import norm
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
CV = pd.read_csv('step8_presso_convergence_summary.csv'); PB = pd.read_csv('step8_plink_both_directions_summary.csv'); R = pd.read_csv('step8_mr_results_v6.csv')
CF = pd.read_csv('step8_mvmr_v6_conditionalF.csv'); LDc = pd.read_csv('step8_ld_haplotype_vs_dosage.csv'); FA = pd.read_csv('step8_finngen_harmonisation_audit_v6.csv')
cal = json.load(open('source_data/kp_v4_round2.json'))['cal']; SPV = pd.read_csv('step6_codd_source_vs_portal.csv')
fig = plt.figure(figsize=(W, 10.0)); gs = fig.add_gridspec(3, 2, hspace=0.62, wspace=0.45)
COL = {'S_mvpEU': BLUE, 'D_mvpEU': AQUA, 'H_mvpEU': ORANGE}; LAB = {'S_mvpEU': 'SBP', 'D_mvpEU': 'DBP', 'H_mvpEU': 'Hypertension (log OR)'}
# a
ax = fig.add_subplot(gs[0, 0]); Bs = sorted(CV.B.unique()); xs = np.arange(len(Bs))
for i, code in enumerate(['S_mvpEU', 'D_mvpEU', 'H_mvpEU']):
    g = CV[CV.code == code].set_index('B').loc[Bs]; off = (i - 1) * 0.22
    for k, (B, r) in enumerate(g.iterrows()):
        ax.add_patch(Rectangle((k + off - 0.09, r.b_min), 0.18, max(r.b_max - r.b_min, 0.0015), color=COL[code], alpha=0.85 if r.identical_across_5_seeds else 0.4, lw=0))
        ax.text(k + off, r.b_max + 0.006, f'{int(r.outliers_min)}' + ('' if r.outliers_min == r.outliers_max else f'–{int(r.outliers_max)}'), ha='center', fontsize=4.4, color=COL[code], rotation=90, va='bottom')
ax.set_xticks(xs); ax.set_xticklabels([f'{B:,}' for B in Bs], fontsize=5.4); ax.set_xlabel('simulations (B)'); ax.set_ylabel('outlier-corrected estimate (range over 5 seeds)')
ax.set_ylim(0.08, 0.235); ax.set_xlim(-0.6, len(Bs) - 0.4); ax.legend(handles=[Rectangle((0, 0), 1, 1, color=COL[c], label=LAB[c]) for c in COL] + [Rectangle((0, 0), 1, 1, color=INK2, alpha=0.4, label='outlier set differs between seeds')], fontsize=4.8, loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2)
ax.text(0.99, 0.98, 'numbers: outliers (range over seeds)', transform=ax.transAxes, ha='right', va='top', fontsize=4.8, color=INK2)
ax.set_title('MR-PRESSO convergence (Python re-implementation; agrees with MRPRESSO 1.0 within tolerance on shared data, S70)'); tag(ax, 'a')
# b
ax = fig.add_subplot(gs[0, 1]); dirs = ['forward', 'SBP', 'DBP', 'Hypertension']
for i, d in enumerate(dirs):
    r = PB[PB.direction == d].iloc[0]; tot = r.custom_only + r.shared + r.plink_only; y = -i
    for left, wdt, col, lab in [(0, r.custom_only, BLUE, 'custom only'), (r.custom_only, r.shared, VIOLET, 'shared'), (r.custom_only + r.shared, r.plink_only, ORANGE, 'PLINK only')]:
        ax.add_patch(Rectangle((left / tot, y - 0.3), wdt / tot, 0.6, color=col, alpha=0.75, lw=0))
    ax.text(r.custom_only / tot + r.shared / tot / 2, y, f'{int(r.shared)} shared · Jaccard {r.jaccard:.2f}', ha='center', va='center', fontsize=5, color='white', fontweight='bold')
    ax.text(-0.02, y, {'forward': 'LTL → BP'}.get(d, f'{d} → LTL'), ha='right', va='center', fontsize=5.6)
    ax.text(1.02, y + 0.12, f'custom {int(r.custom_kept)} / PLINK {int(r.plink_kept)}', ha='left', va='center', fontsize=4.8, color=INK2)
    code = {'forward': 'S_mvpEU', 'SBP': 'SBP->L_codd', 'DBP': 'DBP->L_codd', 'Hypertension': 'Hypertension->L_codd'}[d]
    e1 = R[(R.code == code) & (R.analysis == 'r2<0.001 (primary)') & (R.method == 'IVW (MRE, t)')].iloc[0]; e2 = R[(R.code == code) & (R.analysis == 'PLINK 1.9 clumping') & (R.method == 'IVW (MRE, t)')].iloc[0]
    ax.text(1.02, y - 0.18, f'IVW {e1.b:+.3f} vs {e2.b:+.3f}'.replace('-', '−'), ha='left', va='center', fontsize=4.8, color=INK2)
ax.set_xlim(0, 1); ax.set_ylim(-3.6, 0.6); ax.axis('off')
ax.legend(handles=[Rectangle((0, 0), 1, 1, color=c, alpha=0.75, label=l) for c, l in [(BLUE, 'custom only'), (VIOLET, 'shared'), (ORANGE, 'PLINK only')]], fontsize=4.8, loc='lower center', ncol=3, bbox_to_anchor=(0.5, -0.14))
ax.set_title('Clumping: phased-haplotype r² vs PLINK 1.9 --clump (both directions)', pad=4); tag(ax, 'b', x=-0.2)
# c
ax = fig.add_subplot(gs[1, 0]); names = []
for i, (k, recs) in enumerate(cal.items()):
    X = pd.DataFrame(recs).iloc[:, :4]; X.columns = ['v', 'b', 'se', 'p']; X = X[(X.p > 1e-300) & (X.b != 0)]
    rt = 100 * (X.b.abs() / norm.isf(X.p / 2) / X.se - 1); q = np.percentile(rt, [0, 2.5, 25, 50, 75, 97.5, 100])
    ax.add_patch(Rectangle((q[1], i - 0.3), q[5] - q[1], 0.6, color=SEQ[2], lw=0)); ax.plot([q[0], q[6]], [i, i], color=INK2, lw=0.5); ax.plot([q[3]] * 2, [i - 0.34, i + 0.34], color=INK, lw=1)
    names.append(k.replace('_MVPTraits_EU', ' MVP').replace('_bp_eu', '').replace('_BP_EU', '').replace('_BPtraits_EU', '').replace('GWAS_UKBiobank409k_eu', 'UKB 409k').replace('|', ' '))
i += 1; rt = 100 * (SPV.se_portal_from_P / SPV.se_source - 1).dropna(); q = np.percentile(rt, [2.5, 50, 97.5])
ax.add_patch(Rectangle((q[0], i - 0.3), q[2] - q[0], 0.6, color='#cdc4e8', lw=0)); ax.plot([q[1]] * 2, [i - 0.34, i + 0.34], color=INK, lw=1); names.append('Codd LTL vs source')
ax.set_yticks(range(len(names))); ax.set_yticklabels(names, fontsize=5.2); ax.set_ylim(len(names) - 0.4, -0.6); ax.axvline(0, color=INK2, lw=0.5)
ax.set_xlabel('SE from P [scipy norm.isf(P/2)] vs reported SE (% difference); band 95% of records'); ax.set_xlim(-4.5, 1.5); clean(ax, left=False); ax.tick_params(axis='y', length=0)
ax.set_title('Standard errors reconstructed from P'); tag(ax, 'c')
# d
ax = fig.add_subplot(gs[1, 1])
for e, col, lab in [('LTL', VIOLET, 'LTL'), ('BMI', YELLOW, 'BMI'), ('LYM', AQUA, 'lymphocyte count')]:
    z = CF[CF.exposure == e].sort_values('rho'); ax.plot(z.rho, z.F_reference_strength_mvmr, color=col, lw=1.4, label=lab, marker='o', ms=2.5)
cu = CF[(CF.rho == 0)].F_custom_multistart.iloc[0]; ax.plot([0], [cu], marker='x', color=RED, ms=6, ls=''); ax.text(0.02, cu - 4, f'custom global Q minimum {cu:.1f}\n(identical for all exposures; not used)', fontsize=4.8, color=RED, va='top')
ax.axhline(10, color=RED, lw=0.7, ls='--'); ax.set_ylim(0, 110); ax.set_xlabel('assumed correlation of exposure estimation errors (ρ)'); ax.set_ylabel('conditional F (MVMR::strength_mvmr)')
ax.legend(fontsize=5.2, loc='center right'); ax.set_title('MVMR conditional F: reference implementation (318 variants)'); tag(ax, 'd')
# e
ax = fig.add_subplot(gs[2, 0]); L = LDc.copy(); L['lab'] = L.matrix.str.replace('FG_HYPERTENSION', 'FinnGen').str.replace('HYPERTENSION', 'Hyp.')
L = L[~L.matrix.str.contains(r'\((?:DBP|Hyp)', regex=True) | L.matrix.str.contains('GLS')]
for i, r in enumerate(L.itertuples()):
    ax.add_patch(Rectangle((0, i - 0.32), r.median_abs_diff, 0.3, color=BLUE, lw=0)); ax.add_patch(Rectangle((0, i + 0.02), r.p99_abs_diff, 0.3, color=ORANGE, alpha=0.8, lw=0))
    ax.text(r.p99_abs_diff + 0.003, i + 0.17, f'max {r.max_abs_diff:.2f}', fontsize=4.6, va='center', color=INK2)
ax.set_yticks(range(len(L))); ax.set_yticklabels(L.lab, fontsize=5.0); ax.set_ylim(len(L) - 0.4, -0.6); ax.set_xlim(0, 0.14)
ax.legend(handles=[Rectangle((0, 0), 1, 1, color=BLUE, label='median |r_hap − r_dosage|'), Rectangle((0, 0), 1, 1, color=ORANGE, label='99th percentile')], fontsize=4.8, loc='upper right')
ax.set_xlabel('absolute difference in signed LD r (same 503 individuals)'); clean(ax, left=False); ax.tick_params(axis='y', length=0)
ax.set_title('Haplotype vs genotype-dosage LD'); tag(ax, 'e')
# f
ax = fig.add_subplot(gs[2, 1]); f_ = FA[FA.used == True]
hb = ax.hexbin(f_.fin_1000G_af_kp_alt, f_.finngen_af_kp_alt, gridsize=16, extent=(0, 1, 0, 1), cmap=LinearSegmentedColormap.from_list('h', ['#f4f8fe'] + SEQ[1:]), mincnt=1, lw=0.2, edgecolors='white')
ax.plot([0, 1], [0, 1], color=INK2, lw=0.6, ls='--'); pal = f_[f_.palindromic == True]
amb = pal.palindrome_rule.str.startswith('ambiguous')
ax.scatter(pal.fin_1000G_af_kp_alt[~amb], pal.finngen_af_kp_alt[~amb], marker='x', s=12, color=AQUA, lw=0.8, label=f'palindromic, retained ({int((~amb).sum())})')
ax.scatter(pal.fin_1000G_af_kp_alt[amb], pal.finngen_af_kp_alt[amb], marker='x', s=14, color=RED, lw=0.9, label=f'palindromic, ambiguous ({int(amb.sum())}; excluded in sensitivity)')
ax.set_xlabel('1000 Genomes FIN frequency, effect allele'); ax.set_ylabel('FinnGen R12 frequency, effect allele'); ax.legend(fontsize=4.6, loc='upper left')
ax.text(0.98, 0.04, f'max |FinnGen − 1000G FIN| {f_.af_diff_vs_FIN.max():.2f}\n(vs 1000G EUR: {f_.af_diff.max():.2f})', transform=ax.transAxes, ha='right', fontsize=4.9, color=INK2)
cb = fig.colorbar(hb, ax=ax, fraction=0.04, pad=0.02); cb.outline.set_visible(False); cb.ax.tick_params(labelsize=5)
ax.set_title('FinnGen R12 harmonisation audit (116 instruments)'); tag(ax, 'f')
save(fig, 'FigS6_verification')
