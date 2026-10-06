"""Supplementary Figure S6 (v5): verification analyses requested in the v4 revision guide.
(a) MR-PRESSO Monte Carlo stability across three seeds; (b) custom versus PLINK 1.9 clumping; (c) P-derived versus reported SE;
(d) MVMR conditional F under assumed exposure-error correlation; (e) TERT SuSiE: L = 10 versus L = 20; (f) FinnGen R12 allele-frequency concordance."""
from figstyle import *
import json, numpy as np, pandas as pd
from scipy.stats import norm
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
PS = pd.read_csv('step7_presso_seed_runs.csv'); PC = pd.read_csv('step7_presso_consensus_ivw.csv'); PL = pd.read_csv('step7_plink_clump_comparison.csv'); PI = pd.read_csv('step7_plink_clump_ivw.csv')
CF = pd.read_csv('step7_mvmr_condF_covariance_grid.csv'); D10 = pd.read_csv('step6_coloc_v4_susie_diagnostics.csv'); D20 = pd.read_csv('step7_susie_L20_TERT_diagnostics.csv')
S10 = pd.read_csv('step6_coloc_v4_susie.csv'); S20 = pd.read_csv('step7_susie_L20_TERT_coloc.csv'); FA = pd.read_csv('step7_finngen_harmonisation_audit.csv')
cal = json.load(open('source_data/kp_v4_round2.json'))['cal']; SPV = pd.read_csv('step6_codd_source_vs_portal.csv')
fig = plt.figure(figsize=(W, 9.6)); gs = fig.add_gridspec(3, 2, hspace=0.62, wspace=0.42, height_ratios=[1, 1, 1])
LAB = {'S_mvpEU': 'SBP', 'D_mvpEU': 'DBP', 'H_mvpEU': 'Hypertension (log OR)'}; COL = {'S_mvpEU': BLUE, 'D_mvpEU': AQUA, 'H_mvpEU': ORANGE}
# a PRESSO
ax = fig.add_subplot(gs[0, 0])
for i, code in enumerate(['S_mvpEU', 'D_mvpEU', 'H_mvpEU']):
    g = PS[PS.code == code].sort_values('seed')
    for k, r in enumerate(g.itertuples()):
        y = i + (k - 1) * 0.22
        ax.add_patch(Rectangle((r.b - 1.96 * r.se, y - 0.08), 3.92 * r.se, 0.16, color=COL[code], alpha=0.45, lw=0)); ax.plot([r.b] * 2, [y - 0.1, y + 0.1], color=INK, lw=1)
        ax.text(0.33, y, f'seed {r.seed}: {r.n_outliers} outliers, Jaccard vs v4 {r.jaccard_vs_v4:.2f}', fontsize=4.6, va='center', color=INK2)
    for lab_, ls in [('excl. outliers flagged in all 3 seeds', '-'), ('excl. outliers flagged in any seed', ':')]:
        c = PC[(PC.code == code) & (PC.set == lab_)].iloc[0]; ax.plot([c.b] * 2, [i - 0.42, i + 0.42], color=COL[code], lw=1.2, ls=ls)
ax.set_yticks(range(3)); ax.set_yticklabels([LAB[c] for c in ['S_mvpEU', 'D_mvpEU', 'H_mvpEU']], fontsize=5.8); ax.set_ylim(2.6, -0.6); ax.set_xlim(0, 0.62)
ax.axvline(0, color=INK2, lw=0.5); ax.set_xlabel('outlier-corrected IVW per SD LTL'); clean(ax, left=False); ax.tick_params(axis='y', length=0)
ax.legend(handles=[Line2D([], [], color=INK2, lw=1.2, label='fixed-set IVW, consensus outliers removed'), Line2D([], [], color=INK2, lw=1.2, ls=':', label='fixed-set IVW, union removed')], fontsize=4.8, loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2)
ax.set_title('MR-PRESSO stability (2,400 simulations × 3 seeds)'); tag(ax, 'a')
# b PLINK
ax = fig.add_subplot(gs[0, 1]); r = PL[PL.r2_threshold == 0.001].iloc[0]
both, co, po = int(r.both), int(r.custom_kept - r.both), int(r.plink_kept - r.both)
ax.barh([2], [co], left=[0], color=BLUE, alpha=0.5); ax.barh([2], [both], left=[co], color=VIOLET, alpha=0.7); ax.barh([2], [po], left=[co + both], color=ORANGE, alpha=0.5)
ax.text(co / 2, 2, f'{co}\ncustom only', ha='center', va='center', fontsize=4.8); ax.text(co + both / 2, 2, f'{both} in both', ha='center', va='center', fontsize=5.5, color='white', fontweight='bold'); ax.text(co + both + po / 2, 2, f'{po}\nPLINK only', ha='center', va='center', fontsize=4.8)
for i, code in enumerate(['S_mvpEU', 'D_mvpEU', 'H_mvpEU', 'H_fg12']):
    for k, ins in enumerate(['custom clumping (primary)', 'PLINK 1.9 --clump']):
        x = PI[(PI.code == code) & (PI.instrument == ins)].iloc[0]; y = -0.35 - i * 0.55 - k * 0.2; sc = 400
        ax.add_patch(Rectangle((60 + x.lo * sc, y - 0.07), (x.hi - x.lo) * sc, 0.14, color=BLUE if k == 0 else ORANGE, alpha=0.55, lw=0)); ax.plot([60 + x.b * sc] * 2, [y - 0.09, y + 0.09], color=INK, lw=1)
    ax.text(58, -0.45 - i * 0.55, {'H_fg12': 'Hyp. FinnGen'}.get(code, LAB.get(code, code).replace(' (log OR)', ' MVP')), ha='right', fontsize=5, va='center')
ax.plot([60, 60], [-2.55, 0.0], color=INK2, lw=0.5); ax.text(60, 0.95, 'IVW per SD LTL (0 at line; bar = 95% CI; blue custom, orange PLINK)', fontsize=4.8, color=INK2, ha='left')
for v in [0.1, 0.2]: ax.plot([60 + v * sc] * 2, [-2.55, 0.0], color=GRID, lw=0.5); ax.text(60 + v * sc, -2.65, f'{v}', fontsize=4.8, ha='center', color=INK2)
ax.set_xlim(0, 180); ax.set_ylim(-2.8, 2.6); ax.axis('off')
ax.set_title(f'Custom vs PLINK 1.9 clumping (r² < 0.001, 10 Mb; {int(r.n_candidates)} candidates)'); tag(ax, 'b', x=-0.08)
# c SE validation
ax = fig.add_subplot(gs[1, 0]); names = []
for i, (k, recs) in enumerate(cal.items()):
    X = pd.DataFrame(recs).iloc[:, :4]; X.columns = ['v', 'b', 'se', 'p']; X = X[(X.p > 1e-300) & (X.b != 0)]
    rt = 100 * (X.b.abs() / norm.isf(X.p / 2) / X.se - 1); q = np.percentile(rt, [0, 2.5, 25, 50, 75, 97.5, 100])
    ax.add_patch(Rectangle((q[1], i - 0.3), q[5] - q[1], 0.6, color=SEQ[1], lw=0)); ax.add_patch(Rectangle((q[2], i - 0.3), q[4] - q[2], 0.6, color=SEQ[4], lw=0))
    ax.plot([q[0], q[6]], [i, i], color=INK2, lw=0.5); ax.plot([q[3]] * 2, [i - 0.34, i + 0.34], color=INK, lw=1)
    names.append(k.replace('_MVPTraits_EU', ' MVP').replace('_bp_eu', '').replace('_BP_EU', '').replace('_BPtraits_EU', '').replace('GWAS_UKBiobank409k_eu', 'UKB 409k').replace('|', ' '))
i += 1; rt = 100 * (SPV.se_portal_from_P / SPV.se_source - 1).dropna(); q = np.percentile(rt, [0, 2.5, 25, 50, 75, 97.5, 100])
ax.add_patch(Rectangle((q[1], i - 0.3), q[5] - q[1], 0.6, color='#cdc4e8', lw=0)); ax.add_patch(Rectangle((q[2], i - 0.3), q[4] - q[2], 0.6, color=VIOLET, lw=0)); ax.plot([q[3]] * 2, [i - 0.34, i + 0.34], color=INK, lw=1)
names.append('Codd LTL vs source')
ax.set_yticks(range(len(names))); ax.set_yticklabels(names, fontsize=5.2); ax.set_ylim(len(names) - 0.4, -0.6); ax.axvline(0, color=INK2, lw=0.5)
ax.set_xlabel('SE from P vs reported SE (% difference)'); ax.set_xlim(-4.5, 1.5); clean(ax, left=False); ax.tick_params(axis='y', length=0)
ax.set_title('SE from P vs reported SE (top 1,500 records/dataset)'); tag(ax, 'c')
# d condF
ax = fig.add_subplot(gs[1, 1])
for e, col, lw, ls in [('LYM', AQUA, 1.3, '-'), ('LTL', VIOLET, 2.6, '-'), ('BMI', YELLOW, 1.1, '--')]:
    z = CF[CF.exposure == e].sort_values('rho'); ax.plot(z.rho, z.conditional_F, color=col, lw=lw, ls=ls, label={'LYM': 'lymphocyte count'}.get(e, e), marker='o', ms=2.5)
ax.fill_between(z.rho, 10, CF[CF.exposure == 'LTL'].sort_values('rho').conditional_F, color=VIOLET, alpha=0.08, lw=0)
ax.text(0.98, 0.9, 'Curves coincide where the minimising δ has all non-zero\ncomponents (the minimum of the Sanderson Q is invariant to rescaling)', transform=ax.transAxes, ha='right', fontsize=4.8, color=INK2)
ax.axhline(10, color=RED, lw=0.8, ls='--'); ax.text(0.3, 11, 'F = 10', ha='right', fontsize=5, color=RED)
ax.set_xlabel('assumed correlation of exposure estimation errors (ρ)'); ax.set_ylabel('conditional F'); ax.set_ylim(0, 55); ax.legend(fontsize=5.3, loc='lower left')
ax.set_title('MVMR conditional F under sample overlap'); tag(ax, 'd')
# e SuSiE L
ax = fig.add_subplot(gs[2, 0]); TRL = {'SBP': 'SBP', 'DBP': 'DBP', 'HYPERTENSION': 'Hyp. MVP', 'FG_HYPERTENSION': 'Hyp. FinnGen'}
rows = [(w, t) for w in [100, 250] for t in TRL]
for i, (w, t) in enumerate(rows):
    a10 = D10[(D10.locus == 'TERT') & (D10.window_kb == w) & (D10.pair == t) & (D10.trait == 'LTL')].n_credible_sets.iloc[0]
    a20 = D20[(D20.window_kb == w) & (D20.pair == t) & (D20.trait == 'LTL')].n_credible_sets.iloc[0]
    ax.add_patch(Rectangle((0, i - 0.32), a10, 0.3, color=BLUE, alpha=0.7, lw=0)); ax.add_patch(Rectangle((0, i + 0.02), a20, 0.3, color=RED, alpha=0.7, lw=0))
    ax.text(a10 + 0.3, i - 0.17, f'{a10}', fontsize=4.8, va='center'); ax.text(a20 + 0.3, i + 0.17, f'{a20}', fontsize=4.8, va='center')
    lead = '5:1285974:C:A'; x10 = S10[(S10.locus == 'TERT') & (S10.trait == t) & (S10.window_kb == w) & np.isclose(S10.p12, 1e-5)]
    m10 = x10.PP4.max() if len(x10) and x10.PP4.notna().any() else np.nan; m20 = S20[(S20.trait == t) & (S20.window_kb == w)].maxPP4.iloc[0]
    f = lambda v: 'no CS' if np.isnan(v) else f'{v:.2f}'
    ax.text(24, i, f'max PP4: {f(m10)} → {f(m20)}', fontsize=4.9, va='center', color=INK)
ax.set_yticks(range(len(rows))); ax.set_yticklabels([f'±{w} kb, pair {TRL[t]}' for w, t in rows], fontsize=5.2); ax.set_ylim(len(rows) - 0.5, -0.6); ax.set_xlim(0, 34)
ax.axvline(10, color=BLUE, lw=0.6, ls=':'); ax.axvline(20, color=RED, lw=0.6, ls=':'); ax.set_xlabel('LTL credible sets at TERT (blue L = 10, red L = 20)'); clean(ax, left=False); ax.tick_params(axis='y', length=0)
ax.set_title('TERT SuSiE: L = 10 vs L = 20'); tag(ax, 'e')
# f FinnGen AF
ax = fig.add_subplot(gs[2, 1]); f_ = FA[FA.used == True]
hb = ax.hexbin(f_.eur_af_kp_alt, f_.finngen_af_kp_alt, gridsize=16, extent=(0, 1, 0, 1), cmap=LinearSegmentedColormap.from_list('h', ['#f4f8fe'] + SEQ[1:]), mincnt=1, lw=0.2, edgecolors='white')
ax.plot([0, 1], [0, 1], color=INK2, lw=0.6, ls='--'); pal = f_[f_.palindromic == True]
ax.scatter(pal.eur_af_kp_alt, pal.finngen_af_kp_alt, marker='x', s=12, color=ORANGE, lw=0.8, label=f'palindromic ({len(pal)}; {int((pal.palindromic_af_check == "consistent").sum())} consistent, {int(pal.palindromic_af_check.str.startswith("MAF").sum())} MAF > 0.42)')
ax.set_xlabel('1000G EUR frequency, effect allele'); ax.set_ylabel('FinnGen R12 frequency, effect allele'); ax.legend(fontsize=4.8, loc='upper left')
ax.text(0.98, 0.04, f'{len(f_)} of {len(FA)} instruments matched on identical alleles;\nno strand flips or mismatches; max |ΔAF| {f_.af_diff.max():.2f}', transform=ax.transAxes, ha='right', fontsize=4.9, color=INK2)
cb = fig.colorbar(hb, ax=ax, fraction=0.04, pad=0.02); cb.outline.set_visible(False); cb.ax.tick_params(labelsize=5)
ax.set_title('FinnGen R12 harmonisation audit (116 instruments)'); tag(ax, 'f')
save(fig, 'FigS6_verification')
