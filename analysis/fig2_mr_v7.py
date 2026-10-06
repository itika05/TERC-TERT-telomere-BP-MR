"""Figure 2 (v6): single 95% CI bands (50% bands removed), taller panel c, MR-PRESSO at B = 10^6, reverse MR on source LTL statistics.
Figure 2 (v5): (c) endpoint-specific interval bands replace the estimator heatmap (full grid moved to Supplementary Figure S3); (d) MDE relabelled.
Figure 2 (v4): Mendelian randomization with instruments clumped on 1000G EUR LD and estimates from MendelianRandomization / MRPRESSO."""
from figstyle import *
import numpy as np, pandas as pd
from matplotlib.colors import Normalize, LinearSegmentedColormap
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from scipy.stats import t as tdist, norm

R = pd.read_csv('step8_mr_results_v6.csv'); PR = pd.read_csv('step6_v4_forward_pair_ledger.csv'); MV = pd.read_csv('step8_mvmr_v6_estimates.csv'); MV = MV[MV.set.str.startswith('primary')]
PS = pd.read_csv('step7_presso_seed_runs.csv'); FV = pd.read_csv('step6_v4_forward_variant_ledger.csv')
g = lambda code, an, m: R[(R.code == code) & (R.analysis == an) & (R.method == m)]
fig = plt.figure(figsize=(W, 9.3))

# ---- a  LD among candidate pairs: r2 vs distance (2D density) with decisions
ax = fig.add_axes([0.09, 0.76, 0.34, 0.2])
k = PR[PR.r2.notna()].copy(); k['lr2'] = np.log10(k.r2.clip(lower=1e-6))
hb = ax.hexbin(k.dist / 1e6, k.lr2, gridsize=(22, 14), extent=(0, 10, -6, 0), cmap=LinearSegmentedColormap.from_list('h', ['#f4f8fe'] + SEQ[1:]), mincnt=1, lw=0.2, edgecolors='white')
for thr, lab, col in [(np.log10(0.001), 'r² = 0.001 (clumping threshold)', RED)]:
    ax.axhline(thr, color=col, lw=1, ls='--'); ax.text(9.9, thr + 0.12, lab, ha='right', fontsize=5.4, color=col)
ax.axvline(1, color=INK2, lw=0.7, ls=':'); ax.text(1.08, -0.3, '1 Mb', fontsize=6.1, color=INK2)
ax.set_ylim(-6.3, 1.7); ax.set_yticks([-6, -5, -4, -3, -2, -1, 0]); ax.set_xlabel('distance between candidate LTL instruments (Mb)'); ax.set_ylabel('log$_{10}$ r² (1000G EUR)')
cb = fig.colorbar(hb, ax=ax, fraction=0.04, pad=0.02); cb.outline.set_visible(False); cb.ax.tick_params(labelsize=5); cb.set_label('pairs', fontsize=5.5)
nk, nu = len(k), int(PR.r2.isna().sum()); nhi = int((k.r2 >= 0.001).sum()); n05 = int((k.r2 >= 0.05).sum())
nX = int(PR.ld_status.str.contains('chrX').sum()); nM = int(PR.ld_status.str.contains('MHC').sum())
ax.text(0.99, 0.03, f'{len(PR)} pairs within 10 Mb among {len(FV)} candidates;\n{nk} shown ({nhi} with r² ≥ 0.001); {nu} involving an\nexcluded candidate (chrX {nX}, MHC {nM}) not shown', transform=ax.transAxes, ha='right', va='bottom', fontsize=6.1, color=INK, linespacing=1.15, bbox=dict(fc='white', ec=GRID, alpha=1.0, pad=1.5), zorder=6)
ax.set_title('LD among candidate instruments (pre-exclusion pair set)'); tag(ax, 'a', x=-0.16)

# ---- b  SNP scatter with R-package fits (MVP EUR SBP)
s = pd.read_csv('mr_v4_inputs/snp_level_S_mvpEU.csv'); sg = np.sign(s.bx); s['x'] = s.bx * sg; s['y'] = s.by * sg
ax = fig.add_axes([0.6, 0.76, 0.31, 0.16]); axt = fig.add_axes([0.6, 0.922, 0.31, 0.04], sharex=ax); axr = fig.add_axes([0.915, 0.76, 0.06, 0.16], sharey=ax)
ax.errorbar(s.x, s.y, xerr=1.96 * s.sx, yerr=1.96 * s.sy, fmt='none', ecolor='#e6e3db', elinewidth=0.5, zorder=1)
sc = ax.scatter(s.x, s.y, c=np.clip(s.steiger_z, 0, 60), cmap=LinearSegmentedColormap.from_list('z', SEQ[1:]), vmin=0, vmax=60, s=11, edgecolor='white', lw=0.3, zorder=3)
xx = np.linspace(0, s.x.max() * 1.05, 50); iv = g('S_mvpEU', 'r2<0.001 (primary)', 'IVW (MRE, t)').iloc[0]
ax.fill_between(xx, iv.lo * xx, iv.hi * xx, color=BLUE, alpha=0.13, lw=0)
for m, col, ls, lab in [('IVW (MRE, t)', BLUE, '-', 'IVW'), ('Weighted median', AQUA, '--', 'W-median'), ('Weighted mode (MBE)', YELLOW, '-.', 'W-mode')]:
    r = g('S_mvpEU', 'r2<0.001 (primary)', m).iloc[0]; ax.plot(xx, r.b * xx, color=col, ls=ls, lw=1.1, label=lab)
e = g('S_mvpEU', 'r2<0.001 (primary)', 'MR-Egger (t)').iloc[0]; ax.plot(xx, e.intercept + e.b * xx, color=ORANGE, lw=1.1, label='MR-Egger')
pr = R[(R.code == 'S_mvpEU') & R.method.str.startswith('MR-PRESSO')]
if len(pr): ax.plot(xx, pr.b.iloc[0] * xx, color=VIOLET, ls=':', lw=1.1, label='MR-PRESSO')
ax.axhline(0, color=AXIS, lw=0.6); ax.set_xlabel('SNP → LTL (Codd 2021), oriented positive'); ax.set_ylabel('SNP → SBP (MVP EUR)')
ax.legend(*ax.get_legend_handles_labels(), loc='upper left', fontsize=5.0, ncol=2, handlelength=1.5, columnspacing=0.7, frameon=True, framealpha=0.9, edgecolor='none')
axt.hist(s.x, bins=30, color=SEQ[2], lw=0); axt.axis('off'); axr.hist(s.y, bins=30, orientation='horizontal', color=SEQ[2], lw=0); axr.axis('off')
cax = ax.inset_axes([0.6, 0.1, 0.36, 0.03]); cb = fig.colorbar(sc, cax=cax, orientation='horizontal'); cb.outline.set_visible(False); cb.ax.tick_params(labelsize=5, pad=1); cax.set_title('Steiger z', fontsize=5, color=INK2, fontweight='normal', pad=1, loc='left')
axt.text(-0.25, 1.25, 'b', transform=axt.transAxes, fontsize=10, fontweight='bold'); axt.text(0.0, 1.25, f'SBP, {int(iv.nsnp)} instruments (r² < 0.001)', transform=axt.transAxes, fontsize=7.5, fontweight='bold')

# ---- c  endpoint-specific interval bands (primary IVW highlighted; everything else secondary)
rowsc = [('IVW, r² < 0.001', ('r2<0.001 (primary)', 'IVW (MRE, t)'), 'primary'),
         ('MR-Egger', ('r2<0.001 (primary)', 'MR-Egger (t)'), 'est'), ('Weighted median', ('r2<0.001 (primary)', 'Weighted median'), 'est'),
         ('Weighted mode (MBE)', ('r2<0.001 (primary)', 'Weighted mode (MBE)'), 'est'), ('MR-PRESSO (B = 10⁶)', 'PRESSO', 'est'),
         ('Steiger-filtered', ('r2<0.001, Steiger-filtered', 'IVW (MRE, t)'), 'set'), ('PLINK 1.9 clumping', 'PLINK', 'set'),
         ('Published Codd instrument', ('Codd 2021 published MR instrument', 'IVW (MRE, t)'), 'set'),
         ('Without TERC/TERT regions', ('r2<0.001, excluding TERC/TERT regions', 'IVW (MRE, t)'), 'set'),
         ('GLS, r² < 0.1 (correlated)', ('r2<0.1, covariance-aware', 'IVW, correlated instruments (GLS)'), 'set'),
         ('Phenome-screened set', ('r2<0.001, excluding instruments with confounder-trait associations', 'IVW (MRE, t)'), 'pl'),
         ('MVMR | BMI, lymphocytes', 'MVMR', 'pl')]
endp = [('S_mvpEU', 'SBP, MVP', 'SD per SD LTL', False), ('D_mvpEU', 'DBP, MVP', 'SD per SD LTL', False),
        ('H_mvpEU', 'Hypertension, MVP', 'OR per SD LTL', True), ('H_fg12', 'Hypertension, FinnGen R12 (secondary)', 'OR per SD LTL', True)]
GC = {'primary': BLUE, 'est': VIOLET, 'set': AQUA, 'pl': ORANGE}
def getrow(code, spec):
    if spec == 'MVMR':
        x = MV[(MV.outcome == code) & (MV.exposure == 'LTL') & (MV.method == 'MV-IVW (random, t)')]
        return None if x.empty else x.iloc[0]
    if spec == 'PLINK':
        x = R[(R.code == code) & (R.analysis == 'PLINK 1.9 clumping') & (R.method == 'IVW (MRE, t)')]; return None if x.empty else x.iloc[0]
    if spec == 'PRESSO':
        x = R[(R.code == code) & (R.method == 'MR-PRESSO (outlier-corrected)')]; return None if x.empty else x.iloc[0]
    x = R[(R.code == code) & (R.analysis == spec[0]) & (R.method == spec[1])]; return None if x.empty else x.iloc[0]
x0, x1 = 0.215, 0.925; gap = 0.032; wpan = (x1 - x0 - 3 * gap) / 4
for j, (code, lab, unit, orr) in enumerate(endp):
    ax = fig.add_axes([x0 + j * (wpan + gap), 0.33, wpan, 0.30])
    tr = (lambda v: np.exp(v)) if orr else (lambda v: v)
    prim = getrow(code, rowsc[0][1])
    ax.axvspan(tr(prim.lo), tr(prim.hi), color=BLUE, alpha=0.07, lw=0)
    ax.axvline(tr(0), color=INK2, lw=0.6)
    for i, (rl, spec, grp) in enumerate(rowsc):
        r = getrow(code, spec); y = i
        if r is None:
            ax.text(0.62, y, 'not estimated', transform=ax.get_yaxis_transform(), fontsize=5.6, color=MUTED, va='center', ha='left', zorder=6, bbox=dict(fc='white', ec='none', pad=0.6)); continue
        col = GC[grp]; lo, hi, b = tr(r.lo), tr(r.hi), tr(r.b)
        ax.add_patch(Rectangle((lo, y - 0.3), hi - lo, 0.6, color=col, alpha=0.45 if grp != 'primary' else 0.7, lw=0))
        ax.plot([b, b], [y - 0.38, y + 0.38], color=INK, lw=1.1, solid_capstyle='butt')
        if spec == 'PRESSO':
            ps = PS[PS.code == code]
            if len(ps):
                for bb in ps.b: ax.plot([tr(bb)] * 2, [y - 0.42, y - 0.3], color=INK, lw=0.6)
        k_ = int(r.nsnp)
        star = '**' if (grp == 'primary' and code != 'H_fg12' and r.p < 0.05 / 6) else ''
        ax.text(1.02, y, f'{k_}', transform=ax.get_yaxis_transform(), fontsize=6.6, color=INK, va='center', ha='left', clip_on=False)   # v9: k printed alone (the ** marker had made 116 read as 119)
        if star: ax.annotate('**', (hi, y), xytext=(2.5, 0), textcoords='offset points', fontsize=7.0, color=INK, va='center', ha='left', fontweight='bold', zorder=7)   # marker beside the interval, clear of the k column
    if orr:
        ax.set_xscale('log'); ax.set_xlim(0.8, 1.75); ax.set_xticks([0.8, 1.0, 1.2, 1.4, 1.6]); ax.set_xticklabels(['0.8', '1.0', '1.2', '1.4', '1.6'])
    else:
        hmax = max(getrow(code, sp).hi for _, sp, _ in rowsc if getrow(code, sp) is not None); ax.set_xlim(-0.08, max(0.36, hmax + 0.015)); ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4] if hmax > 0.39 else [0, 0.1, 0.2, 0.3])
    ax.minorticks_off()
    ax.set_ylim(len(rowsc) - 0.4, -0.7); ax.set_yticks(range(len(rowsc)))
    ax.set_yticklabels([r_[0] for r_ in rowsc] if j == 0 else [], fontsize=7.2); ax.tick_params(axis='y', length=0)
    for yy in [0.5, 4.5, 9.5]: ax.axhline(yy, color=GRID, lw=0.6)
    ax.set_title(lab, fontsize=7.4, pad=3, color=INK); ax.set_xlabel(unit, fontsize=6.6); ax.tick_params(axis='x', labelsize=6.4); clean(ax, left=False)
    iv_ = getrow(code, rowsc[0][1]); e_ = getrow(code, rowsc[1][1])
    ax.text(0.02, -0.175, f'I² {iv_.I2:.0f}% · Egger int. P {e_.intercept_p:.3f}', transform=ax.transAxes, fontsize=5.9, color=RED if e_.intercept_p < 0.05 else INK2)
    if j == 0: tag(ax, 'c', x=-0.95, y=1.07)
fig.text(x0, 0.69, 'Primary outcomes: 95% CI band and estimate (notch); k at right (** primary test, P < 0.0083)', fontsize=7.0, fontweight='bold', color=INK)
lg = [Rectangle((0, 0), 1, 1, color=GC[g_], alpha=0.6, label=l_) for g_, l_ in [('primary', 'primary IVW'), ('est', 'other estimators'), ('set', 'other instrument sets'), ('pl', 'pleiotropy sensitivity')]]
fig.legend(handles=lg + [Line2D([0], [0], color=INK, lw=0.6, label='MRPRESSO 1.0, B = 2,400, 3 seeds (ticks)')], loc='lower left', bbox_to_anchor=(x0, 0.652), ncol=5, fontsize=6.0, handlelength=1.4, columnspacing=1.0)

# ---- e  reverse MR: estimate CI against detectable range
ax = fig.add_axes([0.2, 0.075, 0.5, 0.15])
rev = [('SBP->L_codd', 'SBP → LTL'), ('DBP->L_codd', 'DBP → LTL'), ('Hypertension->L_codd', 'Hypertension liability → LTL'), ('SBP->L_nkSA', 'SBP → LTL (SAS)'), ('DBP->L_nkSA', 'DBP → LTL (SAS)'), ('Hypertension->L_nkSA', 'Hypertension → LTL (SAS)')]
for i, (c, lab) in enumerate(rev):
    x = g(c, 'r2<0.001 (primary)', 'IVW (MRE, t)').iloc[0]; k_ = int(x.nsnp)
    m05 = (tdist.ppf(0.975, k_ - 1) + norm.ppf(0.8)) * x.se; mb = (tdist.ppf(1 - 0.05 / 12, k_ - 1) + norm.ppf(0.8)) * x.se
    y = -i - (0.6 if i >= 3 else 0)
    ax.add_patch(Rectangle((-mb, y - 0.36), 2 * mb, 0.72, color=MID, lw=0)); ax.add_patch(Rectangle((-m05, y - 0.36), 2 * m05, 0.72, color='#e3e0d6', lw=0))
    ax.add_patch(Rectangle((x.lo, y - 0.14), x.hi - x.lo, 0.28, color=BLUE if 'codd' in c else ORANGE, lw=0, alpha=0.85))
    ax.plot([x.b, x.b], [y - 0.26, y + 0.26], color=INK, lw=1.4)
    ax.text(0.21, y, f'{x.b:+.3f}   P = {x.p:.2g}   k = {k_}'.replace('-', '−'), va='center', fontsize=6.2, color=INK, transform=ax.transData)
ax.axvline(0, color=INK2, lw=0.6); ax.axhline(-2.8, color=AXIS, lw=0.6, ls='--'); ax.text(-0.195, -2.95, 'exploratory: South Asian LTL outcome (Nakao 2026)', fontsize=5.8, color=INK2, va='top', zorder=6, bbox=dict(fc='white', ec='none', pad=0.8))
ax.set_yticks([-i - (0.6 if i >= 3 else 0) for i in range(len(rev))]); ax.set_yticklabels([r_[1] for r_ in rev], fontsize=6.6)
ax.set_xlim(-0.2, 0.2); ax.set_xlabel('IVW estimate: SD LTL per SD BP (SBP, DBP) or per unit log odds of hypertension liability'); ax.tick_params(axis='y', length=0); clean(ax, left=False)
ax.legend(handles=[Rectangle((0, 0), 1, 1, color='#e3e0d6', label='80%-power threshold, α = 0.05'), Rectangle((0, 0), 1, 1, color=MID, label='80%-power threshold, α = 0.05/6'), Rectangle((0, 0), 1, 1, color=BLUE, label='95% CI (Codd LTL source statistics, primary)'), Rectangle((0, 0), 1, 1, color=ORANGE, label='95% CI (South Asian LTL)')],
          fontsize=6.1, loc='upper center', ncol=2, bbox_to_anchor=(0.62, -0.3), frameon=False)
ax.set_title('Reverse MR: 95% CIs with 80%-power sensitivity thresholds (not equivalence bounds)'); tag(ax, 'd', x=-0.3)

# v8: figure-level footnote moved to the figure legend in the manuscript (print-width proof, V7-10)
save(fig, 'Fig2_mendelian_randomization')
