"""Figure 2 (v3): Mendelian randomization with instruments clumped on 1000G EUR LD and estimates from MendelianRandomization / MRPRESSO."""
from figstyle import *
import numpy as np, pandas as pd
from matplotlib.colors import Normalize, LinearSegmentedColormap
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from scipy.stats import t as tdist, norm

R = pd.read_csv('step3v3_mr_results_R.csv'); PR = pd.read_csv('step3v3_forward_pair_ledger.csv'); P = pd.read_csv('step3v3_mr_results_python.csv')
V = pd.read_csv('validation/mr_python_vs_R.csv')
g = lambda code, an, m: R[(R.code == code) & (R.analysis == an) & (R.method == m)]
fig = plt.figure(figsize=(W, 9.0))

# ---- a  LD among candidate pairs: r2 vs distance (2D density) with decisions
ax = fig.add_axes([0.09, 0.705, 0.34, 0.235])
k = PR[PR.r2.notna()].copy(); k['lr2'] = np.log10(k.r2.clip(lower=1e-6))
hb = ax.hexbin(k.dist / 1e6, k.lr2, gridsize=(22, 14), extent=(0, 10, -6, 0), cmap=LinearSegmentedColormap.from_list('h', ['#f4f8fe'] + SEQ[1:]), mincnt=1, lw=0.2, edgecolors='white')
for thr, lab, col in [(np.log10(0.001), 'r² = 0.001 (clumping threshold)', RED), (np.log10(0.05), 'r² = 0.05 (Ensembl reporting floor, v2)', YELLOW)]:
    ax.axhline(thr, color=col, lw=1, ls='--'); ax.text(9.9, thr + 0.12, lab, ha='right', fontsize=5.4, color=col)
ax.axvline(1, color=INK2, lw=0.7, ls=':'); ax.text(1.08, -5.7, '1 Mb', fontsize=5.4, color=INK2)
ax.set_xlabel('distance between candidate LTL instruments (Mb)'); ax.set_ylabel('log$_{10}$ r² (1000G EUR)')
cb = fig.colorbar(hb, ax=ax, fraction=0.04, pad=0.02); cb.outline.set_visible(False); cb.ax.tick_params(labelsize=5); cb.set_label('pairs', fontsize=5.5)
nk, nu = len(k), int(PR.r2.isna().sum()); nhi = int((k.r2 >= 0.001).sum()); n05 = int((k.r2 >= 0.05).sum())
ax.text(0.99, 0.97, f'{len(PR)} pairs within 10 Mb: {nk} computed, {nu} unknown (chrX)\n{nhi} with r² ≥ 0.001; {n05} with r² ≥ 0.05', transform=ax.transAxes, ha='right', va='top', fontsize=5.6, color=INK)
ax.set_title('LD among candidate instruments'); tag(ax, 'a', x=-0.16)

# ---- b  SNP scatter with R-package fits (MVP EUR SBP)
s = pd.read_csv('mr_v3_inputs/snp_level_S_mvpEU.csv'); sg = np.sign(s.bx); s['x'] = s.bx * sg; s['y'] = s.by * sg
ax = fig.add_axes([0.6, 0.705, 0.31, 0.19]); axt = fig.add_axes([0.6, 0.9, 0.31, 0.045], sharex=ax); axr = fig.add_axes([0.915, 0.705, 0.06, 0.19], sharey=ax)
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

# ---- c  primary outcomes: analysis × outcome matrix
ax = fig.add_axes([0.3, 0.335, 0.55, 0.28])
rows = [('r2<0.001 (primary)', 'IVW (MRE, t)', 'IVW, r² < 0.001 (primary)'), ('r2<0.001 (primary)', 'MR-Egger (t)', 'MR-Egger'), ('r2<0.001 (primary)', 'Weighted median', 'Weighted median'),
        ('r2<0.001 (primary)', 'Weighted mode (MBE)', 'Weighted mode (MBE)'), ('r2<0.001 (primary)', 'PRESSO', 'MR-PRESSO (outlier-corrected)'),
        ('r2<0.001, Steiger-filtered', 'IVW (MRE, t)', 'IVW, Steiger-filtered'), ('r2<0.001, excluding TERC/TERT regions', 'IVW (MRE, t)', 'IVW, without TERC/TERT regions'),
        ('r2<0.001, excluding ambiguous palindromic', 'IVW (MRE, t)', 'IVW, without ambiguous palindromes'),
        ('r2<0.001, second-order weights, exposure SE from P', 'IVW (MRE, t, second-order weights)', 'IVW 2nd-order, exposure SE from P'),
        ('r2<0.001, second-order weights, exposure SE from N and allele frequency', 'IVW (MRE, t, second-order weights)', 'IVW 2nd-order, exposure SE from N, AF'),
        ('r2<0.1, covariance-aware', 'IVW, correlated instruments (GLS)', 'IVW GLS, r² < 0.1 (covariance-aware)'), ('r2<0.1, covariance-aware', 'MR-Egger, correlated instruments (GLS)', 'MR-Egger GLS, r² < 0.1'),
        ('r2<0.1, ignoring correlation (for contrast)', 'IVW (MRE, t)', 'IVW, r² < 0.1, correlation ignored')]
cols = [('S_mvpEU', 'SBP'), ('D_mvpEU', 'DBP'), ('H_mvpEU', 'Hypertension')]
Z = np.full((len(rows), 3), np.nan); T = [[''] * 3 for _ in rows]
for i, (an, m, _) in enumerate(rows):
    for j, (c, _) in enumerate(cols):
        x = R[(R.code == c) & (R.analysis == an) & (R.method.str.startswith('MR-PRESSO') if m == 'PRESSO' else (R.method == m))]
        if x.empty: continue
        x = x.iloc[0]; Z[i, j] = x.b / x.se
        est = f'OR {np.exp(x.b):.2f}' if c.startswith('H_') else f'{x.b:.3f}'
        T[i][j] = f'{est} ({int(x.nsnp)})' + ('**' if x.p < 0.05 / 6 else ('*' if x.p < 0.05 else ''))
im = ax.imshow(Z, cmap=DIV, vmin=-8, vmax=8, aspect='auto')
for i in range(len(rows)):
    for j in range(3): ax.text(j, i, T[i][j], ha='center', va='center', fontsize=5.4, color='white' if abs(Z[i, j]) > 5.5 else INK)
ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[2] for r in rows], fontsize=5.6); ax.set_xticks(range(3)); ax.set_xticklabels([c[1] + ' (MVP EUR)' for c in cols], fontsize=6); ax.xaxis.tick_top(); ax.tick_params(length=0)
for yy in [4.5, 7.5, 9.5]: ax.axhline(yy, color='white', lw=2)
for sp in ax.spines.values(): sp.set_visible(False)
for j, (c, _) in enumerate(cols):
    iv_ = g(c, 'r2<0.001 (primary)', 'IVW (MRE, t)').iloc[0]; e_ = g(c, 'r2<0.001 (primary)', 'MR-Egger (t)').iloc[0]
    ax.text(j, len(rows) - 0.25, f'I² {iv_.I2:.0f}%;  Egger int. P {e_.intercept_p:.3f}', ha='center', va='top', fontsize=5.2, color=RED if e_.intercept_p < 0.05 else INK2)
ax.set_ylim(len(rows) + 0.5, -0.5)
ax.set_title('Primary outcomes: estimator and sensitivity analyses (per SD longer LTL; n instruments)', pad=16); tag(ax, 'c', x=-0.52, y=1.07)
cax = ax.inset_axes([1.04, 0.25, 0.025, 0.5]); cb = fig.colorbar(im, cax=cax); cb.outline.set_visible(False); cb.ax.tick_params(labelsize=5); cb.set_label('z', fontsize=5.5)

# ---- e  reverse MR: estimate CI against detectable range
ax = fig.add_axes([0.2, 0.075, 0.5, 0.2])
rev = [('SBP->L_codd', 'SBP → LTL'), ('DBP->L_codd', 'DBP → LTL'), ('Hypertension->L_codd', 'HTN liability → LTL'), ('SBP->L_nkSA', 'SBP → LTL (SAS)'), ('DBP->L_nkSA', 'DBP → LTL (SAS)'), ('Hypertension->L_nkSA', 'HTN → LTL (SAS)')]
for i, (c, lab) in enumerate(rev):
    x = g(c, 'r2<0.001 (primary)', 'IVW (MRE, t)').iloc[0]; k_ = int(x.nsnp)
    m05 = (tdist.ppf(0.975, k_ - 1) + norm.ppf(0.8)) * x.se; mb = (tdist.ppf(1 - 0.05 / 12, k_ - 1) + norm.ppf(0.8)) * x.se
    y = -i
    ax.add_patch(Rectangle((-mb, y - 0.36), 2 * mb, 0.72, color=MID, lw=0)); ax.add_patch(Rectangle((-m05, y - 0.36), 2 * m05, 0.72, color='#e3e0d6', lw=0))
    ax.add_patch(Rectangle((x.lo, y - 0.14), x.hi - x.lo, 0.28, color=BLUE if 'codd' in c else ORANGE, lw=0, alpha=0.85))
    ax.plot([x.b, x.b], [y - 0.26, y + 0.26], color=INK, lw=1.4)
    ax.text(0.21, y, f'{x.b:+.3f}   P = {x.p:.2g}   k = {k_}', va='center', fontsize=5.4, color=INK, transform=ax.transData)
ax.axvline(0, color=INK2, lw=0.6); ax.set_yticks([-i for i in range(len(rev))]); ax.set_yticklabels([r_[1] for r_ in rev], fontsize=5.8)
ax.set_xlim(-0.2, 0.2); ax.set_xlabel('IVW estimate, SD LTL per unit exposure'); ax.tick_params(axis='y', length=0); clean(ax, left=False)
ax.legend(handles=[Rectangle((0, 0), 1, 1, color='#e3e0d6', label='MDE, α = 0.05'), Rectangle((0, 0), 1, 1, color=MID, label='MDE, α = 0.05/6'), Rectangle((0, 0), 1, 1, color=BLUE, label='95% CI (Codd LTL, primary)'), Rectangle((0, 0), 1, 1, color=ORANGE, label='95% CI (South Asian LTL)')],
          fontsize=5, loc='upper left', ncol=1, bbox_to_anchor=(1.36, 1.02))
ax.set_title('Reverse MR: estimates against minimum detectable effects (80% power)'); tag(ax, 'd', x=-0.3)

fig.text(0.01, 0.0, 'Estimates from MendelianRandomization 0.10.0 (IVW multiplicative random effects, t(k−1); MR-Egger t(k−2); weighted median; MBE mode) and MRPRESSO 1.0 (2,400 simulations; primary forward outcomes only). '
         'Instruments: PLINK-style clumping at r² < 0.001 within 10 Mb using 1000 Genomes phase 3 EUR haplotypes. ** P < 0.0083 (primary family), * P < 0.05.', fontsize=5.4, color=MUTED, wrap=True)
save(fig, 'Fig2_mendelian_randomization')
