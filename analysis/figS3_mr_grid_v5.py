"""Supplementary Figure S3 (v5): full estimator x instrument-set grid for the primary outcomes and FinnGen (moved from main Figure 2c).
Cells: estimate (95% CI); colour: z. Primary-test markers only on the three primary IVW cells (** P < 0.0083); other cells are secondary and unmarked."""
from figstyle import *
import numpy as np, pandas as pd
R = pd.read_csv('step6v4_mr_results_R.csv'); MV = pd.read_csv('step7_mvmr_results.csv'); MV = MV[MV.instrument_set == 'all 330 variants']
g = lambda code, an, m: R[(R.code == code) & (R.analysis == an) & (R.method == m)]
fig = plt.figure(figsize=(W, 5.2))
# ---- c  primary outcomes: analysis x outcome matrix with estimates and 95% CIs in the cells
ax = fig.add_axes([0.27, 0.08, 0.66, 0.8])
rows = [('r2<0.001 (primary)', 'IVW (MRE, t)', 'IVW, r² < 0.001 (primary)'), ('r2<0.001 (primary)', 'MR-Egger (t)', 'MR-Egger'), ('r2<0.001 (primary)', 'Weighted median', 'Weighted median'),
        ('r2<0.001 (primary)', 'Weighted mode (MBE)', 'Weighted mode (MBE)'), ('r2<0.001 (primary)', 'PRESSO', 'MR-PRESSO (outlier-corrected)'),
        ('r2<0.001, Steiger-filtered', 'IVW (MRE, t)', 'IVW, Steiger-filtered'), ('r2<0.001, excluding TERC/TERT regions', 'IVW (MRE, t)', 'IVW, without TERC/TERT regions'),
        ('r2<0.001, excluding ambiguous palindromic', 'IVW (MRE, t)', 'IVW, without ambiguous palindromes'),
        ('r2<0.1, covariance-aware', 'IVW, correlated instruments (GLS)', 'IVW GLS, r² < 0.1 (correlated)'),
        ('Codd 2021 published MR instrument', 'IVW (MRE, t)', 'IVW, published Codd instrument'),
        ('r2<0.001, excluding instruments with confounder-trait associations', 'IVW (MRE, t)', 'IVW, phenome-screened set (42)'),
        ('MVMR', 'LTL + BMI + LYM', 'MVMR: LTL | BMI, lymphocytes')]
cols = [('S_mvpEU', 'SBP, MVP'), ('D_mvpEU', 'DBP, MVP'), ('H_mvpEU', 'Hypertension, MVP'), ('H_fg12', 'Hypertension, FinnGen R12')]
Z = np.full((len(rows), len(cols)), np.nan); T = [[''] * len(cols) for _ in rows]
for i, (an, m, _) in enumerate(rows):
    for j, (c, _) in enumerate(cols):
        if an == 'MVMR':
            x = MV[(MV.outcome == c) & (MV.model == m) & (MV.exposure == 'LTL') & (MV.method == 'MV-IVW (random, t)')]
            x = x.assign(nsnp=x.nsnp)
        else:
            x = R[(R.code == c) & (R.analysis == an) & (R.method.str.startswith('MR-PRESSO') if m == 'PRESSO' else (R.method == m))]
        if x.empty: T[i][j] = '–'; continue
        x = x.iloc[0]; Z[i, j] = x.b / x.se
        if c.startswith('H_'): txt = f'OR {np.exp(x.b):.2f}\n({np.exp(x.lo):.2f} to {np.exp(x.hi):.2f})'
        else: txt = f'{x.b:.3f}\n({x.lo:.3f} to {x.hi:.3f})'
        T[i][j] = txt.replace('-', '−') + ('**' if (i == 0 and c != 'H_fg12' and x.p < 0.05 / 6) else '')
im = ax.imshow(Z, cmap=DIV, vmin=-8, vmax=8, aspect='auto')
for i in range(len(rows)):
    for j in range(len(cols)): ax.text(j, i, T[i][j], ha='center', va='center', fontsize=5.6, linespacing=1.0, color='white' if abs(np.nan_to_num(Z[i, j])) > 5.5 else INK)
ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[2] for r in rows], fontsize=5.5); ax.set_xticks(range(len(cols))); ax.set_xticklabels([c[1] for c in cols], fontsize=5.8); ax.xaxis.tick_top(); ax.tick_params(length=0)
for yy in [4.5, 8.5]: ax.axhline(yy, color='white', lw=2)
for sp in ax.spines.values(): sp.set_visible(False)
for j, (c, _) in enumerate(cols):
    iv_ = g(c, 'r2<0.001 (primary)', 'IVW (MRE, t)').iloc[0]; e_ = g(c, 'r2<0.001 (primary)', 'MR-Egger (t)').iloc[0]
    ax.text(j, len(rows) - 0.3, f'I² {iv_.I2:.0f}%; Egger intercept P {e_.intercept_p:.3f}', ha='center', va='top', fontsize=4.9, color=RED if e_.intercept_p < 0.05 else INK2)
ax.set_ylim(len(rows) + 0.3, -0.5)
ax.set_title('All estimators and instrument sets: estimate (95% CI) per SD longer LTL; colour = z', pad=14)
cax = ax.inset_axes([1.015, 0.25, 0.018, 0.5]); cb = fig.colorbar(im, cax=cax); cb.outline.set_visible(False); cb.ax.tick_params(labelsize=5); cb.set_label('z', fontsize=5.5)

fig.text(0.01, 0.0, 'SBP and DBP: SD units (MVP inverse-normal scale); hypertension: odds ratio. ** marks only the primary IVW tests below 0.05/6; all other cells are secondary analyses without multiplicity correction. MVMR: 330 variants.', fontsize=5.4, color=MUTED)
save(fig, 'FigS3_mr_full_grid')
