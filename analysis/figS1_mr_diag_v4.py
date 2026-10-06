"""Supplementary Figure S1 (v4): MR diagnostics for the primary SBP analysis and leave-one-out for all primary outcomes."""
from figstyle import *
import numpy as np, pandas as pd
from scipy import stats
from matplotlib.lines import Line2D
R = pd.read_csv('step6v4_mr_results_R.csv'); LOO = pd.read_csv('step6_v4_forward_leave_one_out.csv'); MV = pd.read_csv('step6_mvmr_results.csv')
g = lambda code, an, m: R[(R.code == code) & (R.analysis == an) & (R.method == m)].iloc[0]
fig = plt.figure(figsize=(W, 9.4)); gs = fig.add_gridspec(3, 2, hspace=0.5, wspace=0.42, height_ratios=[1, 1, 1.15])
s = pd.read_csv('mr_v4_inputs/snp_level_S_mvpEU.csv'); iv = g('S_mvpEU', 'r2<0.001 (primary)', 'IVW (MRE, t)')
# a radial
ax = fig.add_subplot(gs[0, 0]); sg = np.sign(s.bx); ratio = s.by / s.bx; wj = s.bx ** 2 / s.sy ** 2; xr = np.sqrt(wj); yr = ratio * xr
Qj = wj * (ratio - iv.b) ** 2; thr = stats.chi2.ppf(0.95, 1)
ax.scatter(xr, yr, s=11, c=np.where(Qj > thr, ORANGE, BLUE), edgecolor='white', lw=0.3, zorder=3)
xm = xr.max() * 1.05; xx = np.linspace(0, xm, 50); ax.plot(xx, iv.b * xx, color=INK, lw=1.1); ax.fill_between(xx, iv.lo * xx, iv.hi * xx, color=MUTED, alpha=0.2, lw=0)
for r, xv, yv, q in zip(s.itertuples(), xr, yr, Qj):
    if q > 40: ax.annotate(r.gene if isinstance(r.gene, str) else r.rsid, (xv, yv), xytext=(3, 2), textcoords='offset points', fontsize=5.3, color=ORANGE)
ax.axhline(0, color=AXIS, lw=0.5); ax.set_xlabel('√ weight (|β$_X$| / SE$_Y$)'); ax.set_ylabel('ratio estimate × √ weight')
ax.text(0.02, 0.97, f'orange: Q$_j$ > χ²$_{{1,0.95}}$ ({int((Qj > thr).sum())} of {len(Qj)}); I² = {iv.I2:.0f}%', transform=ax.transAxes, va='top', fontsize=5.6, color=INK2)
ax.set_title('Radial IVW, SBP (MVP EUR)'); tag(ax, 'a')
# b funnel
ax = fig.add_subplot(gs[0, 1]); se_r = np.abs(s.sy / s.bx)
ax.scatter(ratio, 1 / se_r, s=10, color=BLUE, alpha=0.8, lw=0)
for m, col, ls in [('IVW (MRE, t)', INK, '-'), ('MR-Egger (t)', ORANGE, '--'), ('Weighted median', AQUA, '--')]:
    ax.axvline(g('S_mvpEU', 'r2<0.001 (primary)', m).b, color=col, lw=1, ls=ls)
ax.set_xlim(-1.2, 1.6); ax.set_yscale('log'); ax.minorticks_off(); ax.set_xlabel('per-variant ratio (SD SBP / SD LTL)'); ax.set_ylabel('precision (1/SE)')
ax.legend(handles=[Line2D([], [], color=INK, label='IVW'), Line2D([], [], color=ORANGE, ls='--', label='MR-Egger'), Line2D([], [], color=AQUA, ls='--', label='weighted median')], fontsize=5.3, loc='upper left')
ax.set_title('Funnel, SBP'); tag(ax, 'b')
# c LOO
ax = fig.add_subplot(gs[1, 0])
for code, col, lab in [('S_mvpEU', BLUE, 'SBP'), ('D_mvpEU', AQUA, 'DBP'), ('H_mvpEU', ORANGE, 'Hypertension (log OR)')]:
    l = LOO[LOO.code == code].sort_values('b').reset_index(drop=True); x = np.linspace(0, 1, len(l))
    ax.plot(x, l.b - 1.96 * l.se, color=col, lw=0.5, ls=(0, (1, 1))); ax.plot(x, l.b + 1.96 * l.se, color=col, lw=0.5, ls=(0, (1, 1))); ax.plot(x, l.b, color=col, lw=1.4, label=lab)
    ax.axhline(g(code, 'r2<0.001 (primary)', 'IVW (MRE, t)').b, color=col, lw=0.6, ls=':')
ax.axhline(0, color=AXIS, lw=0.6); ax.set_xticks([0, 0.5, 1]); ax.set_xticklabels(['min', 'median', 'max']); ax.set_xlabel('leave-one-variant-out IVW estimates, ranked'); ax.set_ylabel('estimate (95% CI)')
ax.set_ylim(top=ax.get_ylim()[1] + 0.05); ax.legend(fontsize=5.5, loc='upper left', ncol=3); ax.set_title('Leave-one-out stability'); tag(ax, 'c')
# d F
ax = fig.add_subplot(gs[1, 1]); Fv = (s.bx / s.sx) ** 2
bins = np.exp(np.linspace(np.log(20), np.log(Fv.max() * 1.1), 25)); ax.hist(Fv, bins=bins, color=SEQ[3], lw=0.3, edgecolor='white'); ax.set_xscale('log'); ax.minorticks_off()
ax.axvline(10, color=RED, lw=0.8, ls='--'); ax.set_xlim(8, Fv.max() * 1.2); ax.set_xticks([10, 30, 100, 300, 1000]); ax.set_xticklabels(['10', '30', '100', '300', '1000'])
r2 = (s.bx ** 2 / (s.bx ** 2 + 464716 * s.sx ** 2)).sum()
ax.text(0.97, 0.95, f'{len(Fv)} instruments\nmedian F = {np.median(Fv):.0f}, min {Fv.min():.0f}\nΣr² (LTL) = {r2*100:.1f}%', transform=ax.transAxes, ha='right', va='top', fontsize=5.8)
ax.set_xlabel('per-instrument F statistic'); ax.set_ylabel('instruments'); ax.set_title('Instrument strength'); tag(ax, 'd')

from matplotlib.colors import LinearSegmentedColormap
V = pd.read_csv('validation/mr_python_vs_R.csv')
# e secondary outcomes
ax = fig.add_subplot(gs[2, 0])
sec = [('H_fg12', 'Hypertension FinnGen R12'), ('S_bbj', 'SBP BBJ'), ('D_bbj', 'DBP BBJ'), ('S_mvpAA', 'SBP MVP AFR-Am'), ('D_mvpAA', 'DBP MVP AFR-Am'), ('H_mvpAA', 'Hypertension MVP AFR-Am'),
       ('S_mvpHS', 'SBP MVP Hisp'), ('D_mvpHS', 'DBP MVP Hisp'), ('H_mvpHS', 'Hypertension MVP Hisp'), ('S_gh', 'SBP G&H (SAS)'), ('D_gh', 'DBP G&H (SAS)')]
ms = [('IVW (MRE, t)', 'IVW'), ('MR-Egger (t)', 'Egger'), ('Weighted median', 'WM'), ('Weighted mode (MBE)', 'Mode')]
Z = np.full((len(sec), len(ms)), np.nan); T = [[''] * len(ms) for _ in sec]
for i, (c, _) in enumerate(sec):
    for j, (m, _) in enumerate(ms):
        x = g(c, 'r2<0.001 (primary)', m); Z[i, j] = x.b / x.se
        T[i][j] = (f'{np.exp(x.b):.2f}' if c.startswith('H_') else f'{x.b:.2f}').replace('-', '−') + ('*' if x.p < 0.05 else '')
ax.imshow(Z, cmap=DIV, vmin=-8, vmax=8, aspect='auto')
for i in range(len(sec)):
    for j in range(len(ms)): ax.text(j, i, T[i][j], ha='center', va='center', fontsize=5.3, color='white' if abs(Z[i, j]) > 5.5 else INK)
for i, (c, _) in enumerate(sec):
    e = g(c, 'r2<0.001 (primary)', 'MR-Egger (t)'); ax.text(len(ms) - 0.35, i, ('<0.001' if e.intercept_p < 0.001 else f'{e.intercept_p:.3f}'), va='center', fontsize=5.2, color=RED if e.intercept_p < 0.05 else INK2)
ax.text(len(ms) - 0.35, -0.8, 'Egger\nint. P', fontsize=5.2, va='bottom')
ax.set_xlim(-0.5, len(ms) + 0.5)
ax.set_yticks(range(len(sec))); ax.set_yticklabels([s_[1] for s_ in sec], fontsize=5.6); ax.set_xticks(range(len(ms))); ax.set_xticklabels([m[1] for m in ms], fontsize=5.8); ax.xaxis.tick_top(); ax.tick_params(length=0)
for sp in ax.spines.values(): sp.set_visible(False)
ax.set_title('Secondary outcomes (* P < 0.05, uncorrected)', pad=18, fontsize=6.8); tag(ax, 'e', x=-0.45, y=1.08)
# f multivariable MR
ax = fig.add_subplot(gs[2, 1])
oc = [('S_mvpEU', 'SBP'), ('D_mvpEU', 'DBP'), ('H_mvpEU', 'Hypertension')]
ex = [('LTL', 'LTL'), ('BMI', 'BMI'), ('LYM', 'Lymphocyte count')]
Z = np.full((len(ex), len(oc)), np.nan); T = [[''] * len(oc) for _ in ex]
for i, (e, _) in enumerate(ex):
    for j, (o, _) in enumerate(oc):
        x = MV[(MV.outcome == o) & (MV.model == 'LTL + BMI + LYM') & (MV.exposure == e) & (MV.method == 'MV-IVW (random, t)')].iloc[0]
        Z[i, j] = x.b / x.se
        T[i][j] = (f'OR {np.exp(x.b):.2f}\n({np.exp(x.lo):.2f} to {np.exp(x.hi):.2f})' if o.startswith('H_') else f'{x.b:.3f}\n({x.lo:.3f} to {x.hi:.3f})').replace('-', '−') + ('*' if x.p < 0.05 else '')
ax.imshow(Z, cmap=DIV, vmin=-8, vmax=8, aspect='auto')
for i in range(len(ex)):
    for j in range(len(oc)): ax.text(j, i, T[i][j], ha='center', va='center', fontsize=5.2, color='white' if abs(Z[i, j]) > 5.5 else INK)
cf = MV[(MV.model == 'LTL + BMI + LYM') & (MV.method == 'MV-IVW (random, t)')].groupby('exposure').condF.first()
ax.set_yticks(range(len(ex))); ax.set_yticklabels([f'{l}\n(conditional F {cf[e]:.0f})' for e, l in ex], fontsize=5.6)
ax.set_xticks(range(len(oc))); ax.set_xticklabels([l + ', MVP' for _, l in oc], fontsize=5.8); ax.xaxis.tick_top(); ax.tick_params(length=0)
for sp in ax.spines.values(): sp.set_visible(False)
ax.set_title(f'Multivariable MR ({int(MV.nsnp.iloc[0])} variants; per SD exposure)', pad=18, fontsize=6.8); tag(ax, 'f', x=-0.35, y=1.08)
save(fig, 'FigS1_mr_diagnostics')
