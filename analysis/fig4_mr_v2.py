from figstyle import *
from matplotlib.colors import LinearSegmentedColormap
import numpy as np, pandas as pd
from scipy import stats
from matplotlib.colors import Normalize
from matplotlib.patches import Polygon, Rectangle
from matplotlib.lines import Line2D

F = pd.read_csv('step3_mr_forward_results.csv'); S = pd.read_csv('step3_mr_forward_snp_level.csv'); LOO = pd.read_csv('step3_mr_forward_leave_one_out.csv')
RV = pd.read_csv('step3_mr_reverse_results.csv'); RS = pd.read_csv('step3_mr_reverse_snp_level.csv'); IV = pd.read_csv('step3_instruments_LTL.csv')
LED = pd.read_csv('step3_reverse_instrument_ledger.csv'); W_ = pd.read_csv('kp_ltl_leads_wide.csv')

fig = plt.figure(figsize=(W, 9.4))
gs = fig.add_gridspec(4, 6, height_ratios=[0.62, 1.15, 1.0, 1.05], hspace=0.72, wspace=1.05)

# ---- a  instrument ledger (shrinking bands)
ax = fig.add_subplot(gs[0, :]); ax.axis('off')
wide = W_[W_.varId.str.count(':') == 3]
n0 = len(wide); n1 = int((wide.L_codd_p < 5e-8).sum()); n2 = int(((wide.L_codd_p < 5e-8) & (np.minimum(wide.maf, 1 - wide.maf) >= 0.01) & wide.rsid.notna()).sum()); n3 = len(IV)
outs = [('MVP EUR', 'S_mvpEU'), ('FinnGen', 'H_fg'), ('Genes & Health', 'S_gh'), ('Biobank Japan', 'S_bbj')]
stages = [('KP European LTL\nlead variants', n0), ('P < 5×10$^{-8}$\nin Codd 2021', n1), ('MAF ≥ 1%', n2), ('clumped: >1 Mb and\nr² < 0.05 (1000G EUR)', n3)]
xs = np.linspace(0.02, 0.62, len(stages)); maxn = n0; H = 0.8
for i, (lab, n) in enumerate(stages):
    h = H * n / maxn
    ax.add_patch(Rectangle((xs[i], 0.5 - h / 2), 0.07, h, color=BLUE if i == len(stages) - 1 else SEQ[1 + i], lw=0, transform=ax.transAxes))
    ax.text(xs[i] + 0.035, 0.5 + h / 2 + 0.04, f'{n}', ha='center', va='bottom', fontsize=7.5, fontweight='bold', color=INK, transform=ax.transAxes)
    ax.text(xs[i] + 0.035, 0.5 - H / 2 - 0.05, lab, ha='center', va='top', fontsize=5.8, color=INK2, transform=ax.transAxes)
    if i < len(stages) - 1:
        h2 = H * stages[i + 1][1] / maxn
        ax.add_patch(Polygon([[xs[i] + 0.07, 0.5 - h / 2], [xs[i + 1], 0.5 - h2 / 2], [xs[i + 1], 0.5 + h2 / 2], [xs[i] + 0.07, 0.5 + h / 2]], closed=True, color=SEQ[0], alpha=0.8, lw=0, transform=ax.transAxes))
hN = H * n3 / maxn; x0 = xs[-1] + 0.07
for k, (lab, code) in enumerate(outs):
    n = int(F[(F.code == code) & (F.method == 'IVW (MRE)')].nsnp.iloc[0]); yk = 0.88 - k * 0.25; h = H * n / maxn * 0.22
    ax.add_patch(Polygon([[x0, 0.5 - hN / 2 + k * hN / 4], [0.8, yk - h / 2], [0.8, yk + h / 2], [x0, 0.5 - hN / 2 + (k + 1) * hN / 4]], closed=True, color=SEQ[0], alpha=0.8, lw=0, transform=ax.transAxes))
    ax.add_patch(Rectangle((0.8, yk - h / 2), 0.03, h, color=CAT[k], lw=0, transform=ax.transAxes))
    ax.text(0.84, yk, f'{lab}: {n} harmonised', va='center', fontsize=6.2, color=INK, transform=ax.transAxes)
ax.text(0.0, 1.12, 'a', fontsize=10, fontweight='bold', transform=ax.transAxes); ax.text(0.03, 1.12, 'Instrument ledger (forward MR)', fontsize=7.5, fontweight='bold', transform=ax.transAxes)

# ---- b  SNP scatter with marginals (MVP EUR SBP)
s = S[S.code == 'S_mvpEU'].copy(); sg = np.sign(s.bx); s['x'] = s.bx * sg; s['y'] = s.by * sg
sub = gs[1, :3].subgridspec(4, 4, hspace=0.05, wspace=0.05)
ax = fig.add_subplot(sub[1:, :3]); axt = fig.add_subplot(sub[0, :3], sharex=ax); axr = fig.add_subplot(sub[1:, 3], sharey=ax)
zc = np.clip(s.steiger_z, 0, 60)
ax.errorbar(s.x, s.y, xerr=1.96 * s.sx, yerr=1.96 * s.sy, fmt='none', ecolor='#e6e3db', elinewidth=0.5, zorder=1)
sc = ax.scatter(s.x, s.y, c=zc, cmap=LinearSegmentedColormap.from_list('z', SEQ[1:]), vmin=0, vmax=60, s=np.where(s.presso_outlier, 22, 12), edgecolor=np.where(s.presso_outlier, RED, 'white'), lw=np.where(s.presso_outlier, 0.9, 0.3), zorder=3)
xx = np.linspace(0, s.x.max() * 1.05, 50)
ivw = F[(F.code == 'S_mvpEU') & (F.method == 'IVW (MRE)')].iloc[0]
ax.fill_between(xx, ivw.lo * xx, ivw.hi * xx, color=BLUE, alpha=0.13, lw=0)
for mth, col, ls in [('IVW (MRE)', BLUE, '-'), ('Weighted median', AQUA, '--'), ('Weighted mode', YELLOW, '-.'), ('MR-PRESSO (outlier-corrected)', VIOLET, ':')]:
    r = F[(F.code == 'S_mvpEU') & (F.method == mth)].iloc[0]; ax.plot(xx, r.b * xx, color=col, ls=ls, lw=1.2, label=mth.replace(' (MRE)', '').replace(' (outlier-corrected)', ''))
e = F[(F.code == 'S_mvpEU') & (F.method == 'MR-Egger')].iloc[0]; ax.plot(xx, e.intercept + e.b * xx, color=ORANGE, lw=1.2, label='MR-Egger')
for g, off in [('TERT', (5, -9)), ('MYNN', (4, 5)), ('TERC', (4, -9))]:
    for r in s[s.gene.str.contains(g, na=False)].itertuples(): ax.annotate(f'{r.gene}', (r.x, r.y), xytext=off, textcoords='offset points', fontsize=5.8, color=INK, fontweight='bold')
ax.axhline(0, color=AXIS, lw=0.6); ax.set_xlabel('SNP → LTL (Codd 2021), oriented positive'); ax.set_ylabel('SNP → SBP (MVP European)')
axt.legend(*ax.get_legend_handles_labels(), loc='upper right', fontsize=5.3, ncol=2, handlelength=1.6, columnspacing=0.8, bbox_to_anchor=(1.35, 1.05))
axt.hist(s.x, bins=30, color=SEQ[2], lw=0); axt.axis('off'); axr.hist(s.y, bins=30, orientation='horizontal', color=SEQ[2], lw=0); axr.axis('off')
cax = ax.inset_axes([0.58, 0.1, 0.38, 0.03]); cb = fig.colorbar(sc, cax=cax, orientation='horizontal'); cb.outline.set_visible(False); cb.ax.tick_params(labelsize=5, pad=1); cax.set_title('Steiger z; red ring = outlier', fontsize=5.0, color=INK2, fontweight='normal', pad=1, loc='left')
axt.set_title('SNP-level evidence, SBP', loc='left'); axt.text(-0.13, 1.1, 'b', transform=axt.transAxes, fontsize=10, fontweight='bold')

# ---- c  radial MR
ax = fig.add_subplot(gs[1, 3:])
ratio = s.by / s.bx; wj = s.bx ** 2 / s.sy ** 2; xr = np.sqrt(wj); yr = ratio * xr
Qj = wj * (ratio - ivw.b) ** 2; thr = stats.chi2.ppf(0.95, 1)
ax.scatter(xr, yr, s=12, c=np.where(Qj > thr, ORANGE, BLUE), edgecolor='white', lw=0.3, zorder=3, alpha=0.9)
xm = xr.max() * 1.05; xx = np.linspace(0, xm, 50)
ax.plot(xx, ivw.b * xx, color=INK, lw=1.1); ax.fill_between(xx, ivw.lo * xx, ivw.hi * xx, color=MUTED, alpha=0.2, lw=0)
tt = np.linspace(-np.pi / 2, np.pi / 2, 200); R = xm * 0.98
for bb, lab in [(0, '0'), (ivw.b, f'IVW {ivw.b:.2f}'), (0.5, '0.5')]:
    ang = np.arctan(bb); ax.plot([R * np.cos(ang) * 0.97, R * np.cos(ang)], [R * np.sin(ang) * 0.97, R * np.sin(ang)], color=INK2, lw=0.7)
    ax.text(R * np.cos(ang) * 1.02, R * np.sin(ang) * 1.02, lab, fontsize=5.6, color=INK2, va='center')
ax.plot(R * np.cos(tt[80:120]) , R * np.sin(tt[80:120]), color=AXIS, lw=0.6)
for r, xv, yv, q in zip(s.itertuples(), xr, yr, Qj):
    if q > 40 or (isinstance(r.gene, str) and ('MYNN' in r.gene or 'TERT' in r.gene)): ax.annotate(r.gene if isinstance(r.gene, str) else r.rsid, (xv, yv), xytext=(3, 2), textcoords='offset points', fontsize=5.5, color=ORANGE)
ax.set_ylim(-11, 15); ax.axhline(0, color=AXIS, lw=0.5); ax.set_xlabel('√ weight  (|β$_{X}$| / SE$_{Y}$)'); ax.set_ylabel('ratio estimate × √ weight')
ax.text(0.02, 0.97, f'orange: Q$_j$ > χ²$_{{1,0.95}}$ ({int((Qj > thr).sum())} of {len(Qj)})\nCochran Q = {ivw.Q:.0f}, I² = {ivw.I2:.0f}%', transform=ax.transAxes, va='top', fontsize=5.8, color=INK2)
ax.set_title('Radial MR (IVW), SBP'); tag(ax, 'c')

# ---- d  funnel of ratio estimates
ax = fig.add_subplot(gs[2, :2])
se_r = np.sqrt(s.sy ** 2 / s.bx ** 2); prec = 1 / se_r
ax.scatter(ratio, prec, s=10, c=np.where(s.presso_outlier, RED, BLUE), alpha=0.8, lw=0)
for mth, col in [('IVW (MRE)', INK), ('MR-Egger', ORANGE), ('Weighted median', AQUA)]:
    r = F[(F.code == 'S_mvpEU') & (F.method == mth)].iloc[0]; ax.axvline(r.b, color=col, lw=1, ls='-' if 'IVW' in mth else '--')
ax.set_xlim(-1.2, 1.6); ax.set_yscale('log'); ax.minorticks_off(); ax.set_xlabel('per-SNP ratio (SD SBP / SD LTL)'); ax.set_ylabel('precision (1/SE)')
ax.legend(handles=[Line2D([],[],color=INK,label='IVW'),Line2D([],[],color=ORANGE,ls='--',label='Egger'),Line2D([],[],color=AQUA,ls='--',label='W-median'),Line2D([],[],marker='o',ls='',color=RED,ms=3,label='outlier')],fontsize=5.2,loc='upper left',handlelength=1.4)
ax.set_title('Funnel (asymmetry → pleiotropy)'); tag(ax, 'd')

# ---- e  leave-one-out trajectory
ax = fig.add_subplot(gs[2, 2:4])
for code, col, lab in [('S_mvpEU', BLUE, 'SBP'), ('D_mvpEU', AQUA, 'DBP'), ('H_mvpEU', ORANGE, 'HTN (log OR)')]:
    l = LOO[LOO.code == code].sort_values('b').reset_index(drop=True); x = np.linspace(0, 1, len(l))
    ax.plot(x, l.b - 1.96 * l.se, color=col, lw=0.5, ls=(0,(1,1))); ax.plot(x, l.b + 1.96 * l.se, color=col, lw=0.5, ls=(0,(1,1))); ax.plot(x, l.b, color=col, lw=1.4, label=lab)
    full = F[(F.code == code) & (F.method == 'IVW (MRE)')].iloc[0].b; ax.axhline(full, color=col, lw=0.6, ls=':')
    ex = F[(F.code == code) & (F.method == 'IVW, excluding TERC & TERT loci')].iloc[0].b; ax.scatter(1.04, ex, marker='<', color=col, s=18, clip_on=False)
ax.axhline(0, color=AXIS, lw=0.6); ax.set_xlim(0, 1.08); ax.set_xticks([0, 0.5, 1]); ax.set_xticklabels(['min', 'median', 'max'])
ax.set_xlabel('leave-one-SNP-out estimates, ranked'); ax.set_ylabel('IVW estimate (95% CI)'); ax.set_ylim(-0.06, 0.40); ax.legend(fontsize=5.6, loc='upper left', ncol=3, columnspacing=0.8, handlelength=1.2)
ax.text(0.02, 0.87, '◄ IVW without TERC/TERT loci; dotted: 95% CI', transform=ax.transAxes, fontsize=5.0, color=INK2, ha='left', va='top')
ax.set_title('Leave-one-out stability'); tag(ax, 'e')

# ---- e2 instrument strength
ax = fig.add_subplot(gs[2, 4:])
Fv = IV.F.values; bins = np.exp(np.linspace(np.log(20), np.log(Fv.max() * 1.1), 25))
ax.hist(Fv, bins=bins, color=SEQ[3], lw=0.3, edgecolor='white'); ax.set_xscale('log'); ax.minorticks_off()
ax.axvline(10, color=RED, lw=0.8, ls='--'); ax.text(10.5, ax.get_ylim()[1] * 0.9, 'F = 10', color=RED, fontsize=5.5)
ax.set_xlim(8, Fv.max() * 1.2); ax.set_xticks([10, 30, 100, 300, 1000]); ax.set_xticklabels(['10', '30', '100', '300', '1000'])
r2 = (IV.bx ** 2 / (IV.bx ** 2 + 464716 * IV.sx ** 2)).sum()
ax.text(0.97, 0.95, f'{len(IV)} instruments\nmedian F = {np.median(Fv):.0f}, min {Fv.min():.0f}\nΣr² (LTL) = {r2*100:.1f}%', transform=ax.transAxes, ha='right', va='top', fontsize=5.8, color=INK)
ax.set_xlabel('per-instrument F statistic (log)'); ax.set_ylabel('instruments'); ax.set_title('Instrument strength'); tag(ax, 'f')

# ---- f  method × outcome heatmap
ax = fig.add_subplot(gs[2:, 4:]) if False else fig.add_subplot(gs[3, :4])
rows = [('S_mvpEU', 'SBP · MVP EUR'), ('D_mvpEU', 'DBP · MVP EUR'), ('H_mvpEU', 'HTN · MVP EUR'), ('H_fg', 'HTN · FinnGen'), ('S_bbj', 'SBP · BBJ'), ('D_bbj', 'DBP · BBJ'),
        ('S_mvpAA', 'SBP · MVP AFR-Am'), ('S_mvpHS', 'SBP · MVP Hisp'), ('S_gh', 'SBP · G&H (SAS)'), ('D_gh', 'DBP · G&H (SAS)')]
meths = [('IVW (MRE)', 'IVW'), ('MR-Egger', 'Egger'), ('Weighted median', 'W-median'), ('Weighted mode', 'W-mode'), ('MR-PRESSO (outlier-corrected)', 'PRESSO'), ('IVW, Steiger-filtered', 'Steiger'), ('IVW, excluding TERC & TERT loci', 'no TERC/\nTERT')]
Z = np.full((len(rows), len(meths)), np.nan); TXT = [[''] * len(meths) for _ in rows]
for i, (c, _) in enumerate(rows):
    for j, (m, _) in enumerate(meths):
        r = F[(F.code == c) & (F.method == m)].iloc[0]; Z[i, j] = r.b / r.se
        TXT[i][j] = (f'{np.exp(r.b):.2f}' if c.startswith('H_') else f'{r.b:.2f}') + ('**' if r.p < 0.05 / 6 else ('*' if r.p < 0.05 else ''))
im = ax.imshow(Z, cmap=DIV, vmin=-8, vmax=8, aspect='auto')
for i in range(len(rows)):
    for j in range(len(meths)): ax.text(j, i, TXT[i][j], ha='center', va='center', fontsize=5.6, color='white' if abs(Z[i, j]) > 5 else INK)
ax.set_xticks(range(len(meths))); ax.set_xticklabels([m[1] for m in meths], fontsize=6); ax.xaxis.tick_top(); ax.tick_params(length=0)
ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[1] for r in rows], fontsize=6)
for sp in ax.spines.values(): sp.set_visible(False)
ax.axhline(2.5, color='white', lw=2); ax.axhline(3.5, color='white', lw=2); ax.axhline(7.5, color='white', lw=2)
# side columns: egger intercept P and I2
for i, (c, _) in enumerate(rows):
    e = F[(F.code == c) & (F.method == 'MR-Egger')].iloc[0]; iv_ = F[(F.code == c) & (F.method == 'IVW (MRE)')].iloc[0]
    ax.text(len(meths) - 0.3, i, ('<0.001' if e.intercept_p < 0.001 else f'{e.intercept_p:.3f}'), va='center', fontsize=5.6, color=RED if e.intercept_p < 0.05 else INK2)
    ax.text(len(meths) + 0.75, i, f'{iv_.I2:.0f}%', va='center', fontsize=5.6, color=INK2)
ax.text(len(meths) - 0.3, -0.7, 'Egger\nint. P', fontsize=5.6, color=INK, va='bottom'); ax.text(len(meths) + 0.75, -0.7, 'I²', fontsize=5.6, color=INK, va='bottom')
ax.set_xlim(-0.5, len(meths) + 1.5)
cax = ax.inset_axes([0.0, -0.1, 0.35, 0.03]); cb = fig.colorbar(im, cax=cax, orientation='horizontal'); cb.outline.set_visible(False); cb.ax.tick_params(labelsize=5.3)
cb.set_label('z = estimate / SE; cell = estimate (OR for HTN); ** P < 0.0083 (Bonferroni), * P < 0.05', fontsize=5.5)
ax.set_title('Evidence matrix: estimator × outcome (per SD longer LTL)', pad=22); tag(ax, 'g', x=-0.2, y=1.12)

# ---- g  reverse MR density
ax = fig.add_subplot(gs[3, 4:])
rs = RS[(RS.exposure == 'SBP') & (RS.outcome.str.startswith('LTL (Codd'))].copy(); sg = np.sign(rs.bx); rs['x'] = rs.bx * sg; rs['y'] = rs.by * sg
ax.errorbar(rs.x, rs.y, yerr=1.96*rs.sy, fmt='none', ecolor='#ebe8e1', elinewidth=0.4, zorder=1); ax.scatter(rs.x, rs.y, s=6, color=SEQ[3], alpha=0.7, lw=0, zorder=2)
xx = np.linspace(0, rs.x.max(), 20)
for ex, col in [('SBP', BLUE), ('DBP', AQUA), ('Hypertension', ORANGE)]:
    r = RV[(RV.exposure == ex) & (RV.outcome.str.startswith('LTL (Codd')) & (RV.method == 'IVW (MRE)')].iloc[0]
    from scipy.stats import t as _t; q = _t.ppf(0.975, int(r.nsnp) - 1)
    if ex == 'SBP': ax.fill_between(xx, (r.b - q * r.se) * xx, (r.b + q * r.se) * xx, color=col, alpha=0.2, lw=0)
    ax.plot(xx, r.b * xx, color=col, lw=1.1, label=f'{ex}: {r.b:+.3f} ({r.b-q*r.se:+.3f}, {r.b+q*r.se:+.3f})')
ax.axhline(0, color=AXIS, lw=0.6); ax.set_xlabel('SNP → SBP (MVP European)'); ax.set_ylabel('SNP → LTL (Codd 2021)')
ax.legend(fontsize=5.2, loc='upper left', title='IVW, SD LTL per unit exposure', title_fontsize=5.2, framealpha=0.9, frameon=True, edgecolor='none')
ax.set_ylim(-0.03, 0.04); ax.set_title('Reverse MR: BP → LTL'); tag(ax, 'h')
fig.text(0.01, 0.0, 'Exposure and outcome samples do not overlap (forward: UK Biobank LTL → MVP/FinnGen/Genes & Health/Biobank Japan; reverse: MVP → UK Biobank). IVW uses multiplicative random effects with t(k−1) inference. '
         'Clumping: >1 Mb and r² < 0.05 in 1000G EUR (Ensembl reports r² ≥ 0.05 only); stricter r² < 0.001 requires a local reference.', fontsize=5.5, color=MUTED, wrap=True)
save(fig, 'Fig4_mendelian_randomization')
