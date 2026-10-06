"""Figure 4 (v5): (a) ancestry-specific associations of the LTL-shortening alleles as interval bands with N and SE source (no pooled row);
(b) transportability of the forward MR estimate (European SNP-LTL weights) across outcome ancestries. Former panels b, c moved to Supplementary Figure S5."""
from figstyle import *
import numpy as np, pandas as pd
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
A = pd.read_csv('kp_snp_lookup_by_ancestry_checked.csv'); R = pd.read_csv('step6v4_mr_results_R.csv')
ancs = ['EU', 'EA', 'SA', 'AF', 'AA', 'HS']; traits = [('SBP', 'SBP (SD)'), ('DBP', 'DBP (SD)'), ('HYPERTENSION', 'Hypertension (log OR)')]
snps = [('rs10936599', 'TERC rs10936599-T'), ('rs2736100', 'TERT rs2736100-A')]
fig = plt.figure(figsize=(W, 7.4))
X0, X1, Y0, Y1 = 0.14, 0.985, 0.45, 0.9; gx, gy = 0.03, 0.07; pw = (X1 - X0 - 2 * gx) / 3; ph = (Y1 - Y0 - gy) / 2
lims = {'SBP': (-0.065, 0.035), 'DBP': (-0.065, 0.035), 'HYPERTENSION': (-0.12, 0.1)}
for r_, (s, slab) in enumerate(snps):
    for c_, (t, tl) in enumerate(traits):
        ax = fig.add_axes([X0 + c_ * (pw + gx), Y1 - (r_ + 1) * ph - r_ * gy, pw, ph])
        for i, a in enumerate(ancs):
            x = A[(A.rs == s) & (A.ancestry == a) & (A.trait == t)]
            if x.empty: ax.text(0.03, i, 'no record', transform=ax.get_yaxis_transform(), ha='left', va='center', fontsize=4.8, color=MUTED); continue
            x = x.iloc[0]; b, se = x.beta_alt, x.se_p; lo, hi = b - 1.96 * se, b + 1.96 * se; col = ANC[a]
            if not bool(x.qc_pass):
                ax.add_patch(Rectangle((lo, i - 0.28), hi - lo, 0.56, fc='none', ec=MUTED, lw=0.6, hatch='////'))
                ax.text(0.98, i, 'excluded*', transform=ax.get_yaxis_transform(), ha='right', va='center', fontsize=4.6, color=RED); continue
            ax.add_patch(Rectangle((lo, i - 0.28), hi - lo, 0.56, color=col, alpha=0.3, lw=0))
            ax.add_patch(Rectangle((b - 0.674 * se, i - 0.28), 1.348 * se, 0.56, color=col, alpha=0.75, lw=0))
            ax.plot([b, b], [i - 0.36, i + 0.36], color=INK, lw=1.0)
            nn = f'{x.n / 1e6:.1f}M' if x.n >= 1e6 else f'{x.n / 1e3:.0f}k'
            ax.text(0.99, i - 0.02, nn + ('†' if not bool(x.consistent) else ''), transform=ax.get_yaxis_transform(), ha='right', va='center', fontsize=4.5, color=INK2)
        ax.axvline(0, color=INK2, lw=0.6); ax.set_xlim(*lims[t]); ax.set_ylim(len(ancs) - 0.5, -0.6)
        ax.set_yticks(range(len(ancs))); ax.set_yticklabels([ANCLAB[a] for a in ancs] if c_ == 0 else [], fontsize=5.6); ax.tick_params(axis='y', length=0)
        for yy in np.arange(0.5, len(ancs) - 0.5): ax.axhline(yy, color=GRID, lw=0.4)
        ax.set_title(tl if r_ == 0 else '', fontsize=6.3); clean(ax, left=False); ax.tick_params(axis='x', labelsize=5.5)
        if c_ == 0: ax.text(-0.52, 0.5, slab, transform=ax.transAxes, rotation=90, ha='center', va='center', fontsize=6.5, fontweight='bold', color=INK)
fig.text(X0, 0.945, 'Per-allele associations of the LTL-shortening allele by ancestry: 95% CI (light), 50% interval (dark), estimate (notch); right: N in the record', fontsize=6.6, fontweight='bold', color=INK)
fig.text(0.02, 0.945, 'a', fontsize=10, fontweight='bold')
# ---- b  transportability of the forward MR estimate
ax = fig.add_axes([0.2, 0.09, 0.62, 0.27])
rows = [('MVP European (primary)', {'SBP': 'S_mvpEU', 'DBP': 'D_mvpEU', 'HYP': 'H_mvpEU'}, 'EU'), ('FinnGen R12 (secondary)', {'HYP': 'H_fg12'}, 'EU'),
        ('Biobank Japan', {'SBP': 'S_bbj', 'DBP': 'D_bbj'}, 'EA'), ('Genes & Health (South Asian)', {'SBP': 'S_gh', 'DBP': 'D_gh'}, 'SA'),
        ('MVP African American', {'SBP': 'S_mvpAA', 'DBP': 'D_mvpAA', 'HYP': 'H_mvpAA'}, 'AA'), ('MVP Hispanic', {'SBP': 'S_mvpHS', 'DBP': 'D_mvpHS', 'HYP': 'H_mvpHS'}, 'HS')]
off = {'SBP': -0.26, 'DBP': 0.0, 'HYP': 0.26}; mk = {'SBP': 'SBP', 'DBP': 'DBP', 'HYP': 'Hyp.'}
for i, (lab, codes, a) in enumerate(rows):
    for tr, code in codes.items():
        x = R[(R.code == code) & (R.analysis == 'r2<0.001 (primary)') & (R.method == 'IVW (MRE, t)')]
        if x.empty: continue
        x = x.iloc[0]; y = i + off[tr]; col = ANC[a]
        ax.add_patch(Rectangle((x.lo, y - 0.1), x.hi - x.lo, 0.2, color=col, alpha=0.35 if tr != 'HYP' else 0.2, lw=0, hatch=None if tr != 'HYP' else '....'))
        ax.plot([x.b] * 2, [y - 0.12, y + 0.12], color=INK, lw=1.0)
        ax.text(-0.33, y, mk[tr], fontsize=4.8, va='center', ha='left', color=INK2)
        ax.text(0.62, y, f'{x.b:+.3f} ({x.lo:+.3f}, {x.hi:+.3f})  k={int(x.nsnp)}  P={x.p:.1e}'.replace('-', '−'), fontsize=4.7, va='center', color=INK2)
ax.axvline(0, color=INK2, lw=0.6); ax.set_xlim(-0.34, 0.6); ax.set_ylim(len(rows) - 0.5, -0.6)
ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0] for r in rows], fontsize=5.8); ax.tick_params(axis='y', length=0); clean(ax, left=False)
for yy in np.arange(0.5, len(rows) - 0.5): ax.axhline(yy, color=GRID, lw=0.5)
ax.set_xlabel('IVW estimate per SD longer LTL (SBP, DBP: SD; hypertension: log OR)', fontsize=6)
ax.set_title('Transportability: forward MR with European SNP–LTL weights applied to other outcome ancestries (95% CI)', fontsize=6.6); tag(ax, 'b', x=-0.29, y=1.03)
fig.text(0.01, -0.012, 'Look-ups: Knowledge Portal ancestry-specific meta-analyses (accessed 29 Sep 2026), which can include UK Biobank-based studies; SE derived from P. † reported SE and P in the record disagree '
         '(P-derived SE shown). * rs10936599 Hispanic SBP record excluded because its reported β, SE and P are mutually inconsistent. No cross-ancestry pooled estimate is shown. '
         'Panel b: non-European estimates test transportability, not ancestry-matched causation.', fontsize=5.3, color=MUTED, wrap=True)
save(fig, 'Fig4_cross_ancestry')
