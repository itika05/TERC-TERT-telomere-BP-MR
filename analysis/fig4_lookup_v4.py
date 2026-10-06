"""Figure 4 (v4): cross-ancestry look-up of rs10936599 and rs2736100 (Knowledge Portal ancestry-specific meta-analyses) and 1000 Genomes allele frequencies."""
from figstyle import *
import numpy as np, pandas as pd, json
from scipy.stats import norm
from matplotlib.patches import Circle, Wedge, Patch
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize

A = pd.read_csv('kp_snp_lookup_by_ancestry_checked.csv'); P = pd.read_csv('step2_cross_ancestry_pooled.csv')
RF = pd.read_csv('ensembl_allele_frequencies.csv')
F = pd.read_csv('step6v4_mr_results_R.csv')
ancs = ['EU', 'SA', 'EA', 'AF', 'AA', 'HS']; traits = [('SBP', 'SBP'), ('DBP', 'DBP'), ('HYPERTENSION', 'Hyper-\ntension')]; snps = [('rs10936599', 'TERC\nrs10936599 T'), ('rs2736100', 'TERT\nrs2736100 A')]

fig = plt.figure(figsize=(W, 6.6))
gs = fig.add_gridspec(2, 6, height_ratios=[1.25, 0.7], hspace=0.62, wspace=1.1)

# ---- a  bubble matrix
ax = fig.add_subplot(gs[0, :4]); norm_ = Normalize(-8, 8)
cols = [(s, t) for s, _ in snps for t, _ in traits]
for j, (s, t) in enumerate(cols):
    for i, a in enumerate(ancs + ['POOL']):
        y = -i
        if a == 'POOL':
            r = P[(P.rs == s) & (P.trait == t)].iloc[0]; b, se, p = r.b_random, r.se_random, r.p_random; ok = True
        else:
            x = A[(A.rs == s) & (A.ancestry == a) & (A.trait == t)]
            if x.empty:
                ax.text(j, y, '–', ha='center', va='center', color=AXIS, fontsize=8); continue
            r = x.iloc[0]; b, se, p, ok = r.beta_alt, r.se_p, r.p, bool(r.qc_pass)
        z = b / se; size = 14 + 22 * min(-np.log10(max(p, 1e-300)), 25) ** 0.8
        col = DIV(norm_(np.clip(z, -8, 8)))
        if a == 'POOL':
            ax.scatter(j, y, s=size, marker='D', color=col, edgecolor=INK, lw=0.9, zorder=3)
        else:
            ax.scatter(j, y, s=size, color=col, edgecolor='white' if ok else INK, lw=0.6, zorder=3, hatch=None if ok else '////')
            if not ok: ax.scatter(j, y, s=size, facecolor='none', edgecolor=INK, lw=0.8, zorder=4); ax.text(j + 0.2, y + 0.33, 'excluded', fontsize=5, color=INK)
        ax.text(j + 0.13 + 0.012 * np.sqrt(size), y, f'{b:+.3f}'.replace('-', '−'), ha='left', va='center', fontsize=5.2, color=INK2)
ax.axhline(-5.55, color=AXIS, lw=0.6); ax.axvline(2.72, color=AXIS, lw=0.8)
ax.set_xticks(range(6)); ax.set_xticklabels([t[1] for t in traits] * 2); ax.xaxis.tick_top(); ax.tick_params(axis='x', length=0, pad=2)
ax.set_yticks([-i for i in range(7)]); ax.set_yticklabels([ANCLAB[a] for a in ancs] + ['pooled (RE)']); ax.tick_params(axis='y', length=0)
for k, (s, lab) in enumerate(snps): ax.text(1 + 3 * k, 1.35, lab.replace('\n', '  '), ha='center', fontsize=7, fontweight='bold', color=INK)
ax.set_xlim(-0.5, 5.9); ax.set_ylim(-6.6, 0.6); clean(ax, False, False)
sm = plt.cm.ScalarMappable(cmap=DIV, norm=norm_); cax = ax.inset_axes([0.02, -0.13, 0.4, 0.035]); cb = fig.colorbar(sm, cax=cax, orientation='horizontal'); cb.outline.set_visible(False)
cb.set_label('z = β / SE, LTL-shortening allele (blue = lower BP)', fontsize=6); cb.ax.tick_params(labelsize=5.8)
for pv, lab in [(1e-2, 'P = 10⁻²'.replace('⁻²', '$^{-2}$')), (1e-8, 'P = 10$^{-8}$'), (1e-20, 'P ≤ 10$^{-20}$')]:
    pass
leg_sizes = [2, 8, 20]
hs = [ax.scatter([], [], s=14 + 22 * v ** 0.8, color=AXIS) for v in leg_sizes]
ax.legend(hs, [f'P = 10$^{{-{v}}}$' for v in leg_sizes], loc='upper left', bbox_to_anchor=(0.5, -0.07), fontsize=5.8, ncol=3, columnspacing=1.5, title='bubble area', title_fontsize=5.8)
ax.set_title('Per-allele associations across ancestries', pad=26); tag(ax, 'a', x=-0.16, y=1.13)

# ---- b  allelic alignment plane
ax = fig.add_subplot(gs[0, 4:])
S = F[(F.code == 'S_mvpEU') & (F.analysis == 'r2<0.001 (primary)') & (F.method == 'IVW (MRE, t)')].iloc[0]
xx = np.linspace(-0.14, 0.02, 10); ax.fill_between(xx, (S.lo) * xx, (S.hi) * xx, color=MUTED, alpha=0.18, lw=0)
ax.plot(xx, S.b * xx, color=INK2, lw=0.9, ls='--'); ax.text(-0.075, -0.045, f'genome-wide MR slope\n({int(S.nsnp)} instruments)', fontsize=5.6, color=INK2, va='top')
for s, mk in [('rs10936599', 'o'), ('rs2736100', 's')]:
    for a in ['EU', 'SA', 'EA', 'AF']:
        l = A[(A.rs == s) & (A.ancestry == a) & (A.trait == 'LTL')]; b = A[(A.rs == s) & (A.ancestry == a) & (A.trait == 'SBP')]
        if l.empty or b.empty: continue
        l, b = l.iloc[0], b.iloc[0]
        ax.errorbar(l.beta_alt, b.beta_alt, xerr=1.96 * l.se_p, yerr=1.96 * b.se_p, fmt='none', ecolor=ANC[a], elinewidth=0.8, alpha=0.8)
        ax.scatter(l.beta_alt, b.beta_alt, marker=mk, s=28, color=ANC[a], edgecolor='white', lw=0.6, zorder=4)
ax.axhline(0, color=AXIS, lw=0.6); ax.axvline(0, color=AXIS, lw=0.6)
ax.set_xlabel('allele effect on LTL (SD)'); ax.set_ylabel('allele effect on SBP (SD)'); ax.set_xlim(-0.15, 0.02); ax.set_ylim(-0.06, 0.03)
hl = [Line2D([], [], marker='o', ls='', color=INK2, label='rs10936599 T'), Line2D([], [], marker='s', ls='', color=INK2, label='rs2736100 A')] + \
     [Line2D([], [], marker='o', ls='', color=ANC[a], label=ANCLAB[a]) for a in ['EU', 'SA', 'EA', 'AF']]
ax.legend(handles=hl, loc='upper center', bbox_to_anchor=(0.5, -0.22), fontsize=5.5, ncol=3, columnspacing=0.6, handletextpad=0.2)
ax.set_title('Allelic alignment: LTL vs SBP'); tag(ax, 'b')

# ---- c  allele frequency strip, 26 populations grouped by super-population
ax = fig.add_subplot(gs[1, :])
sup = {'AFR': ['YRI', 'LWK', 'GWD', 'MSL', 'ESN', 'ASW', 'ACB'], 'AMR': ['MXL', 'PUR', 'CLM', 'PEL'], 'EAS': ['CHB', 'JPT', 'CHS', 'CDX', 'KHV'],
       'EUR': ['CEU', 'TSI', 'FIN', 'GBR', 'IBS'], 'SAS': ['GIH', 'PJL', 'BEB', 'STU', 'ITU']}
supname = {'AFR': 'African', 'AMR': 'Admixed American', 'EAS': 'East Asian', 'EUR': 'European', 'SAS': 'South Asian'}
supcol = {'AFR': VIOLET, 'AMR': YELLOW, 'EAS': AQUA, 'EUR': BLUE, 'SAS': ORANGE}
order = [(s_, p) for s_ in sup for p in sup[s_]]
fq = np.full((2, len(order)), np.nan)
for k, (rs, al) in enumerate([('rs10936599', 'T'), ('rs2736100', 'A')]):
    for i, (s_, p) in enumerate(order):
        a = RF[(RF.rs == rs) & (RF.population == f'1000GENOMES:phase_3:{p}') & (RF.allele == al)]
        if len(a): fq[k, i] = a.freq.iloc[0]
xs = []; x0 = 0.0
for s_ in sup:
    for p in sup[s_]: xs.append(x0); x0 += 1
    x0 += 0.6
xs = np.array(xs); fcm = LinearSegmentedColormap.from_list('f', ['#f4f8fe', '#9ec5f4', '#3987e5', '#184f95'])
for k in range(2):
    for i in range(len(order)):
        v = fq[k, i]; ax.add_patch(plt.Rectangle((xs[i] - 0.46, -k - 0.42), 0.92, 0.84, color=fcm(v) if np.isfinite(v) else MID, lw=0))
        ax.text(xs[i], -k, '' if not np.isfinite(v) else f'{v:.2f}'.lstrip('0'), ha='center', va='center', fontsize=5.3, color='white' if np.nan_to_num(v) > 0.55 else INK)
for i, (s_, p) in enumerate(order): ax.text(xs[i], 0.62, p, ha='center', va='bottom', fontsize=5.3, rotation=90, color=INK2)
for s_ in sup:
    idx = [i for i, (q, _) in enumerate(order) if q == s_]
    ax.plot([xs[idx[0]] - 0.46, xs[idx[-1]] + 0.46], [-1.62, -1.62], color=supcol[s_], lw=3, solid_capstyle='butt')
    ax.text((xs[idx[0]] + xs[idx[-1]]) / 2, -1.85, f'{supname[s_]}\nmean {np.nanmean(fq[0, idx]):.2f} / {np.nanmean(fq[1, idx]):.2f}', ha='center', va='top', fontsize=5.6, color=supcol[s_], fontweight='bold')
ax.set_xlim(xs[0] - 0.8, xs[-1] + 0.8); ax.set_ylim(-2.45, 1.25); ax.set_yticks([0, -1]); ax.set_yticklabels(['rs10936599 T\n(LTL-shortening)', 'rs2736100 A\n(LTL-shortening)'], fontsize=5.8)
ax.set_xticks([]); clean(ax, False, False); ax.tick_params(length=0)
ax.set_title('Allele frequency in 26 populations (1000 Genomes phase 3), grouped by super-population', pad=4); tag(ax, 'c', x=-0.1, y=1.02)
fig.text(0.01, 0.0, 'Source: Knowledge Portal ancestry-specific meta-analyses (accessed 29 Sep 2026); SE derived from P. One internally inconsistent record (rs10936599, Hispanic, SBP) excluded from pooling. '
         'Cross-ancestry pooling is secondary; ancestry-specific estimates are primary. Panel b slope: primary IVW estimate for SBP (MVP European).', fontsize=5.6, color=MUTED, wrap=True)
save(fig, 'Fig4_cross_ancestry')
