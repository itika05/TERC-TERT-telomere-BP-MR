"""Figure 1 (v4): study design. Sample sizes are those recorded in Table 2 of the manuscript; instrument counts are read from the analysis outputs."""
from figstyle import *
import numpy as np, pandas as pd
from matplotlib.patches import FancyBboxPatch, Polygon, PathPatch
from matplotlib.path import Path

LED = pd.read_csv('step3v3_reverse_ledger.csv'); FL = pd.read_csv('step6_v4_forward_variant_ledger.csv')
d = pd.read_csv('kp_ltl_leads_wide.csv')
f1 = len(d); f2 = int((d.L_codd_p < 5e-8).sum()); e = d[(d.L_codd_p < 5e-8) & d.rsid.notna()]
f3 = int((np.minimum(e.maf, 1 - e.maf) >= 0.01).sum()); f3a = int((~FL.mhc).sum() - (FL.chr.astype(str) == 'X').sum()) if 'chr' in FL else f3 - 14
f3b = int(FL.usable.sum()); f4 = int(FL['selected_r2_0.001'].sum())
MV = pd.read_csv('mr_v4_inputs/mvmr_instruments.csv')

# resource, n (for bar), label n, colour, module targets
RES = [('Codd 2021 LTL GWAS (UK Biobank)', 464716, '464,716', BLUE, [0, 1, 2]),
       ('MVP BP, hypertension (European)', 425740, '425,740', BLUE, [0, 1, 2, 3]),
       ('FinnGen R12 hypertension', 500264, '154,630 cases', BLUE, [0, 2]),
       ('MVP African American, Hispanic', 119331, '≤119,331', MAGENTA, [0, 3]),
       ('Biobank Japan SBP, DBP', 136597, '136,597', AQUA, [0, 3]),
       ('Genes & Health SBP, DBP', 18536, '18,536', ORANGE, [0, 3]),
       ('Nakao 2026 LTL (South Asian)', 11277, '11,277', ORANGE, [0]),
       ('BMI (GIANT + UKB); lymphocytes (BCX)', 694649, '≤694,649', YELLOW, [1]),
       ('GTEx v8 eQTL (per tissue)', 670, '73–670', AQUA, [2])]

MOD = [('1  Bidirectional MR (primary)', f'forward: {f4} LTL instruments → SBP, DBP, hypertension\nreverse: BP instruments → LTL; 6 primary tests', 'Fig. 2', BLUE),
       ('2  Pleiotropy and mediation checks', f'phenome screen; published Codd instrument;\nmultivariable MR with BMI, lymphocytes ({len(MV)} variants)', 'Fig. 2, S1', YELLOW),
       ('3  Locus architecture', 'coloc.abf (3 windows × 4 priors); SuSiE;\nMVP and FinnGen outcomes; GTEx eQTLs', 'Fig. 3', AQUA),
       ('4  Cross-ancestry look-up', 'rs10936599, rs2736100 → SBP, DBP,\nhypertension in 6 ancestry groups', 'Fig. 4', MAGENTA)]

fig = plt.figure(figsize=(W, 6.3))
# ---------- a: resource ladder (log sample size)
axa = fig.add_axes([0.25, 0.42, 0.24, 0.52])
y = np.arange(len(RES))[::-1]
for yi, (lab, n, nl, col, _) in zip(y, RES):
    axa.barh(yi, np.log10(n), color=col, alpha=0.85, height=0.62, lw=0)
    axa.text(np.log10(n) + 0.08, yi, nl, va='center', fontsize=5.6, color=INK)
    axa.text(-0.12, yi, lab, va='center', ha='right', fontsize=5.8, color=INK, transform=axa.get_yaxis_transform() if False else axa.transData)
axa.set_xlim(0, 7.2); axa.set_xticks([1, 2, 3, 4, 5, 6]); axa.set_xticklabels(['10', '10²', '10³', '10⁴', '10⁵', '10⁶'])
axa.set_yticks([]); clean(axa, left=False); axa.set_xlabel('participants (log scale)', fontsize=6)
for xg in range(1, 7): axa.axvline(xg, color=GRID, lw=0.4, zorder=0)
axa.set_title('Data resources', pad=6, x=-0.9)
axa.text(-1.02, 1.07, 'a', transform=axa.transAxes, fontsize=10, fontweight='bold')

# ---------- b: modules with bezier connectors
axb = fig.add_axes([0.52, 0.42, 0.48, 0.52]); axb.set_xlim(0, 10); axb.set_ylim(-0.75, len(RES) - 0.25); axb.axis('off')
ymod = np.linspace(len(RES) - 1.1, 0.1, len(MOD))
for (t, body, fr, col), ym in zip(MOD, ymod):
    axb.add_patch(FancyBboxPatch((2.2, ym - 0.7), 6.2, 1.4, boxstyle='round,pad=0.02,rounding_size=0.12', fc='white', ec=col, lw=1.1, zorder=3))
    axb.add_patch(FancyBboxPatch((2.2, ym - 0.7), 0.16, 1.4, boxstyle='square,pad=0', fc=col, ec='none', zorder=4))
    axb.text(2.5, ym + 0.55, t, fontsize=6.4, fontweight='bold', va='top', color=INK, zorder=5)
    axb.text(2.5, ym + 0.12, body, fontsize=5.2, va='top', color=INK2, zorder=5, linespacing=1.15)
    axb.add_patch(FancyBboxPatch((8.65, ym - 0.22), 1.1, 0.44, boxstyle='round,pad=0.02,rounding_size=0.2', fc=col, ec='none', zorder=3))
    axb.text(9.2, ym, fr, ha='center', va='center', fontsize=5.8, color='white', fontweight='bold', zorder=4)
for yi, (lab, n, nl, col, tg) in zip(y, RES):
    for t in tg:
        ym = ymod[t]
        verts = [(0.0, yi), (1.2, yi), (1.0, ym), (2.2, ym)]
        axb.add_patch(PathPatch(Path(verts, [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]), fc='none', ec=col, lw=0.8, alpha=0.55, zorder=1))
    axb.scatter([0.0], [yi], s=8, color=col, zorder=2)
axb.set_title('Analysis modules', loc='left', x=0.22, pad=6)
axb.text(0.0, 1.07, 'b', transform=axb.transAxes, fontsize=10, fontweight='bold')

# ---------- c: instrument funnels
def funnel(ax, x0, steps, col, title):
    mx = steps[0][1]; w = 1.72
    for i, (lab, n) in enumerate(steps):
        h = n / mx
        ax.add_patch(Polygon([(x0 + i * w, 0.5 - h / 2), (x0 + i * w, 0.5 + h / 2), (x0 + i * w + w * 0.78, 0.5 + h / 2), (x0 + i * w + w * 0.78, 0.5 - h / 2)], fc=col, alpha=0.25 + 0.6 * (i + 1) / len(steps), lw=0))
        if i < len(steps) - 1:
            h2 = steps[i + 1][1] / mx
            ax.add_patch(Polygon([(x0 + i * w + w * 0.78, 0.5 - h / 2), (x0 + i * w + w * 0.78, 0.5 + h / 2), (x0 + (i + 1) * w, 0.5 + h2 / 2), (x0 + (i + 1) * w, 0.5 - h2 / 2)], fc=col, alpha=0.12, lw=0))
        ax.text(x0 + i * w + w * 0.39, 0.5, f'{n}', ha='center', va='center', fontsize=6, fontweight='bold', color=INK)
        ax.text(x0 + i * w + w * 0.39, -0.08, lab, ha='center', va='top', fontsize=4.9, color=INK2, linespacing=1.05)
    ax.text(x0, 1.1, title, fontsize=6.3, fontweight='bold', color=col)

axc = fig.add_axes([0.02, 0.08, 0.96, 0.24]); axc.set_xlim(0, 26); axc.set_ylim(-0.6, 1.25); axc.axis('off')
funnel(axc, 0.0, [('Knowledge Portal\nLTL leads', f1), ('Codd\nP<5×10$^{-8}$', f2), ('rsID,\nMAF ≥ 1%', f3), ('autosomal, non-\nMHC, in panel', f3b), ('clumped\nr² < 0.001', f4)], BLUE, 'Forward MR: LTL instruments')
for k, (ex, col) in enumerate([('SBP', ORANGE), ('DBP', YELLOW), ('Hypertension', MAGENTA)]):
    r = LED[LED.exposure == ex].iloc[0]
    funnel(axc, 9.2 + k * 5.6, [('MVP EUR\ngenome-wide', int(r.candidates)), ('in EUR\npanel', int(r.usable)), ('clumped\nr² < 0.001', int(r.selected))], col, f'Reverse MR: {ex}')
fig.text(0.02, 0.345, 'c', fontsize=10, fontweight='bold')
fig.text(0.02, 0.005, 'Primary test family: forward MR on MVP European SBP, DBP, hypertension and reverse MR on Codd LTL (6 tests; Bonferroni α = 0.0083); all other outcomes and ancestries secondary. '
         'None of the outcome samples includes UK Biobank. Forward candidates: 9 chrX and 5 MHC (chr6:25–34 Mb) variants excluded. FinnGen bar: total N (500,264). '
         'Clumping: r² < 0.001 within 10 Mb; LD from 1000 Genomes phase 3 European haplotypes (503 individuals). Reverse candidates: MVP European P < 5×10$^{-8}$, MAF ≥ 1%, rsID, MHC excluded.',
         fontsize=5.3, color=MUTED, wrap=True)
save(fig, 'Fig1_workflow')
