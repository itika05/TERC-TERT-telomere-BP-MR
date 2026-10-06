"""Figure 3 (v5): (c) simplified to stacked H0-H4 posteriors at p12 = 1e-5 per window with the prior range, plus SuSiE at ±100 and ±250 kb
(full prior grid in Supplementary Figure S4); (e) unit-labelled text rows with SE, no common quantitative axis.
Figure 3 (v4): locus architecture at TERC and TERT using single-cohort data: LTL from the original UK Biobank summary statistics
(Codd 2021), SBP, DBP and hypertension from MVP (European) and hypertension from FinnGen R12; LD from 1000 Genomes EUR haplotypes; coloc.abf PP4
across windows and priors; SuSiE credible sets; GTEx eQTLs; parallel associations of the rs10936599 T allele."""
from figstyle import *
import numpy as np, pandas as pd
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize, BoundaryNorm, ListedColormap
from ld_ref import load, match, hap_matrix

A = pd.read_csv('step6_coloc_v4_abf.csv'); SUS = pd.read_csv('step6_coloc_v4_susie.csv'); L20 = pd.read_csv('step7_susie_L20_TERT_coloc.csv'); WL = pd.read_csv('step6_coloc_v4_window_ledger.csv'); CS = pd.read_csv('step6_coloc_v4_susie_credible_sets.csv')
DG = pd.read_csv('step6_coloc_v4_susie_diagnostics.csv'); EQ = pd.read_csv('gtex_v8_eqtl_lookup.csv'); ST = pd.read_csv('step4_variant_identity_sign_table.csv')
G = pd.read_csv('ensembl_grch37_genes.csv').set_index('gene')
loci = {'TERC': dict(snp='rs10936599', chr=3, genes=['TERC', 'ACTRT3', 'MYNN', 'LRRC34', 'LRRIQ4'], band='3q26.2'),
        'TERT': dict(snp='rs2736100', chr=5, genes=['TERT', 'CLPTM1L', 'SLC6A3'], band='5p15.33')}
bnd = [0, 0.05, 0.2, 0.4, 0.6, 0.8, 1.0001]; ldc = ListedColormap(['#d9d7cf', '#9ec5f4', '#3987e5', AQUA, YELLOW, RED]); ldn = BoundaryNorm(bnd, ldc.N)
TR = ['SBP', 'DBP', 'HYPERTENSION', 'FG_HYPERTENSION']
TRL = {'LTL': 'LTL (UKB)', 'SBP': 'SBP (MVP)', 'DBP': 'DBP (MVP)', 'HYPERTENSION': 'Hypertension (MVP)', 'FG_HYPERTENSION': 'Hypertension (FinnGen)'}
TRS = {'LTL': 'LTL', 'SBP': 'SBP', 'DBP': 'DBP', 'HYPERTENSION': 'Hyp. MVP', 'FG_HYPERTENSION': 'Hyp. FinnGen'}
TRC = {'LTL': VIOLET, 'SBP': BLUE, 'DBP': AQUA, 'HYPERTENSION': ORANGE, 'FG_HYPERTENSION': MAGENTA}

fig = plt.figure(figsize=(W, 10.8))
outer = fig.add_gridspec(3, 1, height_ratios=[1.45, 0.62, 0.95], hspace=0.3, top=0.92, bottom=0.085)
top = outer[0].subgridspec(1, 2, wspace=0.4)
for j, (loc, L) in enumerate(loci.items()):
    g = top[j].subgridspec(4, 1, height_ratios=[1.1, 1.1, 0.75, 0.55], hspace=0.07)
    w = WL[(WL.locus == loc) & (WL.window_kb == 250)].iloc[0]; cpos = w.centre_pos; lo, hi = cpos - 250000, cpos + 250000
    d = {t: pd.read_csv(f'coloc_v4_inputs/{loc}_250kb_{t}.csv') for t in TR}
    x = d['SBP'][d['SBP'].in_panel.astype(str) == 'True'].reset_index(drop=True)
    ref = load(f'ldref/eur_region_{loc}.json'); m = match(x.varId, ref); H = hap_matrix(m, ref)
    ci = int(np.where(x.rsid.values == w.centre)[0][0]); Hc = H - H.mean(1, keepdims=True)
    x['r2'] = ((Hc @ Hc[ci]) / (np.sqrt((Hc ** 2).sum(1)) * np.sqrt((Hc[ci] ** 2).sum()))) ** 2
    axes = [fig.add_subplot(g[k]) for k in range(4)]
    for k, (col, lab) in enumerate([('ltl', 'LTL (UKB) |z|'), ('out', 'SBP (MVP) |z|')]):
        ax = axes[k]; y = (x[f'b_{col}'] / x[f'se_{col}']).abs(); xx = x.assign(az=y).sort_values('r2')
        yv = xx.az if k == 0 else -xx.az
        ax.scatter(xx.pos / 1e6, yv, c=xx.r2, cmap=ldc, norm=ldn, s=np.where(xx.r2 >= 0.8, 9, 5), lw=0.15, edgecolor='white', zorder=2)
        for rs, mk, sz, c_ in [(w.centre, 'D', 34, VIOLET), (L['snp'], '*', 60, INK)]:
            y1 = xx[xx.rsid == rs]
            if len(y1): ax.scatter(y1.pos / 1e6, y1.az if k == 0 else -y1.az, marker=mk, s=sz, color=c_, edgecolor='white', lw=0.6, zorder=5)
        ax.set_xlim(lo / 1e6, hi / 1e6); ax.tick_params(labelbottom=False, bottom=False); clean(ax, bottom=False)
        ax.set_ylabel(lab, fontsize=6.0); ax.axhline(0, color=INK2, lw=0.5)
        if k == 1:
            ax.set_ylim(top=0); tk = [v for v in ax.get_yticks() if v < 0]; ax.set_yticks(tk); ax.set_yticklabels([f'{abs(v):.0f}' for v in tk]); ax.set_ylim(top=0)
        else: ax.set_ylim(bottom=0)
    axes[0].set_title(f'{loc} ({L["band"]}): LTL above, SBP below', loc='left', fontsize=6.8)
    ax = axes[2]
    lanes = ['LTL'] + TR
    for q, t in enumerate(lanes):
        if t == 'LTL': cs = CS[(CS.locus == loc) & (CS.window_kb == 250) & (CS.pair == 'SBP') & (CS.trait == 'LTL')]; dd = d['SBP']
        else: cs = CS[(CS.locus == loc) & (CS.window_kb == 250) & (CS.pair == t) & (CS.trait == t)]; dd = d[t]
        for c in cs.itertuples():
            mem = set(c.members.split(';')); pp = dd[dd.rsid.isin(mem)].pos.values / 1e6
            if len(pp): ax.plot([pp.min(), max(pp.max(), pp.min() + 0.0015)], [-q, -q], color=TRC[t], lw=3, solid_capstyle='butt', alpha=0.8)
            lp = dd[dd.rsid == c.lead_rsid].pos.values
            if len(lp): ax.scatter(lp / 1e6, [-q], s=10, color='white', edgecolor=TRC[t], lw=0.8, zorder=4)
        if not len(cs): ax.text(lo / 1e6 + 0.01, -q, 'no credible set', fontsize=4.8, va='center', color=MUTED)
        dg = DG[(DG.locus == loc) & (DG.window_kb == 250) & (DG.pair == ('SBP' if t == 'LTL' else t)) & (DG.trait == t)]
        s_ = dg.s_estimate.iloc[0]
        ax.text(1.005, -q, f's={s_:.2f}', transform=ax.get_yaxis_transform(), fontsize=5, va='center', color=RED if s_ > 0.1 else INK2)
    ax.axvline(cpos / 1e6, color=VIOLET, lw=0.5, ls=':')
    ax.set_yticks([-q for q in range(len(lanes))]); ax.set_yticklabels([TRS[t] for t in lanes], fontsize=5.3); ax.set_ylim(-len(lanes) + 0.4, 0.6)
    ax.set_xlim(lo / 1e6, hi / 1e6); ax.tick_params(labelbottom=False, bottom=False, length=0); clean(ax, False, False)
    ax.text(1.0, 1.03, 'SuSiE 95% credible sets (1000G EUR LD; s = LD-mismatch)', transform=ax.transAxes, fontsize=5.0, color=INK2, ha='right')
    ax = axes[3]
    for q, gn in enumerate(L['genes']):
        s0, e0, strand = G.loc[gn, 'start'], G.loc[gn, 'end'], G.loc[gn, 'strand']; yy = 0.9 - q * 0.36
        ax.add_patch(Rectangle((s0 / 1e6, yy - 0.1), max((e0 - s0) / 1e6, 0.002), 0.2, color=INK, lw=0))
        ax.text(e0 / 1e6 + 0.004, yy, gn + (' ←' if strand < 0 else ' →'), fontsize=5.6, style='italic', va='center', color=INK)
    ax.set_ylim(0.9 - len(L['genes']) * 0.36, 1.1); ax.set_yticks([]); ax.set_xlim(lo / 1e6, hi / 1e6); clean(ax, left=False)
    ax.set_xlabel(f'chr{L["chr"]} (Mb, GRCh37)')
    axes[0].text(-0.22, 1.08, 'ab'[j], transform=axes[0].transAxes, fontsize=10, fontweight='bold')
cax = fig.add_axes([0.33, 0.985, 0.34, 0.006]); cb = fig.colorbar(plt.cm.ScalarMappable(cmap=ldc, norm=ldn), cax=cax, orientation='horizontal', ticks=[0.05, 0.2, 0.4, 0.6, 0.8])
cb.outline.set_visible(False); cb.ax.tick_params(labelsize=5.5, length=1.5); cax.set_title('r² with the lead LTL variant (1000 Genomes European haplotypes)', fontsize=5.7, color=INK2, fontweight='normal', pad=2)
fig.legend(handles=[Line2D([], [], marker='D', ls='', color=VIOLET, label='lead LTL variant (TERC rs2293607; TERT rs7705526)'), Line2D([], [], marker='*', ls='', color=INK, ms=8, label='candidate SNP (rs10936599; rs2736100)')],
           loc='upper center', bbox_to_anchor=(0.5, 0.972), ncol=2, fontsize=5.8)

# ---- c  key comparisons: stacked posteriors (p12 = 1e-5) with prior range; SuSiE columns
ax = fig.add_subplot(outer[1])
rows = [(l, t) for l in ['TERC', 'TERT'] for t in TR]
HC = ['#e1e0d9', '#cfcdc4', '#bdbbb1', ORANGE, VIOLET]
LEAD = {'TERC': '3:169482335:T:C', 'TERT': '5:1285974:C:A'}
xw = [0, 1.15, 2.3]; bw = 1.0
for i, (l, t) in enumerate(rows):
    y = i
    for k, w in enumerate([50, 100, 250]):
        r = A[(A.locus == l) & (A.trait == t) & (A.window_kb == w) & np.isclose(A.p12, 1e-5)].iloc[0]
        left = xw[k]
        for h_, col in zip(['PP0', 'PP1', 'PP2', 'PP3', 'PP4'], HC):
            v = r[h_]; ax.add_patch(Rectangle((left, y - 0.32), v * bw, 0.64, color=col, lw=0)); left += v * bw
        allp = A[(A.locus == l) & (A.trait == t) & (A.window_kb == w)].PP4
        x4a, x4b = xw[k] + bw - allp.max() * bw, xw[k] + bw - allp.min() * bw
        ax.plot([x4a, x4b], [y + 0.42, y + 0.42], color=INK, lw=0.9); ax.plot([x4a] * 2, [y + 0.36, y + 0.48], color=INK, lw=0.6); ax.plot([x4b] * 2, [y + 0.36, y + 0.48], color=INK, lw=0.6)
        v4 = r.PP4; ax.text(xw[k] + bw - 0.02, y, ('<0.01' if v4 < 0.01 else f'{v4:.2f}'), ha='right', va='center', fontsize=5.0, color='white' if v4 > 0.25 else INK, fontweight='bold')
        v3 = r.PP3
        if v3 > 0.25: ax.text(xw[k] + 0.02 + (1 - v3 - v4) * bw, y, f'H3 {v3:.2f}', ha='left', va='center', fontsize=4.4, color='white')
    for k, w in enumerate([100, 250]):
        xs = 3.75 + k * 0.85
        x_ = SUS[(SUS.locus == l) & (SUS.trait == t) & (SUS.window_kb == w) & (SUS.hit1 == LEAD[l]) & np.isclose(SUS.p12, 1e-5)]
        if x_.empty or x_.PP4.isna().all(): txt, c_ = 'no CS', MUTED
        else:
            v = x_.PP4.max(); txt = '<0.01' if v < 0.01 else f'{v:.2f}'; c_ = VIOLET if v > 0.8 else INK
        ax.text(xs, y - 0.05, txt, ha='center', va='center', fontsize=5.6, color=c_, fontweight='bold')
        if l == 'TERT':
            z = L20[(L20.trait == t) & (L20.window_kb == w)]
            if len(z): ax.text(xs, y + 0.3, 'L=20: ' + ('no CS' if z.maxPP4.isna().all() else f'{z.maxPP4.iloc[0]:.2f}'), ha='center', va='center', fontsize=4.2, color=RED)
ax.set_xlim(-0.05, 5.2); ax.set_ylim(len(rows) - 0.4, -1.0)
for k, lab in enumerate(['±50 kb', '±100 kb', '±250 kb']): ax.text(xw[k] + 0.5, -0.72, f'coloc.abf {lab}', ha='center', fontsize=5.6, fontweight='bold', color=INK)
for k, lab in enumerate(['±100 kb', '±250 kb']): ax.text(3.75 + k * 0.85, -0.72, f'SuSiE {lab}', ha='center', fontsize=5.6, fontweight='bold', color=INK)
ax.axvline(3.4, color=AXIS, lw=0.6); ax.axhline(3.5, color=AXIS, lw=0.6)
ax.set_yticks(range(len(rows))); ax.set_yticklabels([f'{l}  {TRL[t]}' for l, t in rows], fontsize=5.6); ax.set_xticks([]); clean(ax, False, False); ax.tick_params(length=0)
ax.set_title('Colocalization with LTL: posterior of H0–H4 at p12 = 10$^{-5}$ (bar), PP4 range across 4 priors (bracket); coloc.susie PP4 for the lead LTL credible set', fontsize=6.0)
lg = [Rectangle((0, 0), 1, 1, color=c, label=l_) for c, l_ in [(HC[1], 'H0–H2 (no/one-trait signal)'), (ORANGE, 'H3 distinct variants'), (VIOLET, 'H4 shared variant')]]
ax.legend(handles=lg, loc='upper left', bbox_to_anchor=(0.0, -0.02), ncol=3, fontsize=5.0)
ax.text(1.0, -0.06, 'TERT SuSiE: LTL fit saturated at L = 10; red = L = 20 refit (not stable; see text)', transform=ax.transAxes, ha='right', fontsize=4.8, color=RED)
tag(ax, 'c', x=-0.2, y=1.06)

# ---- d  eQTL bubbles ; e  parallel associations
bot = outer[2].subgridspec(1, 2, width_ratios=[1.0, 1.25], wspace=0.45)
ax = fig.add_subplot(bot[0])
tis = ['Artery_Aorta', 'Artery_Tibial', 'Artery_Coronary', 'Heart_Left_Ventricle', 'Heart_Atrial_Appendage', 'Kidney_Cortex', 'Adrenal_Gland', 'Whole_Blood']
tl = ['Aorta', 'Tibial', 'Coron.', 'LV', 'Atrium', 'Kidney', 'Adrenal', 'Blood']
rowsE = [('rs10936599', 'LRRC34'), ('rs10936599', 'MYNN'), ('rs10936599', 'ACTRT3'), ('rs10936599', 'TERC'), ('rs3821383', 'LRRC34'), ('rs3821383', 'MYNN'), ('rs2736100', 'TERT'), ('rs2736100', 'CLPTM1L'), ('rs7726159', 'CLPTM1L')]
nn = Normalize(-0.45, 0.45)
for i, (rs, gn) in enumerate(rowsE):
    for k, t in enumerate(tis):
        e = EQ[(EQ.rs == rs) & (EQ.gene == gn) & (EQ.tissue == t)]
        if e.empty or pd.isna(e.nes.iloc[0]): ax.text(k, -i, '×', ha='center', va='center', fontsize=6, color=AXIS); continue
        e = e.iloc[0]; sig = pd.notna(e.p_thresh_gene_level) and e.p < float(e.p_thresh_gene_level)
        ax.scatter(k, -i, s=6 + 14 * min(-np.log10(e.p), 10), color=DIV(nn(e.nes)), edgecolor=INK if sig else 'white', lw=0.9 if sig else 0.3)
ax.set_xticks(range(len(tis))); ax.set_xticklabels(tl, fontsize=5.8); ax.xaxis.tick_top(); ax.tick_params(length=0)
ax.set_yticks([-i for i in range(len(rowsE))]); ax.set_yticklabels([f'{rs}  {gn}' for rs, gn in rowsE], fontsize=5.6)
ax.set_xlim(-0.6, len(tis) - 0.4); ax.set_ylim(-len(rowsE) + 0.4, 0.6); clean(ax, False, False)
for yy in [-3.5, -5.5]: ax.axhline(yy, color=GRID, lw=0.8)
cax = ax.inset_axes([0.25, -0.1, 0.5, 0.035]); cb = fig.colorbar(plt.cm.ScalarMappable(cmap=DIV, norm=nn), cax=cax, orientation='horizontal'); cb.outline.set_visible(False); cb.ax.tick_params(labelsize=5); cb.set_label('NES (effect allele)', fontsize=5.4)
ax.set_title('GTEx v8 eQTLs (ring: gene-level significant; ×: no data)', pad=14, fontsize=6.3); tag(ax, 'd', x=-0.42, y=1.12)

ax = fig.add_subplot(bot[1]); tr = ST.set_index('rsid').loc['rs10936599']
from scipy.stats import norm as _n
def rec(t):
    x = pd.read_csv(f'coloc_v4_inputs/TERC_250kb_{t}.csv'); return x[x.rsid == 'rs10936599'].iloc[0]
r_s, r_d, r_h, r_f = rec('SBP'), rec('DBP'), rec('HYPERTENSION'), rec('FG_HYPERTENSION')
pz = lambda b, se: 2 * _n.sf(abs(b / se))
items = [('LTL', 'UK Biobank', r_s.b_ltl, r_s.se_ltl, 'SD'), ('SBP', 'MVP EUR', r_s.b_out, r_s.se_out, 'SD'), ('DBP', 'MVP EUR', r_d.b_out, r_d.se_out, 'SD'),
         ('Hypertension', 'MVP EUR', r_h.b_out, r_h.se_out, 'log OR'), ('Hypertension', 'FinnGen R12', r_f.b_out, r_f.se_out, 'log OR')]
q0 = EQ[(EQ.rs == 'rs10936599')]
for gn, t, lab in [('LRRC34', 'Artery_Aorta', 'aorta'), ('LRRC34', 'Heart_Left_Ventricle', 'left ventricle'), ('MYNN', 'Artery_Aorta', 'aorta')]:
    e = q0[(q0.gene == gn) & (q0.tissue == t)].iloc[0]; items.append((f'{gn} expr.', f'GTEx {lab}'.replace('left ventricle', 'LV'), e.nes, float(e.se), 'NES'))
hdr = ['trait', 'source', 'β (SE)', 'unit', 'P']; xs = [0.0, 0.27, 0.51, 0.77, 1.03]
for x_, h_ in zip(xs, hdr): ax.text(x_, 0.75, h_, fontsize=5.4, fontweight='bold', color=INK, ha='left' if x_ < 0.95 else 'right')
for i, (lab, src, b, se, unit) in enumerate(items):
    y = -i
    if i % 2 == 0: ax.add_patch(Rectangle((-0.02, y - 0.45), 1.04, 0.9, color=MID, lw=0))
    ax.text(xs[0], y, ('▼ ' if b < 0 else '▲ ') + lab, fontsize=5.5, va='center', color=BLUE if b < 0 else ORANGE)
    ax.text(xs[1], y, src, fontsize=5.3, va='center', color=INK2)
    ax.text(xs[2], y, f'{b:+.4f} ({se:.4f})'.replace('-', '−'), fontsize=5.5, va='center', color=INK)
    ax.text(xs[3], y, unit, fontsize=5.3, va='center', color=INK2)
    p_ = pz(b, se); ax.text(xs[4], y, '<1e-300' if p_ == 0 else f'{p_:.1e}', fontsize=5.3, va='center', ha='right', color=INK2)
ax.set_xlim(-0.03, 1.03); ax.set_ylim(-len(items) + 0.3, 1.1); ax.axis('off')
ax.set_title('rs10936599-T: parallel associations (each in its own unit)', fontsize=6.3)
ax.text(0.5, -len(items) + 0.05, 'Separate datasets and scales; values are not comparable across rows.\nNo pathway between them is implied. P from β/SE (Wald).', fontsize=5, color=INK2, va='top', ha='center')
tag(ax, 'e', x=-0.12, y=1.06)
fig.text(0.01, 0.004, 'LTL: Codd et al. 2021 UK Biobank summary statistics (β, SE as published). SBP, DBP, hypertension: MVP European (Verma et al. 2024; SE = |β|/z(P)). FinnGen R12 I9_HYPTENS (154,630 cases, 345,634 controls). '
         'The MVP and FinnGen blood-pressure data do not include UK Biobank. coloc.abf (R coloc 6.0.3) assumes one causal variant per trait; SuSiE (susieR 0.12.35) uses an out-of-sample 503-person reference (s > 0.1 in red).',
         fontsize=5.3, color=MUTED, wrap=True)
save(fig, 'Fig3_locus_architecture')
