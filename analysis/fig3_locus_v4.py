"""Figure 3 (v4): locus architecture at TERC and TERT using single-cohort data: LTL from the original UK Biobank summary statistics
(Codd 2021), SBP, DBP and hypertension from MVP (European) and hypertension from FinnGen R12; LD from 1000 Genomes EUR haplotypes; coloc.abf PP4
across windows and priors; SuSiE credible sets; GTEx eQTLs; parallel associations of the rs10936599 T allele."""
from figstyle import *
import numpy as np, pandas as pd
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize, BoundaryNorm, ListedColormap
from ld_ref import load, match, hap_matrix

A = pd.read_csv('step6_coloc_v4_abf.csv'); WL = pd.read_csv('step6_coloc_v4_window_ledger.csv'); CS = pd.read_csv('step6_coloc_v4_susie_credible_sets.csv')
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

# ---- c  coloc.abf PP4 grid
ax = fig.add_subplot(outer[1])
rows = [(l, t) for l in ['TERC', 'TERT'] for t in TR]; cols = [(w, p) for w in [50, 100, 250] for p in [1e-6, 5e-6, 1e-5, 5e-5]]
ppcm = LinearSegmentedColormap.from_list('pp', ['#f4f2ec', '#cdc4e8', VIOLET])
M = np.array([[A[(A.locus == l) & (A.trait == t) & (A.window_kb == w) & (A.p12 == p)].PP4.iloc[0] for (w, p) in cols] for (l, t) in rows])
M3 = np.array([[A[(A.locus == l) & (A.trait == t) & (A.window_kb == w) & (A.p12 == p)].PP3.iloc[0] for (w, p) in cols] for (l, t) in rows])
ax.imshow(M, cmap=ppcm, vmin=0, vmax=1, aspect='auto')
for i in range(len(rows)):
    for k in range(len(cols)):
        ax.text(k, i - 0.14, ('0' if M[i, k] < 0.005 else f'{M[i, k]:.2f}'.lstrip('0')), ha='center', va='center', fontsize=5.2, color='white' if M[i, k] > 0.6 else INK, fontweight='bold')
        ax.text(k, i + 0.24, 'H3 ' + ('0' if M3[i, k] < 0.005 else f'{M3[i, k]:.2f}'.lstrip('0')), ha='center', va='center', fontsize=3.9, color='white' if M[i, k] > 0.6 else INK2)
ax.set_yticks(range(len(rows))); ax.set_yticklabels([f'{l}  {TRL[t]}' for l, t in rows], fontsize=5.6)
ax.set_xticks([1.5, 5.5, 9.5]); ax.set_xticklabels(['±50 kb', '±100 kb', '±250 kb'], fontsize=5.8); ax.tick_params(length=0)
for k in [3.5, 7.5]: ax.axvline(k, color='white', lw=1.8)
ax.axhline(3.5, color='white', lw=2.2)
for sp in ax.spines.values(): sp.set_visible(False)
ax.set_title('coloc.abf posterior probability of a shared variant with LTL (PP4, large; PP3 = distinct variants, small)', fontsize=6.3)
tag(ax, 'c', x=-0.2, y=1.06)
ax.text(0.5, -0.17, 'within each window, columns are p12 = 10$^{-6}$, 5×10$^{-6}$, 10$^{-5}$, 5×10$^{-5}$ (p1 = p2 = 10$^{-4}$)', transform=ax.transAxes, ha='center', fontsize=5.2, color=INK2)

# ---- d  eQTL bubbles ; e  parallel associations
bot = outer[2].subgridspec(1, 2, width_ratios=[1.15, 1], wspace=0.7)
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
def rec(t):
    x = pd.read_csv(f'coloc_v4_inputs/TERC_250kb_{t}.csv'); return x[x.rsid == 'rs10936599'].iloc[0]
r_s, r_d, r_h, r_f = rec('SBP'), rec('DBP'), rec('HYPERTENSION'), rec('FG_HYPERTENSION')
items = [('LTL (UK Biobank)', r_s.b_ltl, 'SD', True), ('SBP (MVP)', r_s.b_out, 'SD', r_s.p_out < 5e-8), ('DBP (MVP)', r_d.b_out, 'SD', r_d.p_out < 5e-8),
         ('Hypertension (MVP)', r_h.b_out, 'log OR', r_h.p_out < 5e-8), ('Hypertension (FinnGen)', r_f.b_out, 'log OR', r_f.p_out < 5e-8), ('SBP (South Asian)', tr.SBP_SA_beta, 'SD', False)]
pv = [None, r_s.p_out, r_d.p_out, r_h.p_out, r_f.p_out, tr.SBP_SA_p]
q0 = EQ[(EQ.rs == 'rs10936599')]
for gn, t, lab in [('LRRC34', 'Artery_Aorta', 'aorta'), ('LRRC34', 'Heart_Left_Ventricle', 'LV'), ('MYNN', 'Artery_Aorta', 'aorta')]:
    e = q0[(q0.gene == gn) & (q0.tissue == t)].iloc[0]; items.append((f'{gn} expression, {lab}', e.nes, 'NES', bool(e.p <= e.p_thresh_gene_level))); pv.append(e.p)
for i, (lab, v, unit, ok) in enumerate(items):
    y = -i; col = BLUE if v < 0 else ORANGE
    ax.add_patch(Rectangle((0, y - 0.34), 1, 0.68, color=col, alpha=0.9 if ok else 0.25, lw=0))
    ax.text(0.5, y, ('↓ ' if v < 0 else '↑ ') + f'{v:+.4f} {unit}'.replace('-', '−'), ha='center', va='center', fontsize=5.6, color='white' if ok else INK, fontweight='bold')
    ax.text(-0.05, y, lab, ha='right', va='center', fontsize=5.6, color=INK)
    p_ = pv[i]; ptxt = 'P < 10$^{-300}$' if p_ is None else f'P = {p_:.1e}'
    ax.text(1.05, y, ptxt, ha='left', va='center', fontsize=5.1, color=INK2 if ok else MUTED)
ax.set_xlim(0, 1); ax.set_ylim(-len(items) + 0.4, 0.8); ax.axis('off')
ax.set_title('rs10936599-T: parallel associations', fontsize=6.3)
ax.text(0.5, -len(items) + 0.1, 'Separate datasets; each cell is an association of the same allele.\nNo pathway between them is implied. Solid: P < 5×10$^{-8}$ (GWAS) or\ngene-level significant (GTEx).', fontsize=5, color=INK2, va='top', ha='center')
tag(ax, 'e', x=-0.5, y=1.06)
fig.text(0.01, 0.004, 'LTL: Codd et al. 2021 UK Biobank summary statistics (β, SE as published). SBP, DBP, hypertension: MVP European (Verma et al. 2024; SE = |β|/z(P)). FinnGen R12 I9_HYPTENS (154,630 cases, 345,634 controls). '
         'None of the blood-pressure samples includes UK Biobank. coloc.abf (R coloc 6.0.3) assumes one causal variant per trait; SuSiE (susieR 0.12.35) uses an out-of-sample 503-person reference (s > 0.1 in red).',
         fontsize=5.3, color=MUTED, wrap=True)
save(fig, 'Fig3_locus_architecture')
