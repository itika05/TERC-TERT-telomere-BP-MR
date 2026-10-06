"""Supplementary Figure S2: (a) phenome screen of the 116 forward instruments; (b) colocalization with multi-cohort Knowledge Portal
European meta-analyses (two SE sources) compared with single-cohort data (MVP, FinnGen) and source LTL statistics."""
from figstyle import *
import numpy as np, pandas as pd
from matplotlib.colors import ListedColormap, BoundaryNorm

S = pd.read_csv('step6_phenome_screen_summary.csv'); L = pd.read_csv('step6_phenome_screen_long.csv')
A3 = pd.read_csv('step4v3_coloc_abf_R.csv'); A4 = pd.read_csv('step6_coloc_v4_abf.csv')
cats = ['Adiposity', 'Lipids', 'Glycaemia / T2D', 'Kidney function', 'Leukocyte counts', 'Clonal haematopoiesis / mCA', 'Blood pressure (not excluded)']
catl = ['Adiposity', 'Lipids', 'Glycemia / T2D', 'Kidney', 'Leukocyte counts', 'CHIP / mCA', 'Blood pressure*']
fig = plt.figure(figsize=(W, 8.2)); gs = fig.add_gridspec(2, 1, height_ratios=[1.35, 1], hspace=0.6)

ax = fig.add_subplot(gs[0])
S = S.copy(); S['n'] = S[[c for c in cats if c in S]].notna().sum(axis=1)
S = S.sort_values(['any_excluding_category', 'n'], ascending=[True, False]).reset_index(drop=True)
M = np.zeros((len(cats), len(S)))
for j, c in enumerate(cats):
    if c in S: M[j] = np.where(S[c].notna(), np.clip(-np.log10(S[c].fillna(1).clip(lower=1e-300)), 0, 60), 0)
cm = LinearSegmentedColormap.from_list('ph', ['#f4f2ec', '#f6b79b', ORANGE, '#a8401a'])
ax.imshow(np.where(M > 0, M, np.nan), cmap=cm, vmin=7, vmax=60, aspect='auto', interpolation='nearest')
ax.set_facecolor('#f7f6f2')
ax.set_yticks(range(len(cats))); ax.set_yticklabels(catl, fontsize=5.8)
ax.set_xticks(range(len(S))); ax.set_xticklabels([f'{g if isinstance(g, str) else r}' for g, r in zip(S.gene, S.rsid)], rotation=90, fontsize=3.6)
ax.tick_params(length=0); nclean = int((~S.any_excluding_category).sum())
ax.axvline(nclean - 0.5, color=INK, lw=0.8)
ax.text(nclean / 2, -0.75, f'retained in phenome-screened set ({nclean})', ha='center', fontsize=5.6, color=INK)
ax.text(nclean + (len(S) - nclean) / 2, -0.75, f'associated (P < 5×10$^{{-8}}$) with ≥1 screened trait ({len(S) - nclean})', ha='center', fontsize=5.6, color=INK)
for sp in ax.spines.values(): sp.set_visible(False)
cax = ax.inset_axes([0.78, -0.2, 0.2, 0.03]); cb = fig.colorbar(plt.cm.ScalarMappable(cmap=cm, norm=plt.Normalize(7, 60)), cax=cax, orientation='horizontal')
cb.outline.set_visible(False); cb.ax.tick_params(labelsize=5); cb.set_label('−log$_{10}$ P (strongest record)', fontsize=5.3)
ax.set_title('Phenome-based instrument screen of the 116 LTL instruments (Knowledge Portal datasets queried 29 Sep 2026)', pad=24); tag(ax, 'a', x=-0.1, y=1.07)
ax.text(0.0, -0.22, '* Blood-pressure associations are shown but were not used for exclusion (they are the outcome pathway).', transform=ax.transAxes, fontsize=5.2, color=INK2)

ax = fig.add_subplot(gs[1])
rows = [(l, t) for l in ['TERC', 'TERT'] for t in ['SBP', 'DBP', 'HYPERTENSION']]
cols = [('Portal meta-analysis,\nSE from P', 'P-derived unless floored (primary)'), ('Portal meta-analysis,\nreported SE', 'portal-reported'),
        ('MVP (primary)', 'v4'), ('FinnGen R12', 'fg')]
Mx = np.full((len(rows), len(cols) * 3), np.nan)
for i, (l, t) in enumerate(rows):
    for j, (_, src) in enumerate(cols):
        for k, w in enumerate([50, 100, 250]):
            if src in ('v4', 'fg'):
                tt = t if src == 'v4' else ('FG_HYPERTENSION' if t == 'HYPERTENSION' else None)
                if tt is None: continue
                x = A4[(A4.locus == l) & (A4.trait == tt) & (A4.window_kb == w) & (A4.p12 == 1e-5)]
            else:
                x = A3[(A3.locus == l) & (A3.trait == t) & (A3.window_kb == w) & (A3.p12 == 1e-5) & (A3.se_source == src)]
            if len(x): Mx[i, j * 3 + k] = x.PP4.iloc[0]
ppcm = LinearSegmentedColormap.from_list('pp', ['#f4f2ec', '#cdc4e8', VIOLET])
ax.imshow(Mx, cmap=ppcm, vmin=0, vmax=1, aspect='auto')
for i in range(Mx.shape[0]):
    for k in range(Mx.shape[1]):
        v = Mx[i, k]; ax.text(k, i, '–' if np.isnan(v) else ('0' if v < 0.005 else f'{v:.2f}'.lstrip('0')), ha='center', va='center', fontsize=5.4, color='white' if np.nan_to_num(v) > 0.6 else INK)
ax.set_yticks(range(len(rows))); ax.set_yticklabels([f'{l}  {t.title() if t == "HYPERTENSION" else t}' for l, t in rows], fontsize=5.8)
ax.set_xticks(range(Mx.shape[1])); ax.set_xticklabels(['±50', '±100', '±250'] * len(cols), fontsize=5.3); ax.tick_params(length=0)
for j, (lab, _) in enumerate(cols): ax.text(j * 3 + 1, -0.6, lab, va='bottom', ha='center', fontsize=5.6, fontweight='bold', color=INK)
for k in [2.5, 5.5, 8.5]: ax.axvline(k, color='white', lw=2.5)
ax.axhline(2.5, color='white', lw=2.5)
for sp in ax.spines.values(): sp.set_visible(False)
ax.set_title('Colocalization PP4 by data source (p12 = 10$^{-5}$; window in kb)', pad=26); tag(ax, 'b', x=-0.1, y=1.12)
ax.text(0.0, -0.2, 'Left: Knowledge Portal European meta-analyses for all traits (several contributing studies include UK Biobank). Right: LTL from the original\n'
        'UK Biobank statistics with single-cohort blood-pressure data that exclude UK Biobank (MVP European; FinnGen R12 hypertension). FinnGen has no SBP or DBP data.',
        transform=ax.transAxes, fontsize=5.3, color=INK2, va='top')
save(fig, 'FigS2_phenome_and_coloc_sources')
