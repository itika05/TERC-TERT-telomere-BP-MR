"""Supplementary Figure S4 (v5): complete coloc.abf prior x window grid (moved from main Figure 3c). PP4 (large) and PP3 (small) for every
locus x trait x window x p12; discordant and wider-window cells are retained."""
from figstyle import *
import numpy as np, pandas as pd
A = pd.read_csv('step6_coloc_v4_abf.csv')
TR = ['SBP', 'DBP', 'HYPERTENSION', 'FG_HYPERTENSION']
TRL = {'SBP': 'SBP (MVP)', 'DBP': 'DBP (MVP)', 'HYPERTENSION': 'Hypertension (MVP)', 'FG_HYPERTENSION': 'Hypertension (FinnGen)'}
fig = plt.figure(figsize=(W, 4.2))
# ---- c  coloc.abf PP4 grid
ax = fig.add_axes([0.2, 0.12, 0.78, 0.74])
rows = [(l, t) for l in ['TERC', 'TERT'] for t in TR]; cols = [(w, p) for w in [50, 100, 250] for p in [1e-6, 5e-6, 1e-5, 5e-5]]
ppcm = LinearSegmentedColormap.from_list('pp', ['#f4f2ec', '#cdc4e8', VIOLET])
M = np.array([[A[(A.locus == l) & (A.trait == t) & (A.window_kb == w) & (A.p12 == p)].PP4.iloc[0] for (w, p) in cols] for (l, t) in rows])
M3 = np.array([[A[(A.locus == l) & (A.trait == t) & (A.window_kb == w) & (A.p12 == p)].PP3.iloc[0] for (w, p) in cols] for (l, t) in rows])
ax.imshow(M, cmap=ppcm, vmin=0, vmax=1, aspect='auto')
for i in range(len(rows)):
    for k in range(len(cols)):
        ax.text(k, i - 0.14, ('0' if M[i, k] < 0.005 else f'{M[i, k]:.2f}'.lstrip('0')), ha='center', va='center', fontsize=6.4, color='white' if M[i, k] > 0.6 else INK, fontweight='bold')
        ax.text(k, i + 0.24, 'H3 ' + ('0' if M3[i, k] < 0.005 else f'{M3[i, k]:.2f}'.lstrip('0')), ha='center', va='center', fontsize=5.0, color='white' if M[i, k] > 0.6 else INK2)
ax.set_yticks(range(len(rows))); ax.set_yticklabels([f'{l}  {TRL[t]}' for l, t in rows], fontsize=5.6)
ax.set_xticks([1.5, 5.5, 9.5]); ax.set_xticklabels(['±50 kb', '±100 kb', '±250 kb'], fontsize=5.8); ax.tick_params(length=0)
for k in [3.5, 7.5]: ax.axvline(k, color='white', lw=1.8)
ax.axhline(3.5, color='white', lw=2.2)
for sp in ax.spines.values(): sp.set_visible(False)
ax.set_title('coloc.abf posterior probability of a shared variant with LTL (PP4, large; PP3 = distinct variants, small)', fontsize=6.3)

ax.text(0.5, -0.17, 'within each window, columns are p12 = 10$^{-6}$, 5×10$^{-6}$, 10$^{-5}$, 5×10$^{-5}$ (p1 = p2 = 10$^{-4}$)', transform=ax.transAxes, ha='center', fontsize=5.2, color=INK2)

save(fig, 'FigS4_coloc_full_grid')
