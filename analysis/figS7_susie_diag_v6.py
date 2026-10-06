"""Supplementary Figure S7 (v6): multi-signal (SuSiE/coloc.susie) diagnostics at TERC and TERT. For each locus, BP trait and window: number of LTL credible
sets (L = 10), credible sets led by |z| < 2, LD-mismatch s for LTL, and coloc.susie PP4 for the pair involving the LTL credible set containing the regional
lead (lead pair) and the maximum over all pairs, with haplotype LD and with genotype-dosage LD (same 503 individuals); TERT also with L = 20 (haplotype LD).
Verdict: 'stable' only if the LTL fit is not saturated in either LD version and the lead-pair PP4 is on the same side of 0.5 with both LD definitions."""
from figstyle import *
import numpy as np, pandas as pd
D = pd.read_csv('step8_susie_ld_sensitivity_diagnostics.csv'); S = pd.read_csv('step8_susie_ld_sensitivity_coloc.csv'); L20 = pd.read_csv('step7_susie_L20_TERT_coloc.csv'); D20 = pd.read_csv('step7_susie_L20_TERT_diagnostics.csv')
LEAD = {'TERC': '3:169482335:T:C', 'TERT': '5:1285974:C:A'}; TRL = {'SBP': 'SBP', 'DBP': 'DBP', 'HYPERTENSION': 'Hyp. MVP', 'FG_HYPERTENSION': 'Hyp. FinnGen'}
rows = []
for loc in ['TERC', 'TERT']:
    for t in TRL:
        for w in [100, 250]:
            rec = dict(locus=loc, trait=t, window_kb=w)
            for ld in ['haplotype', 'dosage']:
                d = D[(D.ld == ld) & (D.locus == loc) & (D.window_kb == w) & (D.pair == t) & (D.trait == 'LTL')].iloc[0]
                s = S[(S.ld == ld) & (S.locus == loc) & (S.window_kb == w) & (S.trait == t)]
                rec[f'{ld}_ncs'] = d.n_credible_sets; rec[f'{ld}_lowz'] = d.cs_led_by_abs_z_lt_2; rec[f'{ld}_s'] = d.s_estimate
                rec[f'{ld}_lead'] = s[s.hit1 == LEAD[loc]].PP4.max() if s.PP4.notna().any() else np.nan; rec[f'{ld}_max'] = s.PP4.max()
            if loc == 'TERT':
                rec['L20_ncs'] = D20[(D20.window_kb == w) & (D20.pair == t) & (D20.trait == 'LTL')].n_credible_sets.iloc[0]; rec['L20_max'] = L20[(L20.window_kb == w) & (L20.trait == t)].maxPP4.iloc[0]
            h, dd = rec['haplotype_lead'], rec['dosage_lead']
            stable = (rec['haplotype_ncs'] < 10) and (rec['dosage_ncs'] < 10) and not (np.isnan(h) or np.isnan(dd)) and ((h > 0.5) == (dd > 0.5)) and rec['dosage_lowz'] == 0
            rec['verdict'] = 'stable' if stable else ('no BP credible set' if np.isnan(h) and np.isnan(dd) else 'unstable')
            rows.append(rec)
T = pd.DataFrame(rows); T.to_csv('step8_susie_stability_summary.csv', index=False)
cols = [('haplotype_ncs', 'LTL CS\n(hap)'), ('haplotype_lowz', 'CS |z|<2\n(hap)'), ('haplotype_s', 's\n(hap)'), ('haplotype_lead', 'PP4 lead\n(hap)'), ('haplotype_max', 'PP4 max\n(hap)'),
        ('dosage_ncs', 'LTL CS\n(dos)'), ('dosage_lowz', 'CS |z|<2\n(dos)'), ('dosage_s', 's\n(dos)'), ('dosage_lead', 'PP4 lead\n(dos)'), ('dosage_max', 'PP4 max\n(dos)'),
        ('L20_ncs', 'LTL CS\n(L=20)'), ('L20_max', 'PP4 max\n(L=20)')]
fig = plt.figure(figsize=(W, 6.2)); ax = fig.add_axes([0.2, 0.08, 0.66, 0.8])
pp = LinearSegmentedColormap.from_list('pp', ['#f4f2ec', '#cdc4e8', VIOLET]); cs = LinearSegmentedColormap.from_list('cs', ['#f4f2ec', '#f6b79b', RED])
for i, r in T.iterrows():
    for j, (c, _) in enumerate(cols):
        v = r.get(c, np.nan)
        if pd.isna(v): col, txt = 'white', '—'
        elif 'ncs' in c: col, txt = cs(min(v, 10) / 10), f'{int(v)}'
        elif 'lowz' in c: col, txt = cs(min(v, 3) / 3), f'{int(v)}'
        elif c.endswith('_s'): col, txt = cs(min(v / 0.3, 1)), f'{v:.2f}'
        else: col, txt = pp(v), ('<0.01' if v < 0.005 else f'{v:.2f}')
        ax.add_patch(plt.Rectangle((j, i), 1, 1, color=col, ec='white', lw=0.8)); ax.text(j + 0.5, i + 0.5, txt, ha='center', va='center', fontsize=5.2, color='white' if (not pd.isna(v) and (('PP4' in cols[j][1] and v > 0.7))) else INK)
    ax.text(len(cols) + 0.15, i + 0.5, r.verdict, va='center', fontsize=5.4, color=AQUA if r.verdict == 'stable' else (MUTED if r.verdict.startswith('no') else RED), fontweight='bold')
ax.set_xlim(0, len(cols)); ax.set_ylim(len(T), 0); ax.set_xticks(np.arange(len(cols)) + 0.5); ax.set_xticklabels([c[1] for c in cols], fontsize=5.0); ax.xaxis.tick_top()
ax.set_yticks(np.arange(len(T)) + 0.5); ax.set_yticklabels([f'{r.locus} {TRL[r.trait]} ±{r.window_kb} kb' for r in T.itertuples()], fontsize=5.4); ax.tick_params(length=0)
for x in [5, 10]: ax.axvline(x, color=INK, lw=0.8)
ax.axhline(8, color=INK, lw=0.8)
for sp in ax.spines.values(): sp.set_visible(False)
ax.set_title('Multi-signal colocalization diagnostics (SuSiE, coloc.susie p12 = 10$^{-5}$): not interpreted in the main text', pad=30, fontsize=6.8)
fig.text(0.01, 0.01, 'hap: signed phased-haplotype LD; dos: genotype-dosage LD, both from the same 503 1000 Genomes European individuals. LTL CS = LTL credible sets with L = 10 (10 = all components used). '
         'CS |z|<2: credible sets led by a variant with LTL |z| < 2. s: estimate_s_rss LD-mismatch for LTL. PP4 lead: largest PP4 among pairs involving the LTL credible set containing the regional lead; PP4 max: over all pairs. '
         'Stable: LTL fit not saturated, no low-|z| sets with dosage LD, and lead-pair PP4 on the same side of 0.5 with both LD definitions.', fontsize=5.0, color=MUTED, wrap=True)
save(fig, 'FigS7_susie_diagnostics'); print(T[['locus', 'trait', 'window_kb', 'haplotype_lead', 'dosage_lead', 'verdict']].round(3).to_string())
