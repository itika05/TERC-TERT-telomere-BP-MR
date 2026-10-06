"""Unit check for every summary-statistic dataset (audit major concern 2).
For a quantitative trait, SE(beta) ~= SD_y / sqrt(2 N f (1-f)) (small per-variant R2). Hence the implied trait SD is
SD_y ~= SE * sqrt(2 N f (1-f)). A dataset on a standardised scale gives ~1; a dataset in mmHg gives ~15-20 (SBP) or ~10 (DBP).
SE is |beta|/z(P) from the portal record; f is the 1000 Genomes EUR (or SAS) frequency, which is a rough proxy for non-European samples,
so only the order of magnitude is interpreted. Only instruments with 1e-300 < P < 0.5 are used (stable z)."""
import numpy as np, pandas as pd
from scipy.stats import norm
d = pd.read_csv('kp_ltl_leads_wide.csv'); d = d[d.varId.str.count(':') == 3]
sets = [('L_codd', 'LTL', 'Codd 2021 UK Biobank', 464716, 'afEU'), ('S_mvpEU', 'SBP', 'MVP European', 425740, 'afEU'), ('D_mvpEU', 'DBP', 'MVP European', 425743, 'afEU'),
        ('S_mvpAA', 'SBP', 'MVP African American', 119331, 'afEU'), ('D_mvpAA', 'DBP', 'MVP African American', 119332, 'afEU'),
        ('S_mvpHS', 'SBP', 'MVP Hispanic', 57988, 'afEU'), ('D_mvpHS', 'DBP', 'MVP Hispanic', 57990, 'afEU'),
        ('S_bbj', 'SBP', 'Biobank Japan', 136597, 'afEU'), ('D_bbj', 'DBP', 'Biobank Japan', 136615, 'afEU'),
        ('S_gh', 'SBP', 'Genes & Health', 18536, 'afSA'), ('D_gh', 'DBP', 'Genes & Health', 18536, 'afSA')]
rows = []
for c, t, s, n, afc in sets:
    x = d[[c + '_b', c + '_p', afc]].dropna(); x = x[(x[c + '_p'] > 1e-300) & (x[c + '_p'] < 0.5) & (x[c + '_b'] != 0)]
    f = x[afc].clip(0.01, 0.99); se = x[c + '_b'].abs() / norm.isf(x[c + '_p'] / 2); sd = se * np.sqrt(2 * n * f * (1 - f))
    rows.append(dict(code=c, trait=t, dataset=s, N_used=n, af_source='1000G ' + ('EUR' if afc == 'afEU' else 'SAS'), n_variants=len(x),
                     implied_SD_median=np.median(sd), implied_SD_IQR_low=np.percentile(sd, 25), implied_SD_IQR_high=np.percentile(sd, 75),
                     interpretation='standardised (SD) scale' if np.median(sd) < 3 else 'natural units (e.g. mmHg)'))
U = pd.DataFrame(rows); U.to_csv('step6_unit_check.csv', index=False); print(U.round(3).to_string())
