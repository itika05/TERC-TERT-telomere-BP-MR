"""Step 2 sensitivity analyses requested by review: (a) reported SE instead of P-derived SE; (b) include the record excluded by the
SE-consistency rule; (c) drop 'AF' (possible overlap with 'AA' in the portal's ancestry meta-analyses); (d) European only; (e) ancestry-specific results first."""
import pandas as pd, numpy as np
from step2_lookup_pool import pool
a = pd.read_csv('kp_snp_lookup_by_ancestry_checked.csv')
anc6 = ['EU', 'SA', 'EA', 'AF', 'HS', 'AA']
print('excluded record(s):'); print(a[(~a.qc_pass) & a.ancestry.isin(anc6) & a.trait.isin(['SBP', 'DBP', 'HYPERTENSION'])][['rs', 'ancestry', 'trait', 'beta_alt', 'se', 'p', 'se_p', 'n']].to_string())
rows = []
for rs in ['rs10936599', 'rs2736100']:
    for tr in ['SBP', 'DBP', 'HYPERTENSION']:
        base = a[(a.rs == rs) & (a.trait == tr) & a.ancestry.isin(anc6)]
        for lab, x, secol in [('primary: P-derived SE, QC pass', base[base.qc_pass], 'se_p'),
                              ('reported SE, QC pass', base[base.qc_pass], 'se'),
                              ('P-derived SE, all records', base, 'se_p'),
                              ('reported SE, all records', base, 'se'),
                              ('P-derived SE, QC pass, without AF', base[base.qc_pass & (base.ancestry != 'AF')], 'se_p'),
                              ('P-derived SE, QC pass, without AA', base[base.qc_pass & (base.ancestry != 'AA')], 'se_p'),
                              ('European only', base[base.ancestry == 'EU'], 'se_p')]:
            r = pool(x.beta_alt.values, x[secol].values); r.update(rs=rs, trait=tr, analysis=lab, ancestries=';'.join(x.ancestry)); rows.append(r)
out = pd.DataFrame(rows)[['rs', 'trait', 'analysis', 'ancestries', 'k', 'b_random', 'se_random', 'p_random', 'I2']]
out.to_csv('step2_pooling_sensitivity.csv', index=False)
print(out.round(5).to_string())
