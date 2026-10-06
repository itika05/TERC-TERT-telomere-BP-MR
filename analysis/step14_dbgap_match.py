"""v12 (review of v11, major 3): match each MVP portal record (Verma2024_MVPTraits_*) to the dbGaP phs002453.v1.p1 analysis
with the same N. dbGaP analysis list and descriptions retrieved 3 October 2026 (step14_dbgap_raw.csv); for binary traits the
portal N is compared with 4/(1/cases + 1/controls). Writes step14_dbgap_mvp_match.csv (Supplementary Table S86)."""
import pandas as pd
raw = pd.read_csv('step14_dbgap_raw.csv')
raw['pop'] = raw.analysis.str.split('.').str[1]; raw['pheno'] = raw.analysis.str.split('.').str[0]
raw['N_compare'] = [4 / (1 / a + 1 / b) if pd.notna(a) else n for n, a, b in zip(raw.total_N, raw.cases, raw.controls)]
kp = pd.read_csv('kp_snp_lookup_by_dataset.csv'); kp = kp[kp.dataset.str.match(r'Verma2024_MVPTraits_(EU|AA|HS)$')]
P = kp.groupby(['dataset', 'trait']).n.agg(['min', 'max']).reset_index()
popmap = {'EU': 'EUR', 'AA': 'AFR', 'HS': 'AMR'}
fam = {'SBP': 'Systolic_', 'DBP': 'Diastolic_', 'HYPERTENSION': 'Phe_401', 'EssentialHYPERTENSION': 'Phe_401_1'}
rows = []
for _, r in P.iterrows():
    pop = popmap[r.dataset.split('_')[-1]]; f = fam[r.trait]
    cand = raw[(raw['pop'] == pop) & ((raw.pheno.str.startswith(f)) if r.trait in ('SBP', 'DBP') else (raw.pheno == f))].copy()
    cand['abs_diff'] = (cand.N_compare - r['min']).abs()
    for _, c in cand.sort_values('abs_diff').iterrows():
        rows.append(dict(portal_dataset=r.dataset, portal_trait=r.trait, portal_N=r['min'], portal_N_constant=r['min'] == r['max'],
                         dbgap_analysis=c.analysis, dbgap_accession=c.pha, dbgap_total_N=int(c.total_N),
                         dbgap_cases=c.cases, dbgap_controls=c.controls, dbgap_N_compared=round(c.N_compare, 2),
                         abs_diff=round(c.abs_diff, 2), match=c.abs_diff < 0.5))
O = pd.DataFrame(rows); O.to_csv('step14_dbgap_mvp_match.csv', index=False)
print(O[O.match].to_string(index=False)); print('unique match per record:', O.groupby(['portal_dataset', 'portal_trait']).match.sum().to_dict())
