"""Phenome screen of the 116 forward instruments (audit major concern 4).
Source: Knowledge Portal variant-dataset-associations for each instrument (all datasets), retrieved 29 Sep 2026; every record with P < 5e-8 is kept.
Categories screened as potential confounding or pleiotropic pathways: adiposity, lipids, glycaemia/type 2 diabetes, kidney function,
leukocyte counts, clonal haematopoiesis/mosaic chromosomal alterations, and blood pressure itself (reported, not used for exclusion).
Outputs: step6_phenome_screen_long.csv (every P < 5e-8 record), step6_phenome_screen_summary.csv (instrument x category),
         mr_v4_inputs/forward_phenome_clean_{code}.csv and forward_no_terc_tert_phenome_clean_{code}.csv."""
import json, re, pandas as pd, numpy as np

CATS = {
    'Adiposity': r'^(BMI|WEIGHT|WHR|WHRadjBMI|Obesity|WAIST|HIPC|HIPCadjBMI|BodyFat.*|Ap-LM|TB-LM|LiverFat)$',
    'Lipids': r'^(HDL|LDL|CHOL|TG|nonHDL|ApoA|ApoA1|ApoB|TGtoHDL|Hyperlipidemia|HyperChol|Hyperglyceridemia|Dyslipid.*|RemnantC|IDL.*|TGnonT2D)$',
    'Glycaemia / T2D': r'^(HBA1C|HBA1CadjBMI|T2D.*|BS|BSandFG|BSadjFastingTime|FG|FGadjBMI|PI|PIadjBMI|MetS|MetS_cont)$',
    'Kidney function': r'^(eGFR.*|Creatinine|CystatinC|BUN|SerumUrea|CKD|AcuteRenalFail|CRF)$',
    'Leukocyte counts': r'^(WBC|LymphoCount|LymphoPerc|NeutCount|NeutroPerc|MonoCount|MonoPerc|EosinCount|EosinPerc|BasoCount|BasoPerc)$',
    'Clonal haematopoiesis / mCA': r'^(CHIP.*|mCAY|mCA.*)$',
    'Blood pressure (not excluded)': r'^(SBP|DBP|PulsePress|MAP|HYPERTENSION|EssentialHYPERTENSION)$',
}
EXCL = ['Adiposity', 'Lipids', 'Glycaemia / T2D', 'Kidney function', 'Leukocyte counts', 'Clonal haematopoiesis / mCA']

d = json.load(open('source_data/kp_v4_instrument_phenome_and_mvp_regional.json'))['inst']
S = pd.read_csv('mr_v4_inputs/forward_selected_r2_0.001_S_mvpEU.csv')
rows = []
for v in S.varId:
    for ds, ph, b, p, n in d.get(v, []):
        if p is not None and p < 5e-8 and ph not in ('LTL', 'TL'):
            cat = next((c for c, rx in CATS.items() if re.match(rx, ph)), 'Other')
            rows.append(dict(varId=v, dataset=ds, phenotype=ph, category=cat, beta=b, p=p, n=n))
L = pd.DataFrame(rows).merge(S[['varId', 'rsid', 'gene']], on='varId')
L.to_csv('step6_phenome_screen_long.csv', index=False)
summ = L.groupby(['varId', 'category']).p.min().unstack()
summ = S[['varId', 'rsid', 'gene', 'terc_tert_region']].merge(summ, left_on='varId', right_index=True, how='left')
summ['any_excluding_category'] = summ[[c for c in EXCL if c in summ]].notna().any(axis=1)
summ.to_csv('step6_phenome_screen_summary.csv', index=False)
print('instruments with >=1 P<5e-8 non-LTL association:', L.varId.nunique(), 'of', len(S))
print(summ[[c for c in CATS if c in summ]].notna().sum().to_string())
print('flagged for exclusion:', int(summ.any_excluding_category.sum()))
flag = set(summ.varId[summ.any_excluding_category])
for code in ['S_mvpEU', 'D_mvpEU', 'H_mvpEU']:
    x = pd.read_csv(f'mr_v4_inputs/forward_selected_r2_0.001_{code}.csv')
    x[~x.varId.isin(flag)].to_csv(f'mr_v4_inputs/forward_phenome_clean_{code}.csv', index=False)
    x[~x.varId.isin(flag) & ~x.terc_tert_region].to_csv(f'mr_v4_inputs/forward_no_terc_tert_phenome_clean_{code}.csv', index=False)
    print(code, 'clean', (~x.varId.isin(flag)).sum(), 'clean & no TERC/TERT', (~x.varId.isin(flag) & ~x.terc_tert_region).sum())
