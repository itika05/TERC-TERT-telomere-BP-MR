"""v5 dataset register (Table 1 / Supplementary Table S1): one row per phenotype and cohort, with Knowledge Portal dataset ID and phenotype code,
role, maximum versus analyzed N, scale and UK Biobank overlap."""
import pandas as pd, numpy as np
MV = pd.read_csv('mr_v4_inputs/mvmr_instruments.csv')
rows = [
 ('Codd et al. 2021 (source summary statistics, figshare)', '—', 'LTL', 'European (UK Biobank)', '464,716', 'BOLT-LMM, z-standardised log T/S', 'Univariable forward-MR exposure; European primary reverse-MR outcome; colocalization', 'is UK Biobank'),
 ('Codd et al. 2021 via Knowledge Portal', 'Codd2021_LTL_Mixed', 'LTL', 'European', '472,174 (portal record)', 'as source', 'MVMR exposure (LTL column); comparison with source statistics', 'is UK Biobank'),
 ('MVP, Verma et al. 2024', 'Verma2024_MVPTraits_EU', 'SBP', 'European', '425,740', 'rank-based inverse-normal (SD)', 'Primary outcome; colocalization', 'none'),
 ('MVP, Verma et al. 2024', 'Verma2024_MVPTraits_EU', 'DBP', 'European', '425,743', 'rank-based inverse-normal (SD)', 'Primary outcome; colocalization', 'none'),
 ('MVP, Verma et al. 2024', 'Verma2024_MVPTraits_EU', 'HYPERTENSION', 'European', '318,398 (effective N in record)', 'log odds ratio', 'Primary outcome (hypertension); colocalization', 'none'),
 ('FinnGen R12', 'finngen_R12_I9_HYPTENS (not portal)', 'I9_HYPTENS', 'Finnish', '154,630 cases; 345,634 controls', 'log odds ratio', 'Secondary outcome; colocalization', 'none'),
 ('MVP, Verma et al. 2024', 'Verma2024_MVPTraits_AA', 'SBP; DBP; HYPERTENSION', 'African American', 'SBP 119,331; DBP 119,332 (portal record); hypertension N not in extract', 'inverse-normal (SD); log OR', 'Transportability', 'none'),
 ('MVP, Verma et al. 2024', 'Verma2024_MVPTraits_HS', 'SBP; DBP; HYPERTENSION', 'Hispanic', 'SBP 57,988; DBP 57,990 (portal record); hypertension N not in extract', 'inverse-normal (SD); log OR', 'Transportability', 'none'),
 ('Biobank Japan', 'GWAS_BBJ_ea', 'SBP; DBP', 'East Asian', '136,597; 136,615', 'portal record scale (implied SD: Supplementary Table S2)', 'Transportability', 'none'),
 ('Genes & Health, Huang et al.', 'Huang2021_Cardiometabolic_SA', 'SBP; DBP', 'British Pakistani/Bangladeshi', '18,536 (portal record)', 'portal record scale', 'Transportability', 'none'),
 ('Nakao et al. 2026', 'Nakao2026_LTL_SA', 'LTL', 'South Asian', '11,277', 'SD', 'Exploratory reverse-MR outcome', 'not established'),
 ('Pulit et al. 2019 (GIANT + UK Biobank; Zenodo 10.5281/zenodo.1251813, bmi...combined file)', 'GWAS_UKBiobankGIANT_eu', 'BMI', 'European', f'{int(MV.BMI_n.min()):,}–{int(MV.BMI_n.max()):,} per analyzed variant', 'SD', 'MVMR exposure', 'includes UK Biobank'),
 ('Chen et al. 2020 (Blood Cell Consortium)', 'Chen2020_BCX_eu', 'LymphoCount', 'European', f'{int(MV.LYM_n.min()):,}–{int(MV.LYM_n.max()):,} per analyzed variant', 'SD', 'MVMR exposure', 'includes UK Biobank'),
 ('Knowledge Portal ancestry meta-analyses', 'bottom-line ancestry-specific', 'SBP; DBP; HYPERTENSION', '6 groups', 'aggregate record N (up to 4.9 million; summed across studies, not unique participants)', 'mixed; SE derived from P', 'Descriptive variant look-up; comparison colocalization (Fig. S2b)', 'may include UK Biobank'),
 ('1000 Genomes phase 3', 'release 20130502', 'phased genotypes', 'EUR (503 individuals)', '503', '—', 'LD for clumping, GLS and SuSiE', 'none (public reference genotypes)'),
 ('GTEx v8', 'portal eQTL API', 'cis-eQTL', 'mostly European', '73–670 per tissue', 'NES', 'Annotation only', 'none')]
D = pd.DataFrame(rows, columns=['source', 'portal_dataset_id', 'phenotype_code', 'ancestry', 'N', 'scale', 'role', 'UK_Biobank_overlap'])
D['accessed'] = '2026-09-29'; D.to_csv('step8_dataset_register.csv', index=False); print(D[['source', 'phenotype_code', 'N', 'UK_Biobank_overlap']].to_string())
