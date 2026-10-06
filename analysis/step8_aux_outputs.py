"""v6 auxiliary outputs derived from files already in the package (no new data):
step8_reverse_mde.csv, step8_published_vs_portal_pool.csv, step8_gtex_terc_query.csv, step8_reference_check.csv,
step8_terc_credible_set_tracking.csv, step8_terc_coloc_pairs_tracking.csv.
The GTEx and Europe PMC JSON files were retrieved in a browser on 2026-09-29 (source_data/); this script only tabulates them."""
import json, glob
import numpy as np, pandas as pd
from scipy.stats import norm, t as tdist

# 1. reverse-MR minimum detectable effects (80% power, alpha 0.05 and 0.05/6), IVW SE of the primary reverse analyses
R = pd.read_csv('step8_mr_results_v6.csv')
rows = []
for code in ['SBP->L_nkSA', 'DBP->L_nkSA', 'Hypertension->L_nkSA', 'SBP->L_codd', 'DBP->L_codd', 'Hypertension->L_codd']:
    r = R[(R.code == code) & (R.method == 'IVW (MRE, t)') & (R.analysis == 'r2<0.001 (primary)')].iloc[0]
    for al in (0.05, 0.05 / 6):
        rows.append(dict(code=code, alpha=al, k=int(r.nsnp), se=r.se, mde=(tdist.ppf(1 - al / 2, r.nsnp - 1) + norm.ppf(0.8)) * r.se))
pd.DataFrame(rows).to_csv('step8_reverse_mde.csv', index=False)

# 2. published Codd et al. MR instrument versus the primary (portal-lead) instruments
P = pd.read_csv('step7_published_instrument_reasons.csv')
S = pd.read_csv('source_data/codd2021_sentinels_parsed.csv'); S = S[S.MR == 'Y'][['Sentinel snp', 'MAF']].rename(columns={'Sentinel snp': 'rsid', 'MAF': 'maf'})
sel = pd.read_csv('mr_v4_inputs/snp_level_S_mvpEU.csv')  # the 116 primary instruments
prim = sel.varId.str.split(':', expand=True).iloc[:, :2].astype(int); prim.columns = ['chr', 'pos']
out = P[['rsid', 'chr', 'pos', 'gene']].merge(S.drop_duplicates('rsid'), on='rsid', how='left')
out['maf_lt_1pct'] = out.maf < 0.01
out['chrX'] = out.chr.astype(str).isin(['X', '23'])
out['same_variant'] = out.rsid.isin(sel.rsid)
def dist(r):
    if r.chrX: return np.nan
    d = (prim[prim.chr == int(r.chr)].pos - int(r.pos)).abs()
    return d.min() if len(d) else np.nan
out['dist_to_nearest_primary_instrument'] = out.apply(dist, axis=1)
out['within_500kb'] = out.dist_to_nearest_primary_instrument <= 500000
out.to_csv('step8_published_vs_portal_pool.csv', index=False)

# 3. GTEx v8 TERC query (archived JSON)
g = json.load(open('source_data/gtex_v8_TERC_query_2026-09-29.json'))
gr = g['gene_reference']; gr = gr[0] if isinstance(gr, list) else gr
items = [('gencodeId', gr.get('gencodeId')), ('geneType', gr.get('geneType')),
         ('singleTissueEqtl records returned', g['singleTissueEqtl']['records_returned'])]
med = g['medianGeneExpression_TPM']
for k, v in (med.items() if isinstance(med, dict) else [(m['tissueSiteDetailId'], m['median']) for m in med]):
    if k != 'endpoint': items.append((f'median TPM {k}', v))
items += [('query date', g['query_date']), ('API', g['api'])]
pd.DataFrame(items, columns=['item', 'value']).to_csv('step8_gtex_terc_query.csv', index=False)

# 4. reference metadata from Europe PMC (archived JSON)
e = json.load(open('source_data/europepmc_reference_check_2026-09-29.json'))
pd.DataFrame([dict(key=k, year=v[0], volume=v[1], issue=v[2], pages=v[3]) for k, v in e.items() if not k.startswith('_')]).to_csv('step8_reference_check.csv', index=False)

# 5. TERC credible-set and coloc.susie pair tracking across windows and LD definitions
CS = pd.read_csv('step8_susie_ld_sensitivity_credible_sets.csv')
CS[(CS.locus == 'TERC') & CS.pair.isin(['SBP', 'DBP', 'FG_HYPERTENSION'])].sort_values(['pair', 'trait', 'ld', 'window_kb'], kind='stable').to_csv('step8_terc_credible_set_tracking.csv', index=False)
CO = pd.read_csv('step8_susie_ld_sensitivity_coloc.csv')
co = CO[(CO.locus == 'TERC') & CO.trait.isin(['DBP', 'FG_HYPERTENSION'])].copy()
rs = dict(zip(CS.lead_varId, CS.lead_rsid))
co['hit1_rsid'] = co.hit1.map(rs); co['hit2_rsid'] = co.hit2.map(rs)
co.to_csv('step8_terc_coloc_pairs_tracking.csv', index=False)
print('aux outputs written')
