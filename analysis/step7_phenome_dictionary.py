"""v5: phenome-screen category dictionary and query inventory (E050). Categories were defined by the authors before the sensitivity
analysis was run in v4 but after the primary results were known; they are therefore 'selected', not preregistered."""
import json, re, pandas as pd, importlib.util
src = open('step6_phenome_screen.py').read(); ns = {}; exec(src[src.index('CATS = {'):src.index('EXCL = [')] , {'re': re}, ns)
CATS = ns['CATS']
d = json.load(open('source_data/kp_v4_instrument_phenome_and_mvp_regional.json'))['inst']
S = pd.read_csv('mr_v4_inputs/forward_selected_r2_0.001_S_mvpEU.csv')
recs = [(v, ds, ph, p) for v in S.varId for ds, ph, b, p, n in d.get(v, [])]
R = pd.DataFrame(recs, columns=['varId', 'dataset', 'phenotype', 'p'])
rows = [dict(category=c, regex=rx, used_for_exclusion=c != 'Blood pressure (not excluded)',
             phenotype_codes_matched=';'.join(sorted({p for p in R.phenotype.unique() if re.match(rx, str(p))}))) for c, rx in CATS.items()]
pd.DataFrame(rows).to_csv('step7_phenome_dictionary.csv', index=False)
inv = dict(retrieval_date='2026-09-29', endpoint='bioindex variant-dataset-associations (one query per instrument varId)', instruments_queried=len(S),
           instruments_with_any_record=int(R.varId.nunique()), records_stored=len(R), datasets=int(R.dataset.nunique()), phenotypes=int(R.phenotype.nunique()),
           records_P_lt_5e8_nonLTL=int(((R.p < 5e-8) & ~R.phenotype.isin(['LTL', 'TL'])).sum()),
           note='At retrieval, every returned record was filtered and stored if P < 5e-8 (any dataset) or if it belonged to the BMI, lymphocyte, WBC or outcome datasets; the screen therefore uses stored P < 5e-8 records. Absence of a record is not evidence of no association.')
pd.DataFrame([inv]).to_csv('step7_phenome_query_inventory.csv', index=False); print(inv)
