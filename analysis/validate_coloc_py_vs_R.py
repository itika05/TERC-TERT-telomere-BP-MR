"""Python coloc_abf (mr_lib) vs R coloc 6.0.3 coloc.abf on identical v3 inputs."""
import pandas as pd, numpy as np
from mr_lib import coloc_abf
A = pd.read_csv('step4v3_coloc_abf_R.csv'); rows = []
for r in A.itertuples():
    sec = 'se' if r.se_source.startswith('P') else 'se_reported'
    L = pd.read_csv(f'coloc_v3_inputs/{r.locus}_{r.window_kb}kb_LTL.csv'); B = pd.read_csv(f'coloc_v3_inputs/{r.locus}_{r.window_kb}kb_{r.trait}.csv')
    res, _ = coloc_abf(L.beta.values, L[sec].values, B.beta.values, B[sec].values, sd1=0.15, sd2=0.2 if r.trait == 'HYPERTENSION' else 0.15, p12=r.p12)
    rows.append(dict(locus=r.locus, window_kb=r.window_kb, se_source=r.se_source, trait=r.trait, p12=r.p12, PP4_R=r.PP4, PP4_py=res['PP4'], PP3_R=r.PP3, PP3_py=res['PP3']))
V = pd.DataFrame(rows); V['max_abs_diff'] = np.maximum((V.PP4_R - V.PP4_py).abs(), (V.PP3_R - V.PP3_py).abs())
V.to_csv('validation/coloc_python_vs_R.csv', index=False); print(len(V), 'max |diff|', V.max_abs_diff.max()); print(V.sort_values('max_abs_diff').tail(5).round(4).to_string())
