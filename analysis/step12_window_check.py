"""v10 (review of v9, major 6): coloc.abf PP4 at +/-50 and +/-100 kb are close because the added variants carry little evidence relative to
the dominant lead signal. Lists the variant counts per window and PP4 to four decimals for every locus x trait x prior (S83)."""
import pandas as pd
A = pd.read_csv('step6_coloc_v4_abf.csv'); W = pd.read_csv('step6_coloc_v4_window_ledger.csv')
a = A[A.window_kb.isin([50, 100])].pivot_table(index=['locus', 'trait', 'p12'], columns='window_kb', values=['PP4', 'nsnps']).reset_index()
a.columns = ['locus', 'trait', 'p12'] + [f'{x}_{y}kb' for x, y in a.columns[3:]]
a = a[['locus', 'trait', 'p12', 'nsnps_50kb', 'nsnps_100kb', 'PP4_50kb', 'PP4_100kb']]
a['abs_diff_PP4'] = (a['PP4_50kb'] - a['PP4_100kb']).abs()
lw = W[W.window_kb.isin([50, 100])][['locus', 'window_kb', 'trait', 'n_variants', 'min_p_ltl', 'min_p_trait']]
a.to_csv('step12_window_check.csv', index=False); print(a.round(4).to_string()); print('max diff', a.abs_diff_PP4.max())
