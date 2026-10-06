"""v7: Supplementary Table S21 holds only the reported 318-variant MVMR analysis (the 330-variant run with fallback-SE records is not reported)."""
import pandas as pd
M = pd.read_csv('step8_mvmr_v6_estimates.csv'); P = M[M.set.str.startswith('primary')]
assert (P.nsnp == 318).all(); P.to_csv('step9_mvmr_estimates_primary.csv', index=False); print(len(P), 'rows')
