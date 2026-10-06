"""v6: tests of the SE-from-P implementation. The code uses scipy.stats.norm.isf(P/2) (inverse survival function evaluated directly in the upper
tail; equivalent to R qnorm(P/2, lower.tail = FALSE)), not norm.ppf(1 - P/2), which loses precision when P/2 < ~1e-16."""
import numpy as np, pandas as pd, subprocess
from scipy.stats import norm
P = [0.999, 0.9, 0.5, 0.05, 1e-5, 1e-10, 1e-16, 1e-20, 1e-50, 1e-100, 1e-200, 1e-300]
rcode = 'cat(sprintf("%.15g", qnorm(c(' + ','.join(repr(p) for p in P) + ')/2, lower.tail = FALSE)), sep = "\\n")'
r = subprocess.run(['Rscript', '-e', rcode], capture_output=True, text=True).stdout.split()
rows = []
for p, rq in zip(P, r):
    isf = norm.isf(p / 2); naive = norm.ppf(1 - p / 2)
    rows.append(dict(P=p, z_isf_python=isf, z_qnorm_R_upper_tail=float(rq), rel_diff_python_vs_R=abs(isf - float(rq)) / float(rq),
                     z_naive_ppf_1_minus_P_over_2=naive, naive_fails=not np.isfinite(naive) or abs(naive - isf) / isf > 1e-6,
                     se_for_beta_0_01=0.01 / isf))
T = pd.DataFrame(rows); T.to_csv('step8_se_formula_tests.csv', index=False); print(T.to_string())
assert (T.rel_diff_python_vs_R < 1e-12).all()
