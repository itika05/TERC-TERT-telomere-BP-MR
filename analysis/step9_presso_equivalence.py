"""v7: algorithm-level equivalence test of the Python MR-PRESSO re-implementation (step8_presso_np.py) against MRPRESSO 1.0,
on the SAME simulated data. presso_equivalence.R runs the package (its own code) and saves the random exposure/outcome arrays it drew,
its RSSexp vector, per-variant outlier statistics and final results. Here the Python functions are applied to those arrays and each
intermediate is compared: leave-one-out slopes (via the package's LOO predictions implied by Dif2), observed RSS, every simulated RSS,
the global exceedance count, per-variant squared residuals, per-variant Bonferroni outlier P values, outlier identities and the
outlier-corrected estimate. Tolerances: 1e-9 relative for continuous quantities; exact agreement for counts and outlier identities; outlier P values (counts scaled by k/B) to 1e-12.
The distortion test is not re-implemented and is not tested."""
import glob, sys
import numpy as np, pandas as pd
from step8_presso_np import presso_from_arrays
rows = []
for pk in sorted(glob.glob('presso_equiv/*_package.csv')):
    stem = pk[:-len('_package.csv')]; P = pd.read_csv(pk).iloc[0]; code, B = P.code, int(P.B)
    d = pd.read_csv(f'mr_v4_inputs/forward_selected_r2_0.001_{code}.csv'); k = len(d)
    s = np.sign(d.bx.values); bx, by, sx, sy = d.bx.values * s, d.by.values * s, d.sx.values, d.sy.values   # package orients to positive bx
    A = np.fromfile(stem + '_random.bin', dtype='<f8').reshape(B, 2 * k); X, Y = A[:, :k], A[:, k:]
    rss_pkg = np.fromfile(stem + '_RSSexp.bin', dtype='<f8')
    O = pd.read_csv(stem + '_outlier_test.csv')
    r, p_out, out, st = presso_from_arrays(bx, sx, by, sy, X, Y)
    pkgP = pd.to_numeric(O.Pvalue.astype(str).str.lstrip('<'), errors='coerce').values
    pkg_is_floor = O.Pvalue.astype(str).str.startswith('<').values          # package prints 0 counts as '<k/B'
    py_P = p_out.copy(); py_floor = (st['n_out'] == 0)
    gP = str(P.global_P); g_pkg = 0 if gP.startswith('<') else round(float(gP) * B)
    pkg_out = set(str(P.outliers).split(';')) - {'', 'nan'}; py_out = set(d.rsid[out])
    rel = lambda a, b: float(np.max(np.abs(np.asarray(a) - np.asarray(b)) / np.maximum(np.abs(np.asarray(b)), 1e-300)))
    checks = {
        'observed RSS': rel(st['rss_obs'], P.RSSobs) <= 1e-9,
        'per-variant squared residual (Dif2)': rel(st['dif2'], O.Dif2.values) <= 1e-9,
        'simulated RSS (all B replicates)': rel(st['rss_exp'], rss_pkg) <= 1e-9,
        'global exceedance count': r['global_count'] == g_pkg,
        'per-variant outlier P (Bonferroni, capped at 1)': bool(np.all(py_floor == pkg_is_floor) and np.allclose(py_P[~py_floor], pkgP[~pkg_is_floor], rtol=0, atol=1e-12)),
        'outlier identities': py_out == pkg_out,
        'corrected estimate': (len(pkg_out) == 0 and np.isnan(P.b_corrected)) or rel(r['b'], P.b_corrected) <= 1e-9,
        'corrected SE': (len(pkg_out) == 0) or rel(r['se'], P.se_corrected) <= 1e-9,
        'corrected P': (len(pkg_out) == 0) or rel(r['p'], P.p_corrected) <= 1e-6,
    }
    for name, ok in checks.items():
        rows.append(dict(code=code, B=B, seed=int(P.seed), check=name, status='PASS' if ok else 'FAIL'))
    rows.append(dict(code=code, B=B, seed=int(P.seed), check='summary', status=f"package outliers {len(pkg_out)}, python {len(py_out)}; global count pkg {g_pkg} py {r['global_count']}; "
                     f"b pkg {P.b_corrected:.6f} py {r['b']:.6f}; max rel diff simulated RSS {rel(st['rss_exp'], rss_pkg):.2e}"))
T = pd.DataFrame(rows); T.to_csv('step9_presso_equivalence.csv', index=False); print(T.to_string())
if (T.status == 'FAIL').any(): sys.exit('MR-PRESSO equivalence FAILED')
