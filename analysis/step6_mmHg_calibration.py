"""Scale of the MVP blood-pressure estimates (audit concern 2).
Verma et al. (2024) rank-based inverse-normal transformed quantitative traits, so MVP BP effects are in SD units of the transformed trait.
The Knowledge Portal also stores, for datasets reported in mmHg, the original ('unscaled') beta next to the portal-standardised beta.
Calibration: for variants among the top 1,500 portal associations of both MVP (European) and a mmHg-scale European GWAS (Evangelou 2018; Keaton 2024),
regress the MVP beta on the mmHg beta through the origin (weights 1/SE_MVP^2); the inverse slope is the number of mmHg per MVP unit."""
import json, numpy as np, pandas as pd
d = json.load(open('source_data/kp_v4_round2.json')); cal = d['cal']
rows = []
for ph in ['SBP', 'DBP']:
    mv = {r[0]: r for r in cal[f'Verma2024_MVPTraits_EU|{ph}']}
    for ds in ['Evangelou2018_bp_eu', 'Keaton2024_BP_EU']:
        ref = cal[f'{ds}|{ph}']
        pairs = [(mv[r[0]], r) for r in ref if r[0] in mv]
        b_mv = np.array([p[0][1] for p in pairs]); se_mv = np.array([p[0][2] for p in pairs]); b_mm = np.array([p[1][6] for p in pairs])
        scale = np.median([r[6] / r[1] for r in ref])
        w = 1 / se_mv ** 2; slope = (w * b_mv * b_mm).sum() / (w * b_mm ** 2).sum()
        # bootstrap CI over variants
        rng = np.random.default_rng(1); bs = []
        for _ in range(2000):
            i = rng.integers(0, len(b_mv), len(b_mv)); bs.append((w[i] * b_mv[i] * b_mm[i]).sum() / (w[i] * b_mm[i] ** 2).sum())
        lo, hi = np.percentile(1 / np.array(bs), [2.5, 97.5])
        rows.append(dict(trait=ph, reference=ds, n_shared_variants=len(pairs), portal_scale_mmHg_per_SD=scale, mmHg_per_MVP_unit=1 / slope,
                         boot_lo=lo, boot_hi=hi, r=np.corrcoef(b_mv, b_mm)[0, 1]))
C = pd.DataFrame(rows); C.to_csv('step6_mmHg_calibration.csv', index=False); print(C.round(3).to_string())
