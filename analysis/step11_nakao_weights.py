"""v9 (review M2): the 116 primary instruments re-weighted with the Nakao et al. 2026 European LTL GWAS (Knowledge Portal dataset
Nakao2026_LTL_EU; UK Biobank + All of Us European-like participants). The Nakao records for every instrument were retrieved on
29 September 2026 in the per-instrument portal query used for the phenome screen (source_data/kp_v4_instrument_phenome_and_mvp_regional.json,
key 'inst'); no new data were retrieved. Portal effects are on the portal alternate allele, the same allele as the primary bx, so no
re-alignment is needed (checked: the varId key is identical). SE = |beta| / z(P), as for every portal record.
Instrument selection is NOT redone in Nakao: only the weights change (no genome-wide Nakao statistics are held)."""
import json, numpy as np, pandas as pd, os
from scipy.stats import norm
J = json.load(open('source_data/kp_v4_instrument_phenome_and_mvp_regional.json'))['inst']
os.makedirs('mr_v9_inputs', exist_ok=True); led = []
for code in ['S_mvpEU', 'D_mvpEU', 'H_mvpEU', 'H_fg12']:
    d = pd.read_csv(f'mr_v4_inputs/forward_selected_r2_0.001_{code}.csv')
    nb, np_, nn = [], [], []
    for v in d.varId:
        rec = [r for r in J.get(v, []) if r[0] == 'Nakao2026_LTL_EU' and r[1] == 'LTL']
        assert len(rec) == 1, (v, rec)
        nb.append(rec[0][2]); np_.append(rec[0][3]); nn.append(rec[0][4])
    d['bx_codd'], d['sx_codd'] = d.bx, d.sx
    d['bx'] = nb; d['p_nakao'] = np_; d['n_nakao'] = nn
    d['sx'] = np.abs(d.bx) / norm.isf(d.p_nakao / 2)
    assert (d.p_nakao > 0).all() and (d.bx != 0).all()
    d['same_sign_as_codd'] = np.sign(d.bx) == np.sign(d.bx_codd)
    d.to_csv(f'mr_v9_inputs/forward_nakaoEU_weights_{code}.csv', index=False)
    w = 1 / d.sx_codd ** 2; slope = float((w * d.bx_codd * d.bx).sum() / (w * d.bx_codd ** 2).sum())
    led.append(dict(code=code, k=len(d), same_sign=int(d.same_sign_as_codd.sum()), median_n_nakao=float(d.n_nakao.median()),
                    max_p_nakao=float(d.p_nakao.max()), slope_nakao_on_codd=slope, r_nakao_codd=float(np.corrcoef(d.bx, d.bx_codd)[0, 1]),
                    median_F_nakao=float(np.median((d.bx / d.sx) ** 2)), median_F_codd=float(np.median((d.bx_codd / d.sx_codd) ** 2))))
L = pd.DataFrame(led); L.to_csv('step11_nakao_weights_ledger.csv', index=False); print(L.to_string())
