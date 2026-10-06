"""v6: genotype-dosage LD (correlation of 0/1/2 allele counts in the same 503 individuals) versus phased-haplotype LD, for the GLS matrix
(173 variants) and the TERC/TERT colocalization windows; writes dosage LD matrices for GLS and SuSiE sensitivity analyses."""
import numpy as np, pandas as pd
from ld_ref import load, match, hap_matrix, corr
def dos_corr(H): G = H[:, 0::2] + H[:, 1::2]; return corr(G)
rows = []
REF = load('ldref/eur_instruments.json')
C = pd.read_csv('mr_v4_inputs/forward_r2_0.1_correlation.csv', index_col=0)
m = match(C.index, REF).set_index('varId').loc[C.index].reset_index(); H = hap_matrix(m, REF)
Rh, Rd = corr(H), dos_corr(H); assert np.allclose(Rh, C.values, atol=1e-6)
pd.DataFrame(Rd, index=C.index, columns=C.index).to_csv('mr_v6_inputs/forward_r2_0.1_correlation_dosage.csv')
off = np.triu_indices_from(Rh, 1); dd = np.abs(Rh[off] - Rd[off])
rows.append(dict(matrix='GLS r2<0.1 set', n_variants=len(C), median_abs_diff=np.median(dd), p99_abs_diff=np.percentile(dd, 99), max_abs_diff=dd.max(),
                 sign_disagreements_abs_r_gt_0_05=int(((np.sign(Rh[off]) != np.sign(Rd[off])) & (np.abs(Rh[off]) > 0.05)).sum()), min_eig_dosage=np.linalg.eigvalsh(Rd).min()))
for loc in ['TERC', 'TERT']:
    ref = load(f'ldref/eur_region_{loc}.json')
    for w in [100, 250]:
        for tr in ['SBP', 'DBP', 'HYPERTENSION', 'FG_HYPERTENSION']:
            L = pd.read_csv(f'coloc_v4_inputs/{loc}_{w}kb_{tr}_LD.csv', index_col=0)
            mm = match(L.index, ref).set_index('varId').loc[L.index].reset_index(); H = hap_matrix(mm, ref)
            Rh, Rd = corr(H), dos_corr(H); assert np.allclose(Rh, L.values, atol=1e-5)
            pd.DataFrame(Rd, index=L.index, columns=L.index).to_csv(f'coloc_v6_inputs/{loc}_{w}kb_{tr}_LDdos.csv')
            off = np.triu_indices_from(Rh, 1); dd = np.abs(Rh[off] - Rd[off])
            rows.append(dict(matrix=f'{loc} ±{w} kb ({tr})', n_variants=len(L), median_abs_diff=np.median(dd), p99_abs_diff=np.percentile(dd, 99), max_abs_diff=dd.max(),
                             sign_disagreements_abs_r_gt_0_05=int(((np.sign(Rh[off]) != np.sign(Rd[off])) & (np.abs(Rh[off]) > 0.05)).sum()), min_eig_dosage=np.linalg.eigvalsh(Rd).min()))
D = pd.DataFrame(rows); D.to_csv('step8_ld_haplotype_vs_dosage.csv', index=False); print(D.round(4).to_string())
