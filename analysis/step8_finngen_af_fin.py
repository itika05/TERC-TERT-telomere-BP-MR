"""v6: FinnGen allele-frequency audit refined with 1000 Genomes FIN. The 503 EUR samples in the LD reference are in the order of the phase 3 VCF header
(HG00096 ... NA20832); population labels from integrated_call_samples_v3.20130502.ALL.panel (run-length: GBR 55, FIN 17, GBR 31, FIN 82, GBR 1,
IBS 94, GBR 4, IBS 13, CEU 99, TSI 107). The ordering is checked against Ensembl 1000 Genomes FIN frequencies for rs10936599 and rs2736100.
Rule for palindromic (A/T, C/G) instruments: retained if EUR minor allele frequency <= 0.42 and FinnGen and 1000 Genomes frequencies lie on the same
side of 0.5 (for both EUR and FIN); otherwise 'ambiguous' and excluded in the FinnGen sensitivity analysis. Non-palindromic instruments are aligned by
allele identity; |AF difference| > 0.10 is flagged and compared with 1000 Genomes FIN."""
import numpy as np, pandas as pd
from ld_ref import load, match, hap_matrix
rle = [('GBR', 55), ('FIN', 17), ('GBR', 31), ('FIN', 82), ('GBR', 1), ('IBS', 94), ('GBR', 4), ('IBS', 13), ('CEU', 99), ('TSI', 107)]
pops = np.array([p for p, n in rle for _ in range(n)]); assert len(pops) == 503
fin = np.repeat(pops == 'FIN', 2)
def fin_af(varids, ref):
    m = match(varids, ref); ok = m.idx >= 0; out = pd.Series(np.nan, index=m.varId)
    if ok.any():
        H = hap_matrix(m[ok], ref); out[m.varId[ok].values] = H[:, fin].mean(1)
    return out
# ordering check
E = pd.read_csv('ensembl_allele_frequencies.csv'); chk = []
for rs, loc in [('rs10936599', 'TERC'), ('rs2736100', 'TERT')]:
    ref = load(f'ldref/eur_region_{loc}.json'); x = pd.read_csv(f'coloc_v4_inputs/{loc}_250kb_SBP.csv'); v = x[x.rsid == rs].varId.iloc[0]
    alt = v.split(':')[3]; f = fin_af([v], ref).iloc[0]
    e = E[(E.rs == rs) & (E.population == '1000GENOMES:phase_3:FIN') & (E.allele == alt)].freq
    chk.append(dict(rsid=rs, allele=alt, fin_af_from_reference=f, ensembl_1000G_FIN=e.iloc[0] if len(e) else np.nan))
print(pd.DataFrame(chk))
A = pd.read_csv('step7_finngen_harmonisation_audit.csv'); REF = load('ldref/eur_instruments.json')
A['fin_1000G_af_kp_alt'] = fin_af(A.varId, REF).reindex(A.varId).values
A['af_diff_vs_FIN'] = (A.finngen_af_kp_alt - A.fin_1000G_af_kp_alt).abs()
pal = A.palindromic
A['palindrome_rule'] = np.where(~pal, 'not palindromic', np.where((A.eur_af_kp_alt - .5).abs() <= .08, 'ambiguous (EUR MAF > 0.42): excluded in sensitivity',
                         np.where(((A.finngen_af_kp_alt - .5) * (A.eur_af_kp_alt - .5) > 0) & ((A.finngen_af_kp_alt - .5) * (A.fin_1000G_af_kp_alt - .5) > 0), 'retained: same side of 0.5 in EUR, FIN and FinnGen', 'ambiguous: frequency side disagrees')))
A.to_csv('step8_finngen_harmonisation_audit_v6.csv', index=False)
print(A.palindrome_rule.value_counts()); print(A[A.af_diff > 0.1][['rsid', 'gene', 'palindromic', 'eur_af_kp_alt', 'fin_1000G_af_kp_alt', 'finngen_af_kp_alt', 'af_diff', 'af_diff_vs_FIN']].round(3).to_string())
print('max |FinnGen - 1000G FIN|', A.af_diff_vs_FIN.max().round(3), 'median', A.af_diff_vs_FIN.median().round(4))
pd.DataFrame(chk).to_csv('step8_fin_order_check.csv', index=False)
