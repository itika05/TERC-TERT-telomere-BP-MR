"""v5: FinnGen R12 allele-harmonisation audit for the 116 primary instruments (KP varId GRCh37, alt = effect allele; FinnGen GRCh38,
beta for alt). Records matched by rsID (FinnGen 'rsids' may list several IDs separated by commas; all are considered)."""
import numpy as np, pandas as pd
F = pd.read_csv('source_data/finngen_R12_I9_HYPTENS_subset_hg38.tsv', sep='\t')
F.columns = ['chr38', 'pos38', 'ref38', 'alt38', 'rsids', 'p_fg', 'b_fg', 'se_fg', 'af_fg']
F = F.assign(rsid=F.rsids.astype(str).str.split(',')).explode('rsid')
V = pd.read_csv('step6_v4_forward_variant_ledger.csv'); V = V[V['selected_r2_0.001']].copy()
comp = dict(A='T', T='A', C='G', G='C')
rows = []
for r in V.itertuples():
    c, p, ref, alt = r.varId.split(':'); f = F[F.rsid == r.rsid]
    pal = {ref, alt} in ({'A', 'T'}, {'C', 'G'})
    if len(f) == 0: rows.append(dict(varId=r.varId, rsid=r.rsid, gene=r.gene, finngen_records=0, status='absent from FinnGen extract', palindromic=pal)); continue
    for g in f.itertuples():
        if (g.ref38, g.alt38) == (ref, alt): st, s = 'same alleles', 1
        elif (g.ref38, g.alt38) == (alt, ref): st, s = 'alleles swapped (beta sign flipped)', -1
        elif (g.ref38, g.alt38) == (comp.get(ref), comp.get(alt)): st, s = 'strand flip (not used)', 0
        else: st, s = 'allele mismatch (not used)', 0
        af = g.af_fg if s == 1 else (1 - g.af_fg if s == -1 else np.nan)
        rows.append(dict(varId=r.varId, rsid=r.rsid, gene=r.gene, finngen_records=len(f), chr38=g.chr38, pos38=g.pos38, fg_ref=g.ref38, fg_alt=g.alt38,
                         status=st, used=s != 0, palindromic=pal, eur_af_kp_alt=r.eur_af_alt, finngen_af_kp_alt=af,
                         af_diff=abs(af - r.eur_af_alt) if s != 0 else np.nan, chr_consistent=str(g.chr38) == c,
                         fg_beta_kp_alt=g.b_fg * s if s else np.nan, fg_se=g.se_fg, fg_p=g.p_fg))
A = pd.DataFrame(rows)
A['palindromic_af_check'] = np.where(A.palindromic & A.used.fillna(False).astype(bool),
    np.where((A.eur_af_kp_alt - .5).abs() > .08, np.where((A.finngen_af_kp_alt - .5) * (A.eur_af_kp_alt - .5) > 0, 'consistent', 'INCONSISTENT'), 'MAF > 0.42: ambiguous'), '')
A.to_csv('step7_finngen_harmonisation_audit.csv', index=False)
print(A.status.value_counts().to_string()); print('palindromic:', A.palindromic_af_check.value_counts().to_dict())
print('max AF diff', A.af_diff.max(), 'n AF diff > 0.1:', (A.af_diff > 0.1).sum()); print(A[A.af_diff > 0.1][['rsid', 'gene', 'eur_af_kp_alt', 'finngen_af_kp_alt']].to_string())
print('chr consistent', A.chr_consistent.value_counts().to_dict())
