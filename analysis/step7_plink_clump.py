"""v5: compare the custom clumping (step3_mr_v3.clump, r2 from phased haplotype correlation) with PLINK 1.9 --clump on the same
1000 Genomes phase 3 EUR genotypes (503 individuals) and the same 174 usable non-MHC autosomal candidates, index P = portal Codd P."""
import numpy as np, pandas as pd, subprocess, base64
from ld_ref import load, NHAP
REF = load('ldref/eur_instruments.json')
V = pd.read_csv('step6_v4_forward_variant_ledger.csv')
U = V[V.usable & ~V.mhc].copy(); print('usable non-MHC', len(U))
key = {(r.chr, r.pos, r.ref, r.alt): i for i, r in enumerate(REF.itertuples())}
lines = ['##fileformat=VCFv4.2', '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(f'S{i}' for i in range(NHAP // 2))]
rows = []
for r in U.itertuples():
    c, p, a0, a1 = r.varId.split(':'); p = int(p)
    i = key.get((c, p, a0, a1), key.get((c, p, a1, a0)))
    rr = REF.iloc[i]; b = np.unpackbits(np.frombuffer(base64.b64decode(rr.bits), np.uint8), bitorder='little')[:NHAP]
    gt = '\t'.join(f'{b[2*k]}|{b[2*k+1]}' for k in range(NHAP // 2))
    rows.append((int(c), p, f'{c}\t{p}\t{r.rsid}\t{rr.ref}\t{rr.alt}\t.\tPASS\t.\tGT\t{gt}'))
lines += [x[2] for x in sorted(rows)]
open('plink_clump/cand.vcf', 'w').write('\n'.join(lines) + '\n')
U[['rsid', 'L_codd_p']].rename(columns={'rsid': 'SNP', 'L_codd_p': 'P'}).to_csv('plink_clump/cand.assoc', sep='\t', index=False)
run = lambda a: subprocess.run(a, shell=True, capture_output=True, text=True)
print(run('plink1.9 --vcf plink_clump/cand.vcf --double-id --make-bed --out plink_clump/cand').stdout[-300:])
out = {}
for r2 in ['0.001', '0.1']:
    o = run(f'plink1.9 --bfile plink_clump/cand --clump plink_clump/cand.assoc --clump-p1 1 --clump-p2 1 --clump-r2 {r2} --clump-kb 10000 --out plink_clump/clump_{r2}')
    C = pd.read_csv(f'plink_clump/clump_{r2}.clumped', sep=r'\s+'); out[r2] = set(C.SNP)
res = []
for r2, col in [('0.001', 'selected_r2_0.001'), ('0.1', 'selected_r2_0.1')]:
    mine = set(U.loc[U[col], 'rsid']); pl = out[r2]
    res.append(dict(r2_threshold=r2, window_kb=10000, n_candidates=len(U), custom_kept=len(mine), plink_kept=len(pl), both=len(mine & pl),
                    custom_only=';'.join(sorted(mine - pl)), plink_only=';'.join(sorted(pl - mine))))
R = pd.DataFrame(res); R.to_csv('step7_plink_clump_comparison.csv', index=False); print(R.to_string())

# IVW (multiplicative random effects, t(k-1)) on the PLINK-selected r2 < 0.001 set, using the harmonised r2 < 0.1 inputs (superset)
from scipy import stats
def ivw(d):
    w = 1 / d.sy**2; b = (w * d.bx * d.by).sum() / (w * d.bx**2).sum(); k = len(d)
    phi = max(1, (((d.by - b * d.bx) / d.sy)**2).sum() / (k - 1)); se = np.sqrt(phi / (w * d.bx**2).sum()); t = stats.t.ppf(0.975, k - 1)
    return dict(k=k, b=b, se=se, lo=b - t * se, hi=b + t * se, p=2 * stats.t.sf(abs(b / se), k - 1))
iv = []
for code in ['S_mvpEU', 'D_mvpEU', 'H_mvpEU', 'H_fg12']:
    sup = pd.read_csv(f'mr_v4_inputs/forward_selected_r2_0.1_{code}.csv'); cus = pd.read_csv(f'mr_v4_inputs/forward_selected_r2_0.001_{code}.csv')
    iv.append(dict(code=code, instrument='custom clumping (primary)', **ivw(cus)))
    iv.append(dict(code=code, instrument='PLINK 1.9 --clump', **ivw(sup[sup.rsid.isin(out['0.001'])])))
I = pd.DataFrame(iv); I.to_csv('step7_plink_clump_ivw.csv', index=False); print(I.to_string())
