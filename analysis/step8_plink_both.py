"""v6: PLINK 1.9 --clump in both MR directions on the same 1000 Genomes EUR genotypes (503 individuals; all founders; no missing calls),
with a per-variant discrepancy ledger. PLINK --clump estimates r2 from maximum-likelihood haplotype frequencies (EM on unphased genotypes);
the custom procedure uses the correlation of the phased 1000 Genomes haplotypes. For each pair that decides a discrepant variant, both r2 values
and the genotype-dosage r2 are reported."""
import numpy as np, pandas as pd, subprocess, base64, re
from ld_ref import load, NHAP
REF = load('ldref/eur_instruments.json')
# v7 fix: records are looked up by position AND alleles (ld_ref.match). v6 used a (chr, pos) key, which at multi-allelic sites picked
# the last record at that position (a different allele), so those variants entered PLINK with another allele's haplotypes.
from ld_ref import match
run = lambda a: subprocess.run(a, shell=True, capture_output=True, text=True)
def bits(i): return np.unpackbits(np.frombuffer(base64.b64decode(REF.iloc[i].bits), np.uint8), bitorder='little')[:NHAP].astype(float)
def plink_set(U, pcol, tag):
    U = U.copy(); t = U.varId.str.split(':', expand=True); U['c'] = t[0]; U['p'] = t[1].astype(int)
    mm = match(U.varId.values, REF).set_index('varId'); U['i'] = mm.loc[U.varId, 'idx'].values
    assert (U['i'] >= 0).all(), 'unmatched variant'
    U = U.drop_duplicates('rsid')
    lines = ['##fileformat=VCFv4.2', '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(f'S{i}' for i in range(NHAP // 2))]
    for r in U.sort_values(['c', 'p'], key=lambda s: s.astype(int) if s.name == 'c' else s).itertuples():
        b = bits(r.i).astype(int); rr = REF.iloc[r.i]
        lines.append(f'{r.c}\t{r.p}\t{r.rsid}\t{rr.ref}\t{rr.alt}\t.\tPASS\t.\tGT\t' + '\t'.join(f'{b[2*k]}|{b[2*k+1]}' for k in range(NHAP // 2)))
    open(f'plink_clump/{tag}.vcf', 'w').write('\n'.join(lines) + '\n')
    U[['rsid', pcol]].rename(columns={'rsid': 'SNP', pcol: 'P'}).to_csv(f'plink_clump/{tag}.assoc', sep='\t', index=False)
    run(f'plink1.9 --vcf plink_clump/{tag}.vcf --double-id --make-bed --out plink_clump/{tag}')
    o = run(f'plink1.9 --bfile plink_clump/{tag} --clump plink_clump/{tag}.assoc --clump-p1 1 --clump-p2 1 --clump-r2 0.001 --clump-kb 10000 --clump-verbose --out plink_clump/{tag}_v')
    txt = open(f'plink_clump/{tag}_v.clumped').read()
    kept = set(re.findall(r'^\s*\d+\s+(\S+)\s+\d+\s+\S+', txt, re.M)) if False else None
    o2 = run(f'plink1.9 --bfile plink_clump/{tag} --clump plink_clump/{tag}.assoc --clump-p1 1 --clump-p2 1 --clump-r2 0.001 --clump-kb 10000 --out plink_clump/{tag}_s')
    kept = set(pd.read_csv(f'plink_clump/{tag}_s.clumped', sep=r'\s+').SNP)
    # verbose: index SNP blocks with member SNPs and their r2
    mem = []; idx = None
    for ln in txt.splitlines():
        f = ln.split()
        if len(f) >= 2 and f[0] == '(INDEX)': idx = f[1]; continue
        if idx and len(f) >= 6 and f[0].startswith('rs') and f[0] != idx:
            try: mem.append(dict(index=idx, member=f[0], plink_r2=float(f[2])))
            except ValueError: pass
        if ln.startswith('-----'): idx = None
    return U, kept, pd.DataFrame(mem, columns=['index', 'member', 'plink_r2'])
def phased_r2(U, a, b):
    ia = U.loc[U.rsid == a, 'i'].iloc[0]; ib = U.loc[U.rsid == b, 'i'].iloc[0]; ha, hb = bits(ia), bits(ib)
    ga, gb = ha[0::2] + ha[1::2], hb[0::2] + hb[1::2]
    return np.corrcoef(ha, hb)[0, 1] ** 2, np.corrcoef(ga, gb)[0, 1] ** 2
out, summ = [], []
V = pd.read_csv('step6_v4_forward_variant_ledger.csv'); FU = V[V.usable & ~V.mhc]
jobs = [('forward', FU, 'L_codd_p', set(FU.loc[FU['selected_r2_0.001'], 'rsid']))]
RV = pd.read_csv('step3v3_reverse_variant_ledgers_all.csv')
for ex, pc in [('SBP', 'S_mvp_p'), ('DBP', 'D_mvp_p'), ('Hypertension', 'H_mvp_p')]:
    u = RV[(RV.exposure == ex) & RV.usable]; jobs.append((ex, u, pc, set(u.loc[u['selected_r2_0.001'], 'rsid'])))
for tag, U, pc, custom in jobs:
    U2, kept, mem = plink_set(U, pc, 'ref_' + tag.lower() if tag != 'forward' else 'fwd')
    both = custom & kept
    summ.append(dict(direction=tag, candidates=len(U), custom_kept=len(custom), plink_kept=len(kept), shared=len(both), custom_only=len(custom - kept), plink_only=len(kept - custom),
                     jaccard=len(both) / len(custom | kept)))
    # discrepancy ledger: custom-only variants were removed by PLINK (find the PLINK index that absorbed them); PLINK-only variants were removed by custom
    CL = pd.read_csv('step6_v4_forward_clump_log.csv') if tag == 'forward' else pd.read_csv(f'step3v3_reverse_{tag}_clump_log.csv')
    rsmap = dict(zip(U2.varId, U2.rsid))
    for v in sorted(custom - kept):
        m = mem[mem.member == v]
        for r in m.itertuples():
            ph, dos = phased_r2(U2, r.index, v)
            out.append(dict(direction=tag, variant=v, status='custom only (PLINK removed it)', decided_by=r.index, plink_ml_r2=r.plink_r2, phased_haplotype_r2=ph, genotype_dosage_r2=dos))
        if m.empty: out.append(dict(direction=tag, variant=v, status='custom only (PLINK removed it)', decided_by='not found in verbose output'))
    for v in sorted(kept - custom):
        c = CL[CL.varId.map(rsmap) == v]
        if len(c):
            by = rsmap.get(c.by.iloc[0], c.by.iloc[0]); ph, dos = phased_r2(U2, by, v) if by in set(U2.rsid) else (np.nan, np.nan)
            out.append(dict(direction=tag, variant=v, status='PLINK only (custom removed it)', decided_by=by, plink_ml_r2=np.nan, phased_haplotype_r2=ph, genotype_dosage_r2=dos, custom_reason=c.reason.iloc[0]))
    pd.Series(sorted(kept)).to_csv(f'plink_clump/{tag}_plink_kept.txt', index=False, header=False)
S = pd.DataFrame(summ); S.to_csv('step8_plink_both_directions_summary.csv', index=False); print(S.to_string())
L = pd.DataFrame(out).rename(columns={'plink_ml_r2': 'plink_clump_r2'})
# r2 of each deciding pair from PLINK itself: --r2 (allele-count correlation) and --r2 dprime (maximum-likelihood haplotype frequencies,
# the estimator used by --clump), on the same PLINK filesets
tagof = {'forward': 'fwd', 'SBP': 'ref_sbp', 'DBP': 'ref_dbp', 'Hypertension': 'ref_hypertension'}
def pair_r2(tg, a, b, extra):
    open('plink_clump/pair.txt', 'w').write(f'{a}\n{b}\n')
    run(f'plink1.9 --bfile plink_clump/{tg} --extract plink_clump/pair.txt --r2 {extra} --ld-window 99999 --ld-window-kb 20000 --ld-window-r2 0 --out plink_clump/pair')
    try: return float(pd.read_csv('plink_clump/pair.ld', sep=r'\s+').R2.iloc[0])
    except Exception: return np.nan
ac, ml = [], []
for r in L.itertuples():
    tg = tagof[r.direction]
    if not str(r.decided_by).startswith('rs'): ac.append(np.nan); ml.append(np.nan); continue
    ac.append(pair_r2(tg, r.variant, r.decided_by, '')); ml.append(pair_r2(tg, r.variant, r.decided_by, 'dprime'))
L['plink_r2_command_allele_count'] = ac
L['plink_r2_dprime_ML'] = ml
# selection mechanism: is the deciding variant itself kept by both procedures? If the blocker was retained by one procedure but removed by
# the other, the difference propagates through the greedy selection chain rather than through this pair's own threshold crossing.
kept_sets = {}
for tag, U, pc, custom in jobs:
    kept_sets[tag] = (custom, set(open(f'plink_clump/{tag}_plink_kept.txt').read().split()))
def mech(r):
    if not str(r.decided_by).startswith('rs'): return 'deciding variant not identified'
    cu, pk = kept_sets[r.direction]
    in_cu, in_pk = r.decided_by in cu, r.decided_by in pk
    thr = 0.001
    ml_ = r.plink_r2_dprime_ML; ph = r.phased_haplotype_r2
    if r.status.startswith('custom only'):
        if not in_cu: return f'selection chain: blocker {r.decided_by} kept by PLINK but removed by the custom procedure, so it could not remove {r.variant} there'
        return f'threshold crossing: pair r2 {ml_:.5f} (PLINK ML) >= {thr} > {ph:.5f} (phased haplotype)' if (ml_ >= thr and ph < thr) else 'other (see r2 columns)'
    if not in_pk: return f'selection chain: blocker {r.decided_by} kept by the custom procedure but removed by PLINK, so it could not remove {r.variant} there'
    return f'threshold crossing: pair r2 {ph:.5f} (phased haplotype) >= {thr} > {ml_:.5f} (PLINK ML)' if (ph >= thr and ml_ < thr) else 'other (see r2 columns)'
L['mechanism'] = L.apply(mech, axis=1)
L.to_csv('step8_plink_discrepancy_ledger.csv', index=False)
print(L.groupby(['direction', 'status']).size()); print(L[L.direction == 'forward'].to_string())
