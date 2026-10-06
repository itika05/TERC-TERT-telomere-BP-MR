"""v8 (Supplementary Table S74): correction trail for the LD-panel lookup. For every PLINK candidate (forward and the three reverse
directions) at a position with more than one 1000 Genomes record, lists the panel records at that position, the record a position-only
key selects (the v6 PLINK script: last record at the position) and the record selected by the allele-aware match (ld_ref.match: position
AND both alleles), and whether the variant was kept by PLINK in v6 (audit_v6_plink/) and in the current release, and used in an MR fit."""
import pandas as pd
from ld_ref import load, match
REF = load('ldref/eur_instruments.json')
pos = {}
for i, r in enumerate(REF.itertuples()): pos.setdefault((r.chr, r.pos), []).append(i)
V = pd.read_csv('step6_v4_forward_variant_ledger.csv'); FU = V[V.usable & ~V.mhc]
RV = pd.read_csv('step3v3_reverse_variant_ledgers_all.csv')
S13 = pd.read_csv('step9_reverse_inputs_S13.csv'); FP = pd.read_csv('mr_v6_inputs/forward_plink_S_mvpEU.csv')
sets = [('forward', FU)] + [(ex, RV[(RV.exposure == ex) & RV.usable]) for ex in ['SBP', 'DBP', 'Hypertension']]
rows = []
for tag, U in sets:
    U = U.drop_duplicates('rsid'); m = match(U.varId.values, REF).set_index('varId')
    old = set(open(f'audit_v6_plink/{tag}_plink_kept.txt').read().split()); new = set(open(f'plink_clump/{tag}_plink_kept.txt').read().split())
    for r in U.itertuples():
        c, p = r.varId.split(':')[0], int(r.varId.split(':')[1]); recs = pos.get((c, p), [])
        if len(recs) < 2: continue
        last = REF.iloc[recs[-1]]; aa = m.loc[r.varId, 'idx']; aw = REF.iloc[aa] if aa >= 0 else None
        used = (r.rsid in set(FP.rsid)) if tag == 'forward' else bool(((S13.code == f'{tag}->L_codd') & (S13.rsid == r.rsid) & (S13.in_plink == True)).any())
        rows.append(dict(direction=tag, varId=r.varId, rsid=r.rsid, panel_records_at_position='; '.join(f'{REF.iloc[i].ref}>{REF.iloc[i].alt}' for i in recs),
                         position_only_record=f'{last.ref}>{last.alt}', allele_aware_record=(f'{aw.ref}>{aw.alt}' if aw is not None else 'no matching record'),
                         position_only_record_correct=(aw is not None and recs[-1] == aa), plink_kept_v6=r.rsid in old, plink_kept_current=r.rsid in new,
                         in_current_plink_mr_fit=used))
T = pd.DataFrame(rows); T.to_csv('step10_ld_lookup_correction_S74.csv', index=False)
print(T.to_string()); print('wrong record under position-only key:', int((~T.position_only_record_correct).sum()))
