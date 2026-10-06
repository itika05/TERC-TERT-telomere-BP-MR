"""v8 release checks and manifest (v7 checks plus LD-panel allele-key checks, qhet objective and bootstrap checks, embedded-figure identity). Run after all analyses (before the documents are built): machine-readable PASS/FAIL checks that the
publication tables belong to one run and reconstruct from their named inputs, plus a SHA-256 manifest and a release ID derived from the hashes
of all inputs, scripts and result files. With --render (after the documents are built) it adds checks that key values in the rendered
manuscript equal the result files, and hashes the rendered outputs. The build fails if any check fails."""
import sys, os, glob, json, hashlib, re
import numpy as np, pandas as pd
from scipy import stats
os.makedirs('release_v8', exist_ok=True)
C = []
def add(check, ok, detail=''): C.append(dict(check=check, status='PASS' if ok else 'FAIL', detail=detail))
def ivw(bx, by, sy):
    w = 1 / sy ** 2; b = (w * bx * by).sum() / (w * bx ** 2).sum(); k = len(bx)
    rss = (w * (by - b * bx) ** 2).sum(); se = np.sqrt(max(1.0, rss / (k - 1)) / (w * bx ** 2).sum()); p = 2 * stats.t.sf(abs(b / se), k - 1)
    return b, se, p, rss, k
R = pd.read_csv('step8_mr_results_v6.csv')
key = R.direction.astype(str) + '|' + R.code + '|' + R.analysis + '|' + R.method + '|' + R.source.astype(str)
add('S11 analysis keys unique (direction|code|analysis|method|source)', not key.duplicated().any(), f'{len(R)} rows')
O = pd.read_csv('step9_mr_output_consistency_S46.csv'); add('S46: every S11 row reproduces P and CI from b and SE', (O.status == 'PASS').all() and len(O) == len(R), f'{(O.status == "PASS").sum()}/{len(O)} rows')
T = pd.read_csv('step9_reverse_reconstruction_test.csv'); add('S72: every reverse IVW result reconstructs from S13 rows (k, b, SE, P)', (T.status == 'PASS').all(), f'{len(T)} results; max rel diff {T.max_rel_diff.max():.1e}')
# forward IVW from S12 rows (primary, FinnGen) and PLINK forward inputs
S12 = pd.read_csv('step6_v4_harmonised_inputs_all.csv'); worst = 0; n = 0; okk = True
for code in ['S_mvpEU', 'D_mvpEU', 'H_mvpEU', 'H_fg12']:
    d = S12[(S12.code == code) & (S12.file == f'forward_selected_r2_0.001_{code}.csv')]
    b, se, p, rss, k = ivw(d.bx.values, d.by.values, d.sy.values); st = R[(R.code == code) & (R.analysis == 'r2<0.001 (primary)') & (R.method == 'IVW (MRE, t)')].iloc[0]
    rel = max(abs(b - st.b) / abs(st.b), abs(se - st.se) / st.se); worst = max(worst, rel); n += 1; okk &= (k == st.nsnp) and rel < 1e-8
    d = pd.read_csv(f'mr_v6_inputs/forward_plink_{code}.csv'); b, se, p, rss, k = ivw(d.bx.values, d.by.values, d.sy.values)
    st = R[(R.code == code) & (R.analysis == 'PLINK 1.9 clumping') & (R.method == 'IVW (MRE, t)')].iloc[0]
    rel = max(abs(b - st.b) / abs(st.b), abs(se - st.se) / st.se); worst = max(worst, rel); n += 1; okk &= (k == st.nsnp) and rel < 1e-8
add('Forward IVW (primary from S12 rows; PLINK sets) reconstructs S11', okk, f'{n} results; max rel diff {worst:.1e}')
# MR-PRESSO corrected estimates at B = 1e6 reconstruct from S12 minus the listed outliers
CR = pd.read_csv('step8_presso_convergence_runs.csv'); okp = True; det = []
for code in ['S_mvpEU', 'D_mvpEU', 'H_mvpEU']:
    r6 = CR[(CR.code == code) & (CR.B == 1000000)]; sets = r6.outliers.fillna('').map(lambda x: frozenset(x.split(';')) - {''})
    same = sets.nunique() == 1; okp &= same
    d = S12[(S12.code == code) & (S12.file == f'forward_selected_r2_0.001_{code}.csv')]; keep = ~d.rsid.isin(sets.iloc[0])
    dd = d[keep]; w = 1 / dd.sy ** 2; b = (w * dd.bx * dd.by).sum() / (w * dd.bx ** 2).sum()
    st = R[(R.code == code) & (R.method == 'MR-PRESSO (outlier-corrected)')].iloc[0]; okp &= abs(b - st.b) < 1e-10 and int(keep.sum()) == int(st.nsnp)
    det.append(f'{code}: {len(sets.iloc[0])} outliers, k {int(keep.sum())}')
add('MR-PRESSO B = 1e6: identical outlier sets over 5 seeds; corrected b reconstructs from S12', okp, '; '.join(det))
E = pd.read_csv('step9_presso_equivalence.csv'); E = E[E.check != 'summary']
add('S70: Python MR-PRESSO agrees with MRPRESSO 1.0 on shared simulated data within stated tolerances (all intermediates)', len(E) > 0 and (E.status == 'PASS').all() and E.code.nunique() == 3, f'{(E.status == "PASS").sum()}/{len(E)} checks, {E.code.nunique()} outcomes')
# MVMR: weighted regression on S20 minus S39 reproduces S21 (all three exposures)
M = pd.read_csv('mr_v4_inputs/mvmr_instruments.csv'); fb = set(pd.read_csv('step7_mvmr_fallback_variants.csv').varId); M = M[~M.varId.isin(fb)]
E21 = pd.read_csv('step9_mvmr_estimates_primary.csv'); okm = len(M) == 318; worst = 0
for oc in ['S_mvpEU', 'D_mvpEU', 'H_mvpEU']:
    X = M[['LTL_b', 'BMI_b', 'LYM_b']].values; y = M[oc + '_b'].values; w = 1 / M[oc + '_se'].values ** 2
    beta = np.linalg.solve(X.T @ (X * w[:, None]), X.T @ (w * y))
    for j, ex in enumerate(['LTL', 'BMI', 'LYM']):
        st = E21[(E21.outcome == oc) & (E21.exposure == ex) & (E21.method == 'MV-IVW (random, t)')].iloc[0]
        rel = abs(beta[j] - st.b) / abs(st.b); worst = max(worst, rel); okm &= rel < 1e-8 and int(st.nsnp) == 318
add('S21: MV-IVW coefficients (LTL, BMI, LYM) reconstruct from S20 minus S39 (318 variants)', okm, f'max rel diff {worst:.1e}')
Q = pd.read_csv('step9_mvmr_qhet.csv'); g = Q.groupby(['outcome', 'rho'])
okq = (g.size() == 3).all() and (g.exposure.apply(lambda s: list(s) == ['LTL', 'BMI', 'LYM'])).all() and (g.qhet_estimate.nunique() == 3).all()
V = pd.read_csv('step9_qhet_vectorised_check.csv'); okq &= V.abs_diff.max() < 1e-10
add('S55: qhet_mvmr exported one row per exposure (LTL, BMI, LYM) per outcome and rho; three distinct values; vectorised objective equals reference', bool(okq), f'{len(g)} outcome x rho groups')
# variant accounting for reverse PLINK sets before/after outcome matching
RL = pd.read_csv('step8_reverse_source_ledger.csv'); okv = True; det = []
for r in RL.itertuples():
    an = 'r2<0.001 (primary)' if r.selection == 'custom' else 'PLINK 1.9 clumping'
    st = R[(R.code == f'{r.exposure}->L_codd') & (R.analysis == an) & (R.method == 'IVW (MRE, t)')].iloc[0]
    okv &= (r.instruments - r.absent - r.allele_mismatch) == st.nsnp; det.append(f'{r.exposure}/{r.selection}: {r.instruments} selected, {r.absent} without source record, k {int(st.nsnp)}')
add('Reverse variant accounting: selected - unmatched = k used', okv, '; '.join(det))
L52 = pd.read_csv('step8_plink_discrepancy_ledger.csv')
okl = L52.plink_r2_dprime_ML.notna().all() and (L52.mechanism.str.startswith('threshold') | L52.mechanism.str.startswith('selection chain')).all()
c = L52[L52.plink_clump_r2.notna()]; okl &= ((c.plink_clump_r2 - c.plink_r2_dprime_ML).abs() / c.plink_clump_r2).max() < 0.005
add('S52: every discrepant variant has a PLINK --r2 dprime value and an explained mechanism; dprime r2 matches clump r2 (<0.5%)', bool(okl), f'{len(L52)} rows')
L38 = pd.read_csv('step9_se_fallback_ledger_S38.csv')
add('S38: no reverse fallback record is marked as used; MVMR fallback variants all excluded', (L38[L38.direction == 'reverse'].status_current_release.str.startswith('superseded')).all() and set(L38[L38.direction == 'MVMR'].varId).isdisjoint(set(M.varId)), f'{len(L38)} records')

# ---- v8 checks
# LD-panel allele key: every variant written to a PLINK fileset carries its own alleles (either orientation), and no script that reads the
# LD panel uses a position-only key
okld = True; det = []
for vcf in sorted(glob.glob('plink_clump/*.vcf')):
    if os.path.basename(vcf) in ('pair.vcf',): continue
    bad = n = 0
    ids = {}
    for tag, U in [('fwd', pd.read_csv('step6_v4_forward_variant_ledger.csv'))] + [(f'ref_{e.lower()}', pd.read_csv('step3v3_reverse_variant_ledgers_all.csv').query('exposure == @e')) for e in ['SBP', 'DBP', 'Hypertension']]:
        if os.path.basename(vcf).startswith(tag): ids = dict(zip(U.rsid, U.varId))
    if not ids: continue
    for line in open(vcf):
        if line.startswith('#'): continue
        f = line.split('\t', 5); n += 1; v = ids.get(f[2])
        if v is None: bad += 1; continue
        a0, a1 = v.split(':')[2:4]
        if {f[3], f[4]} != {a0, a1}: bad += 1
    okld &= bad == 0; det.append(f'{os.path.basename(vcf)}: {n} variants, {bad} allele mismatches')
posonly = [f for f in glob.glob('*.py') if 'ld_ref' in open(f).read() and re.search(r"key = \{\(r\.chr, r\.pos\): i", open(f).read()) and not f.startswith(('step9_release', 'step10_release'))]
add('LD panel: PLINK filesets carry each variant\'s own alleles; no current LD script uses a position-only key', okld and not posonly, '; '.join(det) + (f'; position-only key in {posonly}' if posonly else ''))
# every LD-panel lookup through ld_ref.match (primary clumping in both directions, MVMR, GLS set, colocalization windows) returns a record
# whose alleles are the variant's own alleles; palindromic variants are flagged by match() and handled by the documented rule
from ld_ref import load as _load, match as _match
okm2 = True; det = []
REFI = _load('ldref/eur_instruments.json')
sets_ = [('forward candidates', pd.read_csv('step6_v4_forward_variant_ledger.csv').varId), ('reverse candidates', pd.read_csv('step3v3_reverse_variant_ledgers_all.csv').varId.drop_duplicates()),
         ('MVMR instruments', pd.read_csv('mr_v4_inputs/mvmr_instruments.csv').varId)]
REFM = pd.concat([_load('ldref/eur_instruments.json'), _load('ldref/eur_mvmr_candidates.json')]).drop_duplicates(['chr', 'pos', 'ref', 'alt']).reset_index(drop=True)   # as in step6_mvmr_prep.py
for lab, vids in sets_:
    P_ = REFM if lab.startswith('MVMR') else REFI
    mm = _match(list(vids), P_); used = mm[mm.idx >= 0]
    badn = sum(1 for r in used.itertuples() if {P_.iloc[r.idx].ref, P_.iloc[r.idx].alt} != set(r.varId.split(':')[2:4]))
    okm2 &= badn == 0; det.append(f'{lab}: {len(used)} matched records, {badn} allele mismatches, {int(used.palindromic.sum())} palindromic flagged')
for loc in ['TERC', 'TERT']:
    REFR = _load(f'ldref/eur_region_{loc}.json'); vids = pd.concat([pd.read_csv(f) for f in glob.glob(f'coloc_v4_inputs/{loc}_250kb_*.csv') if not f.endswith(('_LD.csv', '_LDdos.csv'))]).varId.drop_duplicates()
    mm = _match(list(vids), REFR); used = mm[mm.idx >= 0]
    badn = sum(1 for r in used.itertuples() if {REFR.iloc[r.idx].ref, REFR.iloc[r.idx].alt} != set(r.varId.split(':')[2:4]))
    okm2 &= badn == 0; det.append(f'{loc} colocalization window: {len(used)} matched records, {badn} allele mismatches')
add('LD panel: every allele-aware lookup (clumping candidates, MVMR, colocalization windows) returns the variant\'s own alleles', okm2, '; '.join(det))
T74 = pd.read_csv('step10_ld_lookup_correction_S74.csv')
add('S74: variants affected by the former position-only lookup are exactly rs113695388 and rs147632758, and neither is in a current PLINK set', set(T74[~T74.position_only_record_correct].rsid) == {'rs113695388', 'rs147632758'} and not T74[~T74.position_only_record_correct].plink_kept_current.any(), f'{len(T74)} multi-allelic candidates')
QO = pd.read_csv('step10_qhet_objective_test.csv')
add('qhet objective: vectorised objective equals the reference PL2_MVMR (extracted from qhet_mvmr) at random coefficient vectors and tau2', QO.rel_diff.max() < 1e-10, f'{len(QO)} points; max rel diff {QO.rel_diff.max():.1e}')
okb = True; det = []
for oc in ['S_mvpEU', 'D_mvpEU', 'H_mvpEU']:
    Rb = pd.read_csv(f'qhet_bootstrap/{oc}_rho0_replicates.csv'); q = Q[(Q.outcome == oc) & (Q.rho == 0)]
    ok_ = Rb.fit_ok == 1
    for ex in ['LTL', 'BMI', 'LYM']:
        r = q[q.exposure == ex].iloc[0]
        okb &= abs(np.quantile(Rb.loc[ok_, ex], 0.025) - r.pct_lo) < 1e-12 and abs(np.quantile(Rb.loc[ok_, ex], 0.975) - r.pct_hi) < 1e-12 and len(Rb) == r.bootstrap_B
    det.append(f'{oc}: {len(Rb)} replicates, {int((~ok_).sum())} failed, seed {int(Rb.seed.iloc[0])}')
add('S55 bootstrap intervals reconstruct from the saved replicates (qhet_bootstrap/)', bool(okb), '; '.join(det))
# no stale sources in the publication workbook builder
bs = open('build_supplement_v8.py').read()
stale = [f for f in ['step6_v4_reverse_inputs_all.csv', 'step7_se_fallback_ledger.csv', 'step7_mr_output_consistency.csv', 'step8_mvmr_v6_qhet.csv', 'step8_mvmr_v6_estimates.csv'] if f"'{f}'" in bs]
add('Workbook uses no superseded source for a current table', not stale, ', '.join(stale) or 'none')

def sha(p): return hashlib.sha256(open(p, 'rb').read()).hexdigest()
pats = ['source_data/*', 'ldref/*.json', 'mr_v4_inputs/*.csv', 'mr_v6_inputs/*.csv', 'ref_MVMR/*.R', '*.py', '*.R', 'run_all.sh', 'step*.csv', 'presso_equiv/*', 'plink_clump/*.log', 'plink_clump/*.clumped', 'plink_clump/*_kept.txt', 'release_v8/environment.txt', 'step10_render_pdfs.sh', 'INDEPENDENT_RERUN_PROTOCOL.md', 'qhet_bootstrap/*', 'audit_v6_plink/*', 'figures/proof_report.csv']
files = sorted({f for p in pats for f in glob.glob(p) if os.path.isfile(f) and not f.startswith('scratch_')})
render = '--render' in sys.argv
man = pd.DataFrame([dict(path=f, sha256=sha(f), bytes=os.path.getsize(f)) for f in files])
rid = 'v8-' + hashlib.sha256(''.join(man.path + man.sha256).encode()).hexdigest()[:12]
if render:
    rid = json.load(open('release_v8/release_summary.json'))['release_id']
    import docx
    txt = '\n'.join(p.text for p in docx.Document('out/Manuscript_TERC_TERT_BP_v8_2026-09-30.docx').paragraphs)
    for t in docx.Document('out/Manuscript_TERC_TERT_BP_v8_2026-09-30.docx').tables: txt += '\n'.join(c.text for row in t.rows for c in row.cells)
    g_ = lambda code, an, m: R[(R.code == code) & (R.analysis == an) & (R.method == m)].iloc[0]
    want = []
    for code in ['S_mvpEU', 'D_mvpEU']:
        r = g_(code, 'r2<0.001 (primary)', 'IVW (MRE, t)'); want.append(f'{r.b:.3f}'); want.append(f'{r.lo:.3f} to {r.hi:.3f}')
    for code in ['SBP->L_codd', 'DBP->L_codd', 'Hypertension->L_codd']:
        r = g_(code, 'r2<0.001 (primary)', 'IVW (MRE, t)'); want.append(f'{r.b:.4f}'.replace('-', '−'))
    r = g_('H_mvpEU', 'r2<0.001 (primary)', 'IVW (MRE, t)'); want.append(f'OR {np.exp(r.b):.2f}')
    mv = E21[(E21.outcome == 'S_mvpEU') & (E21.exposure == 'LTL') & (E21.method == 'MV-IVW (random, t)')].iloc[0]; want.append(f'{mv.b:.3f}')
    for oc in ['S_mvpEU']:
        q = Q[(Q.outcome == oc) & (Q.exposure == 'LTL') & (Q.rho == 0)].iloc[0]; want.append(f'{q.pct_lo:.3f} to {q.pct_hi:.3f}')
    for code in ['S_mvpEU', 'D_mvpEU', 'H_mvpEU']:
        want.append(f'{int(CR[(CR.code == code) & (CR.B == 1000000)].n_outliers.iloc[0])}')
    miss = [w for w in want if w not in txt.replace('-', '−') and w not in txt]
    add('Rendered manuscript contains the key values of the result files (primary, reverse, MVMR, qhet, MR-PRESSO)', not miss, f'{len(want) - len(miss)}/{len(want)} values found' + (f'; missing {miss}' if miss else ''))
    stale_txt = [s for s in ['validated re-implementation', 'prespecified stability', 'explains the TERC', 'all 330 variants', 'Keaton2024 2024', 'v5', 'v6 ', 'reproduced MRPRESSO 1.0 exactly', 'remove Monte Carlo instability', 'MR exposure (all forward', 'with unresolved multi-signal sensitivity; regional'] if s in txt]
    add('Rendered manuscript contains no withdrawn wording or version labels', not stale_txt, ', '.join(stale_txt) or 'none')
    import zipfile
    fig_h = {sha(f): f for f in glob.glob('figures/Fig*.png')}
    det = []; oke = True
    for docf in ['out/Manuscript_TERC_TERT_BP_v8_2026-09-30.docx', 'out/Supplementary_Information_v8_2026-09-30.docx']:
        z = zipfile.ZipFile(docf); media = [n for n in z.namelist() if n.startswith('word/media/')]
        hs = [hashlib.sha256(z.read(n)).hexdigest() for n in media]; found = [fig_h.get(h) for h in hs]
        oke &= all(found) and len(media) > 0; det.append(f'{os.path.basename(docf)}: {sum(1 for x in found if x)}/{len(media)} images identical to release figures')
    add('Figures embedded in the manuscript and Supplementary Information are byte-identical to the release figure files', oke, '; '.join(det))
    PR_ = pd.read_csv('figures/proof_report.csv'); main_ = PR_[PR_.figure.str.match(r'^Fig\d')]
    add('Main figures at 180 mm print width: no text below 6 pt', len(main_) == 4 and (main_.min_effective_pt.astype(float) >= 6.0 - 1e-9).all(), '; '.join(f'{r.figure}: min {float(r.min_effective_pt):.2f} pt, scale {float(r.scale_at_180mm):.3f}' for r in main_.itertuples()))
    outs = sorted(glob.glob('out/*v8*') + glob.glob('figures/Fig*.p*'))
    man = pd.concat([man, pd.DataFrame([dict(path=f, sha256=sha(f), bytes=os.path.getsize(f)) for f in outs])], ignore_index=True)
    prev = pd.read_csv('release_v8/build_checks.csv'); prev = prev[prev.phase == 'pre-build']
    CC = pd.concat([prev, pd.DataFrame(C).assign(release_id=rid, phase='post-build')], ignore_index=True)
else:
    CC = pd.DataFrame(C).assign(release_id=rid, phase='pre-build')
CC = CC[['release_id', 'phase', 'check', 'status', 'detail']]; CC.to_csv('release_v8/build_checks.csv', index=False); man.to_csv('release_v8/release_manifest.csv', index=False)
pre = CC[CC.phase == 'pre-build']
json.dump(dict(release_id=rid, n_checks=int(len(pre)), n_checks_pass=int((pre.status == 'PASS').sum()), n_post=int((CC.phase == 'post-build').sum()),
               n_post_pass=int(((CC.phase == 'post-build') & (CC.status == 'PASS')).sum()), n_files=int(len(man))), open('release_v8/release_summary.json', 'w'), indent=1)
print(CC.to_string()); print('release', rid, 'files', len(man))
if (CC.status != 'PASS').any(): sys.exit('RELEASE CHECKS FAILED')
