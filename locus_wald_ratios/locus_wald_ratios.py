"""Locus-specific Wald ratios at TERC and TERT compared with genome-wide IVW (audit v17, major point 4-5).
Inputs: per-variant MVP and Codd values from step6_codd_published_instrument_lookup.csv, step11_steiger_reverse_removed.csv
and FinnGen values from Table S89 (v17); IVW estimates from step8_mr_results_v6.csv."""
import numpy as np, pandas as pd
from scipy.stats import norm
def se_from_p(b, p): return abs(b) / norm.isf(p / 2)
# (locus, variant, role, outcome, bx, sx, by, sy, comparator code)
rows = [
 ('TERC','rs2293607','LTL lead','SBP (MVP)',      0.0944803,0.00233425, 0.0332,     se_from_p(0.0332,4.422e-30), 'S_mvpEU'),
 ('TERC','rs2293607','LTL lead','DBP (MVP)',      0.0944803,0.00233425, 0.009243,   se_from_p(0.009243,0.001143),'D_mvpEU'),
 ('TERC','rs2293607','LTL lead','Hypertension (MVP)',0.0944803,0.00233425,0.029558802,se_from_p(0.029558802,1.442e-07),'H_mvpEU'),
 ('TERC','rs2293607','LTL lead','Hypertension (FinnGen)',-0.0945,0.0023, -0.0211,0.0057,'H_fg12'),
 ('TERC','rs12638862','SBP lead','SBP (MVP)',    -0.0861226,0.00227875,-0.03307,0.0028525744979717,'S_mvpEU'),
 ('TERT','rs7705526','LTL lead','SBP (MVP)',     -0.0776022,0.00216124,-0.01868,  se_from_p(0.01868,1.162e-11),'S_mvpEU'),
 ('TERT','rs7705526','LTL lead','DBP (MVP)',     -0.0776022,0.00216124,-0.01293,  se_from_p(0.01293,1.459e-06),'D_mvpEU'),
 ('TERT','rs7705526','LTL lead','Hypertension (MVP)',-0.0776022,0.00216124,-0.014504686,0.0048110684828802,'H_mvpEU'),
 ('TERT','rs7705526','LTL lead','Hypertension (FinnGen)',0.0776,0.0022,0.0294,0.0054,'H_fg12'),
]
R = pd.read_csv('../analysis/step8_mr_results_v6.csv')
def ivw(code, analysis):
    x = R[(R.code==code)&(R.analysis==analysis)&(R.method=='IVW (MRE, t)')]
    return (x.b.iloc[0], x.se.iloc[0]) if len(x) else (np.nan, np.nan)
out=[]
for loc,v,role,oc,bx,sx,by,sy,code in rows:
    r = by/bx; se1 = sy/abs(bx); se2 = np.sqrt(sy**2/bx**2 + by**2*sx**2/bx**4)
    gw_b, gw_se = ivw(code,'r2<0.001 (primary)'); ex_b, ex_se = ivw(code,'r2<0.001, excluding TERC/TERT regions')
    comp_b, comp_se, comp = (ex_b, ex_se, 'IVW without TERC/TERT regions') if not np.isnan(ex_b) else (gw_b, gw_se, 'genome-wide IVW')
    z = (r-comp_b)/np.sqrt(se2**2+comp_se**2); p = 2*norm.sf(abs(z))
    expected_by = comp_b*bx
    out.append(dict(locus=loc,variant=v,role=role,outcome=oc,beta_LTL=bx,se_LTL=sx,beta_outcome=by,se_outcome=sy,
        wald_ratio=r,wald_se_2nd=se2,wald_lo=r-1.96*se2,wald_hi=r+1.96*se2,ratio_to_genomewide=r/gw_b,
        ivw_genomewide=gw_b,ivw_genomewide_se=gw_se,comparator=comp,comparator_b=comp_b,comparator_se=comp_se,
        expected_outcome_beta_under_comparator=expected_by,observed_over_expected=by/expected_by,z_diff=z,p_diff=p))
D=pd.DataFrame(out); D.to_csv('locus_wald_ratios_sampling_only.csv',index=False)
pd.set_option('display.width',250); print(D[['locus','variant','outcome','wald_ratio','wald_lo','wald_hi','ivw_genomewide','comparator_b','observed_over_expected','z_diff','p_diff']].round(4).to_string())
