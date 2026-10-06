"""Prior-sensitivity grid for coloc.abf (p12 in 1e-6..5e-5) across windows; also exports per-SNP posterior (H4) for plotting."""
import pandas as pd, numpy as np
from scipy.stats import norm
from mr_lib import coloc_abf
r=pd.read_csv('kp_regional_TERC_TERT.csv')
r['z']=norm.isf(r.p.clip(lower=1e-300)/2)*np.sign(r.beta); r=r[(r.beta!=0)&(r.z!=0)]; r['se_p']=np.abs(r.beta/r.z)
r['pos']=r.varId.str.split(':').str[1].astype(int)
out=[];snp=[]
for loc in ['TERC','TERT']:
    L0=r[(r.locus==loc)&(r.trait=='LTL')&(r.anc=='EU')]; c=L0.loc[L0.p.idxmin(),'pos']
    for tr in ['SBP','DBP','HYPERTENSION']:
        for w in [500,250,100,50]:
            L=L0[(L0.pos-c).abs()<=w*1000].set_index('varId'); B=r[(r.locus==loc)&(r.trait==tr)&(r.anc=='EU')&((r.pos-c).abs()<=w*1000)].set_index('varId')
            m=L.join(B,lsuffix='_l',rsuffix='_b',how='inner')
            for p12 in [1e-6,5e-6,1e-5,5e-5]:
                res,pp=coloc_abf(m.beta_l.values,m.se_p_l.values,m.beta_b.values,m.se_p_b.values,sd1=0.15,sd2=0.2 if tr=='HYPERTENSION' else 0.15,p12=p12)
                res.update(locus=loc,trait=tr,window_kb=w,p12=p12); out.append(res)
                if p12==1e-5 and w==100:
                    for v,x in zip(m.index,pp): snp.append(dict(locus=loc,trait=tr,varId=v,rsid=m.loc[v,'rsid_l'],pos=m.loc[v,'pos_l'],snp_pp4=x))
O=pd.DataFrame(out); O.to_csv('step4_coloc_prior_grid.csv',index=False); pd.DataFrame(snp).to_csv('step4_coloc_snp_pp4.csv',index=False)
pd.set_option('display.width',200)
print(O.pivot_table(index=['locus','trait','window_kb'],columns='p12',values='PP4').round(3))
print(O[O.p12==1e-5].pivot_table(index=['locus','trait'],columns='window_kb',values='PP3').round(6))
