"""Step 4 (part): coloc.abf between LTL and BP traits at TERC (3q26.2) and TERT (5p15.33), +/-500 kb (GRCh37).
Data: CMD Knowledge Portal ancestry-specific meta-analyses (EU; SA). SE derived from P (see Step 2 note)."""
import pandas as pd, numpy as np
from scipy.stats import norm
from mr_lib import coloc_abf
r=pd.read_csv('kp_regional_TERC_TERT.csv')
r['z']=norm.isf(r.p.clip(lower=1e-300)/2)*np.sign(r.beta)
r=r[(r.beta!=0)&r.z.notna()&(r.z!=0)]
r['se_p']=np.abs(r.beta/r.z)
out=[];top=[]
for loc in ['TERC','TERT']:
    for anc,traits in [('EU',['SBP','DBP','HYPERTENSION']),('SA',['SBP','DBP'])]:
        L=r[(r.locus==loc)&(r.trait=='LTL')&(r.anc==anc)].set_index('varId')
        for tr in traits:
            B=r[(r.locus==loc)&(r.trait==tr)&(r.anc==anc)].set_index('varId')
            m=L.join(B,lsuffix='_l',rsuffix='_b',how='inner')
            if len(m)<50: continue
            sd2=0.2 if tr=='HYPERTENSION' else 0.15
            res,pp=coloc_abf(m.beta_l.values,m.se_p_l.values,m.beta_b.values,m.se_p_b.values,sd1=0.15,sd2=sd2)
            i=np.argmax(pp)
            res.update(locus=loc,ancestry=anc,trait2=tr,top_snp_H4=m.rsid_l.iloc[i] if isinstance(m.rsid_l.iloc[i],str) else m.index[i],
                       top_snp_pp=pp[i],min_p_LTL=m.p_l.min(),min_p_trait=m.p_b.min(),
                       lead_LTL=m.rsid_l[m.p_l.idxmin()] ,lead_trait=m.rsid_b[m.p_b.idxmin()] if isinstance(m.rsid_b[m.p_b.idxmin()],str) else m.p_b.idxmin())
            out.append(res)
res=pd.DataFrame(out)[['locus','ancestry','trait2','nsnps','PP0','PP1','PP2','PP3','PP4','top_snp_H4','top_snp_pp','lead_LTL','min_p_LTL','lead_trait','min_p_trait']]
res.to_csv('step4_coloc_results.csv',index=False); pd.set_option('display.width',250); print(res.round(4).to_string())

# ---- Sensitivity: windows centred on the LTL lead (+/-250 kb, +/-100 kb, +/-50 kb)
r[['chr','pos']]=r.varId.str.split(':',expand=True).iloc[:,:2]; r['pos']=r.pos.astype(int)
sens=[]
for loc in ['TERC','TERT']:
    for anc,traits in [('EU',['SBP','DBP','HYPERTENSION']),('SA',['SBP','DBP'])]:
        L0=r[(r.locus==loc)&(r.trait=='LTL')&(r.anc==anc)]
        cpos=L0.loc[L0.p.idxmin(),'pos']
        for w in [250_000,100_000,50_000]:
            L=L0[(L0.pos-cpos).abs()<=w].set_index('varId')
            for tr in traits:
                B=r[(r.locus==loc)&(r.trait==tr)&(r.anc==anc)&((r.pos-cpos).abs()<=w)].set_index('varId')
                m=L.join(B,lsuffix='_l',rsuffix='_b',how='inner')
                sd2=0.2 if tr=='HYPERTENSION' else 0.15
                res,pp=coloc_abf(m.beta_l.values,m.se_p_l.values,m.beta_b.values,m.se_p_b.values,sd1=0.15,sd2=sd2)
                jb=m.p_b.idxmin()
                res.update(locus=loc,ancestry=anc,trait2=tr,window_kb=w//1000,centre=int(cpos),
                           lead_trait=m.rsid_b[jb],lead_trait_pos=int(m.pos_b[jb]),min_p_trait=m.p_b.min())
                sens.append(res)
S=pd.DataFrame(sens)[['locus','ancestry','trait2','window_kb','nsnps','PP3','PP4','lead_trait','lead_trait_pos','min_p_trait','centre']]
S.to_csv('step4_coloc_window_sensitivity.csv',index=False); print(S.round(4).to_string())
