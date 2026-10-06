"""Step 2: cross-ancestry look-up of rs10936599 (TERC/MYNN) and rs2736100 (TERT).
Source: CMD Knowledge Portal BioIndex API (bioindex.hugeamp.org), accessed 2026-09-29.
Effect allele = KP 'alt' allele (T for rs10936599, A for rs2736100) = LTL-shortening allele
(confirmed: Codd 2021 LTL beta < 0 for both). BP betas in KP standardised (SD) units; HTN betas log-OR.
"""
import pandas as pd, numpy as np
from scipy.stats import norm
import json
a=pd.read_csv('kp_snp_lookup_by_ancestry.csv'); d=pd.read_csv('kp_snp_lookup_by_dataset.csv')
# consistency check of KP ancestry records: z from beta/se vs z from p
a['z_se']=a.beta_alt/a.se; a['z_p']=norm.isf(a.p.clip(lower=1e-300)/2)*np.sign(a.beta_alt)
a['consistent']=(np.abs(np.abs(a.z_se)-np.abs(a.z_p))/np.maximum(np.abs(a.z_p),0.5)<0.35)|(a.p<1e-200)
a['se_p']=np.abs(a.beta_alt)/np.abs(a.z_p)
a['qc_pass']=a.se_p>=0.5*a.se  # implausibly precise p-implied SE -> record internally inconsistent
print(a[~a.consistent][['rs','ancestry','trait','beta_alt','se','p','z_se','z_p']])
d['se']=np.abs(d.beta_alt)/norm.isf(d.p.clip(lower=1e-300)/2)
d.to_csv('kp_snp_lookup_by_dataset.csv',index=False)
rsmap={'3:169492101:C:T':'rs10936599','5:1286516:C:A':'rs2736100'}; d['rs']=d.varId_b37.map(rsmap)

# Cross-ancestry pooling (non-overlapping ancestry groups), fixed + DerSimonian-Laird random
def pool(b,s):
    w=1/s**2; bf=np.sum(w*b)/w.sum(); sf=np.sqrt(1/w.sum()); Q=np.sum(w*(b-bf)**2); k=len(b)
    tau2=max(0,(Q-(k-1))/(w.sum()-np.sum(w**2)/w.sum())) if k>1 else 0
    wr=1/(s**2+tau2); br=np.sum(wr*b)/wr.sum(); sr=np.sqrt(1/wr.sum())
    from scipy.stats import chi2
    return dict(k=k,b_fixed=bf,se_fixed=sf,p_fixed=2*norm.sf(abs(bf/sf)),b_random=br,se_random=sr,p_random=2*norm.sf(abs(br/sr)),
                Q=Q,Q_p=chi2.sf(Q,k-1) if k>1 else np.nan,I2=max(0,(Q-(k-1))/Q)*100 if Q>0 else 0,tau2=tau2)
rows=[]
for rs in ['rs10936599','rs2736100']:
    for tr in ['SBP','DBP','HYPERTENSION']:
        x=a[(a.rs==rs)&(a.trait==tr)&a.ancestry.isin(['EU','SA','EA','AF','HS','AA'])&a.qc_pass]
        r=pool(x.beta_alt.values,x.se_p.values); r.update(rs=rs,trait=tr,ancestries=';'.join(x.ancestry)); rows.append(r)
pooled=pd.DataFrame(rows); pooled.to_csv('step2_cross_ancestry_pooled.csv',index=False)
print(pooled.round(5).to_string())
a.to_csv('kp_snp_lookup_by_ancestry_checked.csv',index=False)
