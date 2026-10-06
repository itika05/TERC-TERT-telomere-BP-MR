"""Step 3 (reverse direction): genetically predicted BP / hypertension -> LTL.
Exposure: MVP European (Verma 2024) SBP, DBP, hypertension; instruments = KP European BP/HTN lead variants
  with MVP P < 5e-8 (selection in the exposure sample itself), distance-clumped (strongest MVP P kept per +/-1 Mb;
  LD reference unavailable at this scale - stated limitation), MHC (chr6:25-34 Mb) excluded.
Outcome: Codd 2021 LTL (UK Biobank; no overlap with MVP). Secondary: Nakao 2026 South Asian LTL.
SE from P throughout."""
import pandas as pd, numpy as np
from scipy.stats import norm
from mr_lib import ivw, egger, weighted_median, weighted_mode, mr_presso, fstat, steiger_r2_cont, steiger_direction
d=pd.read_csv('kp_bp_leads_reverse_wide.csv').merge(pd.read_csv('kp_bp_leads_rsid.csv'),on='varId',how='left')
d=d[d.varId.str.count(':')==3].copy()
d[['chr','pos','ref','alt']]=d.varId.str.split(':',expand=True); d['pos']=d.pos.astype(int)
d=d[~((d.chr=='6')&d.pos.between(25_000_000,34_000_000))]
d=d[d.rsid.notna()&(np.minimum(d.maf,1-d.maf)>=0.01)]
LDt=pd.read_csv('ld_eur_pairs_v2.csv'); LDd={}
for a_,b_,r_ in LDt[['a','b','r2']].itertuples(index=False): LDd[(a_,b_)]=r_; LDd[(b_,a_)]=r_
LEDGER=[]
se=lambda b,p: np.abs(b)/norm.isf(np.clip(p,1e-300,1)/2)
res=[];snp=[]
for ex,exname in [('S_mvp','SBP'),('D_mvp','DBP'),('H_mvp','Hypertension')]:
    x=d[(d[ex+'_p']<5e-8)].sort_values(ex+'_p')
    kept=[]
    for r in x.itertuples():
        if any((k.chr==r.chr) and abs(k.pos-r.pos)<1_000_000 for k in kept): continue
        kept.append(r)
    n_dist=len(kept); kept2=[]
    for r in kept:   # kept is ordered by significance: drop any variant in LD (r2>=0.05, 1000G EUR) with a stronger retained one within 10 Mb
        if any(k.chr==r.chr and abs(k.pos-r.pos)<1e7 and (LDd.get((k.rsid,r.rsid),np.nan)>=0.05) for k in kept2): continue
        kept2.append(r)
    LEDGER.append(dict(exposure=exname,candidates=len(x),after_1Mb=n_dist,after_LD=len(kept2)))
    x=pd.DataFrame(kept2)
    for oc,ocname,nout in [('L_codd','LTL (Codd 2021, UKB)',464716),('L_nkSA','LTL South Asian (Nakao 2026)',11277)]:
        m=x.dropna(subset=[oc+'_b']).copy(); m=m[(m[oc+'_b']!=0)&(m[ex+'_b']!=0)]
        bx=m[ex+'_b'].values; sx=se(bx,m[ex+'_p'].values); by=m[oc+'_b'].values; sy=se(by,m[oc+'_p'].values)
        ok=np.isfinite(sx)&np.isfinite(sy)&(sx>0)&(sy>0)
        bx,sx,by,sy=bx[ok],sx[ok],by[ok],sy[ok]; m=m[ok]
        F=fstat(bx,sx)
        ests=[ivw(bx,sx,by,sy),egger(bx,sx,by,sy),weighted_median(bx,sx,by,sy),weighted_mode(bx,sx,by,sy),mr_presso(bx,sx,by,sy)]
        if exname!='Hypertension':
            nx=m[ex+'_n'].median()
            z,_=steiger_direction(steiger_r2_cont(bx,sx,nx),nx,steiger_r2_cont(by,sy,nout),nout); keep=z>0
            e=ivw(bx[keep],sx[keep],by[keep],sy[keep]); e['method']='IVW, Steiger-filtered'; ests.append(e)
            prop=keep.mean()
        else: prop=np.nan
        for e in ests:
            r={k:v for k,v in e.items() if k!='outlier_idx'}; r.update(exposure=exname,outcome=ocname,meanF=F.mean(),minF=F.min(),steiger_prop_correct=prop); res.append(r)
        for i in range(len(bx)): snp.append(dict(exposure=exname,outcome=ocname,varId=m.varId.iloc[i],bx=bx[i],sx=sx[i],by=by[i],sy=sy[i]))
pd.DataFrame(LEDGER).to_csv('step3_reverse_instrument_ledger.csv',index=False); print(pd.DataFrame(LEDGER))
R=pd.DataFrame(res); R.to_csv('step3_mr_reverse_results.csv',index=False); pd.DataFrame(snp).to_csv('step3_mr_reverse_snp_level.csv',index=False)
pd.set_option('display.width',250)
print(R[['exposure','outcome','method','nsnp','b','se','p','Q_p','I2','intercept','intercept_p','global_p','n_outliers','meanF','steiger_prop_correct']].round(4).to_string())
